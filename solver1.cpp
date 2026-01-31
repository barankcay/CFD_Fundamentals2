#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

int main()
{
    // ==================== PROBLEM PARAMETERS ====================
    // Physical and geometric parameters for 1D steady-state heat conduction problem
    
    const double A = 0.1;              // Cross-sectional area, m²
    const double L = 5;                // Total length of the domain, m
    const double Ta = 100;             // Temperature at left boundary (x=0), °C
    const double Tb = 200;             // Temperature at right boundary (x=L), °C
    const double S = 1000;             // Volumetric heat source term, W/m³
    const double thermalCond = 100;    // Thermal conductivity, W/(m·K)
    int iteration = 0;                 // Iteration counter
    const double residualLimit = 1e-3; // Convergence criterion for maximum residual
    double maxResidual;                // Maximum residual value in current iteration
    
    // ==================== DISCRETIZATION PARAMETERS ====================
    const int N = 5;                   // Number of control volumes (nodes)
    const double d = L/N;              // Cell spacing (distance between centroids), m
    const double V = d*A;              // Volume of each control volume, m³
    double coeffSum;                   // Temporary sum for coefficient calculations
    double residualCoeffSum;           // Temporary sum for residual calculations

    // ==================== SOLVER PARAMETERS ====================
    const double alpha = 1.3;          // Successive Over-Relaxation (SOR) coefficient (>1 for over-relaxation)
    
    // ==================== VARIABLE DECLARATIONS ====================
    vector<double> T(N, 0);                           // Temperature at each node centroid, °C (initialized to 0)
    vector<vector<double>> coeffMatrix(N, vector<double>(N, 0)); // Coefficient matrix [A] for system [A]{T}={b}
    vector<double> k(N+1);                            // Thermal conductivity at each face, W/(m·K)
    vector<double> coordFaces(N+1);                   // x-coordinates of cell faces
    vector<double> coordCentroids(N+2);               // x-coordinates of cell centroids (includes ghost cells)
    vector<double> spaceBtwnCentroids(N+1);           // Distance between adjacent centroids, m
    vector<double> diffusiveFlux(N+1);                // Diffusive flux coefficient (k/Δx) at each face, W/(m²·K)
    vector<double> aP(N);                             // Central coefficient (aP) for each node
    vector<vector<double>> aN;                        // Neighbor coefficients: aN[i][0]=aW (west), aN[i][1]=aE (east)
    vector<double> Su(N);                             // Source term vector for each node, W
    vector<double> residual(N);                       // Residual vector for convergence check
    
    // ==================== GEOMETRY INITIALIZATION ====================
    coordFaces[0] = 0;                 // First face at x = 0
    coordCentroids[0] = -d/2;          // Ghost centroid to the left (for boundary condition)

    // Calculate centroid coordinates
    for (int i = 0; i < N+1; i++)
    {
        coordCentroids[i+1] = coordCentroids[i] + d;
    }

    // Calculate face properties and diffusive flux coefficients
    for (int i = 0; i < N+1; i++)
    {
        spaceBtwnCentroids[i] = coordCentroids[i+1] - coordCentroids[i];  // Distance between centroids
        k[i] = thermalCond;                                                // Thermal conductivity at face i
        diffusiveFlux[i] = k[i] / spaceBtwnCentroids[i];                  // Γ/δx where Γ=k and δx=distance
        coordFaces[i+1] = d * (i+1);                                       // Face coordinates
    }

    // ==================== COEFFICIENT MATRIX SETUP ====================
    aN.resize(N, vector<double>(2));   // Resize neighbor coefficient array

    // --- Left Boundary Node (i=0) ---
    // Special treatment for Dirichlet boundary condition at left face
    aP[0] = 2*diffusiveFlux[0]*A + diffusiveFlux[1]*A;  // aP = aW + aE (aW uses factor of 2 for boundary)
    aN[0][0] = 0;                                        // No west neighbor (boundary)
    aN[0][1] = diffusiveFlux[1]*A;                       // aE (east neighbor coefficient)
    Su[0] = Ta*(2*diffusiveFlux[0]*A) + S*V;             // Source includes boundary condition contribution

    // --- Right Boundary Node (i=N-1) ---
    // Special treatment for Dirichlet boundary condition at right face
    aP[N-1] = diffusiveFlux[N-1]*A + 2*diffusiveFlux[N]*A;  // aP = aW + aE (aE uses factor of 2 for boundary)
    aN[N-1][0] = diffusiveFlux[N-1]*A;                      // aW (west neighbor coefficient)
    aN[N-1][1] = 0;                                          // No east neighbor (boundary)
    Su[N-1] = Tb*(2*diffusiveFlux[N]*A) + S*V;              // Source includes boundary condition contribution

    // --- Interior Nodes (i=1 to N-2) ---
    // Standard discretization for internal control volumes
    for (int i = 1; i < N-1; i++)
    {
        aP[i] = diffusiveFlux[i]*A + diffusiveFlux[i+1]*A;  // aP = aW + aE
        aN[i][0] = diffusiveFlux[i]*A;                       // aW (west neighbor)
        aN[i][1] = diffusiveFlux[i+1]*A;                     // aE (east neighbor)
        Su[i] = S*V;                                         // Source term (volumetric generation only)
    }
    
    // Populate the coefficient matrix [A] from aP and aN values
    // Matrix form: aP*T[i] - aW*T[i-1] - aE*T[i+1] = Su[i]
    for (int i = 0; i < N; i++)
    {
        coeffMatrix[i][i-1] = -aN[i][0];  // -aW (west coefficient, negative in matrix)
        coeffMatrix[i][i] = aP[i];         // aP (diagonal/central coefficient)
        coeffMatrix[i][i+1] = -aN[i][1];  // -aE (east coefficient, negative in matrix)
    }

    // ==================== ITERATIVE SOLVER (SOR Method) ====================
    maxResidual = 10;  // Initialize with value larger than residualLimit
    
    while (maxResidual > residualLimit)
    {
        // --- Update Temperature Field ---
        // Gauss-Seidel with Successive Over-Relaxation (SOR)
        for (int i = 0; i < N; i++)
        {
            coeffSum = 0;
            
            // Sum contributions from west neighbors (already updated in this iteration)
            for (int j = 0; j < i; j++)
            {
                coeffSum += (-coeffMatrix[i][j] / coeffMatrix[i][i]) * T[j];
            }
            
            // Sum contributions from east neighbors (old values from previous iteration)
            for (int j = i+1; j < N; j++)
            {
                coeffSum += (-coeffMatrix[i][j] / coeffMatrix[i][i]) * T[j];
            }
            
            // SOR update formula: T_new = (1-α)*T_old + α*T_GaussSeidel
            T[i] = (1-alpha)*T[i] + alpha*(coeffSum + Su[i]/coeffMatrix[i][i]); 
        }

        // --- Calculate Residuals ---
        // Residual = Su - [A]*{T} for each equation
        for (int i = 0; i < N; i++)
        {
            residualCoeffSum = 0;
            for (int j = 0; j < N; j++)
            {
                residualCoeffSum += coeffMatrix[i][j] * T[j];  // Matrix-vector multiplication
            }
            residual[i] = Su[i] - residualCoeffSum;  // Calculate residual for node i
        }
        
        iteration = iteration + 1;  // Increment iteration counter
        maxResidual = *max_element(residual.begin(), residual.end());  // Find maximum residual
    }

    // ==================== OUTPUT RESULTS ====================
    cout << "Solution took " << iteration << " iterations to converge" << endl;
    cout << "Maximum residual value is " << maxResidual << endl; 
    
    // Print temperature distribution
    for (int i = 0; i < N; i++)
    {
        cout << T[i] << " ";
    }
}