#include <iostream>
#include <vector>
#include <fstream>

using namespace std;


int main()
{
    // Problem Parameters
    
    const double A=0.1; //Cross-sectional Area, m^2
    const double L=5;  //Length, m
    const double Ta=100; //Temperature at left face, °C
    const double Tb=200;  //Temperature at right face, °C
    const double S=1000;  //Source term, W/m3
    const double thermalCond=100; //Thermal conductivityü W/mK
    const int iteration=50; //Number of iterations
    
    const int N=5;     //Number of nodes
    const double d=L/N; //Cell spacing, m
    const double V=d*A; //Cell volume, m^3
    double coeffSum;
    vector<double> T(N,0); //Temperature at each node, °C
    vector<vector<double>> coeffMatrix(N,vector<double>(N,0));
    vector<double> k(N+1); //Thermal conductivity, W/m°K
    vector<double> coordFaces(N+1);
    vector<double> coordCentroids(N+2);
    vector<double> spaceBtwnCentroids(N+1);
    vector<double> diffusiveFlux(N+1);
    vector<double> aP(N);
    vector<vector<double>> aN;
    vector<double> Su(N);
    coordFaces[0]=0;
    coordCentroids[0]=-d/2;

    for (int i=0;i<N+1;i++)
    {
        coordCentroids[i+1]=coordCentroids[i]+d;
    }


    for (int i=0;i<N+1;i++)
    {
        spaceBtwnCentroids[i]=coordCentroids[i+1]-coordCentroids[i];
        k[i]=thermalCond;
        diffusiveFlux[i]=k[i]/spaceBtwnCentroids[i];
        coordFaces[i+1]=d*(i+1);
    }





    aN.resize(N,vector<double>(2));

    

    aP[0]=2*diffusiveFlux[0]*A+diffusiveFlux[1]*A;
    aN[0][0]=0;
    aN[0][1]=diffusiveFlux[1]*A;
    Su[0]=Ta*(2*diffusiveFlux[0]*A)+S*V;

    aP[N-1]=diffusiveFlux[N-1]*A+2*diffusiveFlux[N]*A;
    aN[N-1][0]=diffusiveFlux[N-1]*A;
    aN[N-1][1]=0;
    Su[N-1]=Tb*(2*diffusiveFlux[N]*A)+S*V;

    for (int i=1;i<N-1;i++) //Coefficients for interior nodes
    {
        aP[i]=diffusiveFlux[i]*A+diffusiveFlux[i+1]*A;
        aN[i][0]=diffusiveFlux[i]*A;
        aN[i][1]=diffusiveFlux[i+1]*A;
        Su[i]=S*V;
    }
    // for (int i = 0; i < N; i++)
    // {
    //     cout<<aN[i][1]<<endl;
    // }
    for (int i=0;i<N;i++)
    {
        coeffMatrix[i][i-1]=-aN[i][0];
        coeffMatrix[i][i]=aP[i];
        coeffMatrix[i][i+1]=-aN[i][1];
    }


    for (int k=0;k<iteration;k++)
    {
        
        for (int i=0;i<N;i++)
        {
            coeffSum=0;
            for (int j=0;j<i;j++)
            {
                coeffSum+=(-coeffMatrix[i][j]/coeffMatrix[i][i])*T[j];

            }
            for (int j = i+1; j < N; j++)
            {
                coeffSum+=(-coeffMatrix[i][j]/coeffMatrix[i][i])*T[j];
            }
            
            T[i]=coeffSum+Su[i]/coeffMatrix[i][i];
        }
        
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout<<coeffMatrix[i][j] << " ";
        }
        cout<<endl;
         
    }

    for (int i = 0; i < N; i++)
    {
        cout<<Su[i]<<endl;
    }
    
    for (int i = 0; i < N; i++)
    {
        cout<<T[i]<< " ";
    }
    
    
    







}