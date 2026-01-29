#include <iostream>
#include <vector>
#include <fstream>

using namespace std;


int main()
{
    // Problem Parameters
    const double k=100; //Thermal Conductivity, W/mK
    const double A=0.01; //Cross-sectional Area, m^2
    const double L=0.5;  //Length, m
    const double Ta=100; //Temperature at left face, °C
    const double Tb=25;  //Temperature at right face, °C

    // Geometry Discretization
    const int N=10;     //Number of nodes
    const double d=L/N; //Cell spacing, m
    const double V=d*A; //Cell volume, m^3

    vector<double> coordFaces;
    vector<double> coordCentroids;
    vector<double> spaceBtwnCentroids;






}