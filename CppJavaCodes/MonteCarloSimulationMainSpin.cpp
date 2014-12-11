//
//  MonteCarloSimulationMainSpin.cpp
//
// Created orginally by Omid Saremi Sept 2012. Copyright Omid Saremi, MIT license
// Last Edit by Omid Saremi on Sept 2014
//
//Build on my laptop: g++ MonteCarloSimulationMainSpin.cpp  -o MonteCarloSimulationMainSpin -O1 -larmadillo
//qstat -f -ne
//qstat -u omid.saremi

//Build on the Lawrence Berkeley National Lab cluster : g++ -I/share/apps/armadillo/include -L/share/apps/armadillo/lib MonteCarloSimulationMainSpin.cpp -o MonteCarloSimulationMainSpin -O1 -larmadillo
//qsub -b y -cwd myfile
//qsub -q fast2  -V -cwd -b y ./MonteCarloSimulationMainSpin
//

#include <iostream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <complex>
#include <string>
#include <cstring>
#include <sstream>
#include <armadillo>
#include "MonteCarloSimulation.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    clock_t t;
    
    // Simulation parameters
    
    string initType="Hot";
    string boundaryType="Dirichlet";
    double a=systemSize/mySize;
    double myRange=0.03;              
    double myCoupling=0.0;
    double myHopping=0.0;
    int numMeasurements=3000;        
    int numDiscard=20;
    int measurePerFile=10000000;
    
    // Create the appropriate number of output files. Up to 1e7 double precision numbers per file (at most 100 MB stored in each output file)
    
    int numOfFiles=ceil((double) numMeasurements/measurePerFile);
    ofstream outputStream[numOfFiles];
    ofstream poly;
    poly.open("./PolyakovU(N)SpinTheoryMeasurementResults.txt", ios_base::app);
    string Name="./U(N)SpinTheoryMeasurementResults"; 
    string outFullName[numOfFiles];
    
    for(int i=0; i< numOfFiles; ++i){
        stringstream outTemp;
        outTemp << i;
        outFullName[i]=Name+outTemp.str()+".txt";
        outputStream[i].open(outFullName[i].c_str(), ios_base::app);
        cout<< "File name: "<< outFullName[i]<< endl;
    }
    
    // Simulation starts here. Initialize the random generator
    
    srand(time(NULL));
    
    // Get my acceptance rates
    
    int numG=10;               //12
    for(int g=0; g < numG;++g){
        t=clock();
        double printCoupling=0.0*(1.0/2.0/((0.5)*g/(numG-1)+4.1)/a)+1.0/2.0/4.6/a;
        myCoupling =-1.0*N*printCoupling;
        myHopping =-3.0*(g*0.4/(numG-1)+0.5/3.0)*a;
        myRange=myRange/1.1;
        
        MonteCarloSimulation MC=MonteCarloSimulation(initType, boundaryType, mySize, N, myRange, myCoupling, myHopping);
        MC.generateSetX();
        
        // Monte Carlo sweep until equilibrium is reached
        
        for(int i=0; i< MC.EQUILIBRIUM_STEPS; ++i){
            
            MC.MetropolisHastingLocalUpdate();
        }
        
        cout<< "I am measuring now! "<< endl;
        
        double energyDensity,energyDensity2, averageEnergy, averageEnergy2; // Observable is reported per unit color degree of freedom squared namely N^2
        double energySum=0.0;
        double energySum2=0.0;
       
        MC.resetRejection();
        for(int i=0; i < numMeasurements; ++i ){
            
           for(int j=0; j < numDiscard; ++j){  // To reduce autocorrelation, discard intermediate configurations
                
            MC.MetropolisHastingLocalUpdate();
            }
            
            int fileNumber=i/measurePerFile;
            energyDensity=MC.computeL2();
            energyDensity2=MC.computeL1();
            energySum +=energyDensity;
            energySum2+=energyDensity2;
            outputStream[fileNumber] << energyDensity <<" ";
            if (i % 100==0) {
                
            cout<< "num of measurements so far: "<< i << endl;
            cout<< "Is unirary? "<< MC.isUnitary()<< endl;
            }
        }
        
        // Add new line and the parametr value to the already open files
        
        for(int i=0 ; i< numOfFiles; ++i){
            outputStream[i]<< myHopping/a <<endl;
        }
        
        // Compute the averages and write the value of individual measurements into the output files
        
        averageEnergy=energySum/(double) numMeasurements;
        averageEnergy2=energySum2/(double) numMeasurements;
        cout<< "N= "<< N << endl;
        cout<< "The coupling beta is= "<< myHopping/a << endl;
        cout<< "Average energy in units of N^2 is= "<< averageEnergy <<endl;
        cout<< "Acceptance rate is= "<< 100.0-100.0*MC.getRejection()/(double ) (MC.getLatticeSize()-2)/numDiscard/numMeasurements << endl;
        poly<< myHopping/a <<" "<< mySize*(averageEnergy-pow(averageEnergy2,2.0)) << endl;
        
        // Simulation ends. Print the total running time
        
        t=clock()-t;
        
        cout<< "time elapsed was "<< t/(double) CLOCKS_PER_SEC<< endl;
    }
    
    // Close the output files
    
    for(int i=0; i< numOfFiles;++i ){
        outputStream[i].close();
    }
    
    poly.close();
    
    return 0;
}
