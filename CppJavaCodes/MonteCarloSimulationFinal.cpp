//
//  MonteCarloSimulationFinal.cpp
//
// Created orginally by Omid Saremi Sept 2012. Copyright Omid Saremi, MIT license
// Last Edit by Omid Saremi on Sept 2014
//
//Build on my laptop: g++ MonteCarloSimulationFinal.cpp  -o MonteCarloSimulationFinal -O3 -march=native -larmadillo
//Build on the Lawrence Berkeley National Lab cluster : g++ -I/share/apps/armadillo/include -L/share/apps/armadillo/lib MonteCarloSimulationFinal.cpp -o MonteCarloSimulationFinal -O3 -march=native -larmadillo

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

int main(int argc, char* argv[])
{
    // Simulation parameters passed on to main from the command line
    
    int N=atoi(argv[1]);
    double systemSize=atof(argv[2]);
    int mySize=atoi(argv[3]);
    double lambda=atof(argv[4]);
    int numG=atoi(argv[5]);
    cout<< "Entered integer is " << numG << endl;
    char* outputOrder=argv[6];
    char* outputAcceptance=argv[7];
    
    // Other simulation parameters
    
    clock_t t;
    string initType="Hot";
    string boundaryType="Dirichlet";
    double a=systemSize/mySize;
    double myRange;         
    double myHopping, gamma, myCoupling;
    int numMeasurements=200000;
    int numDiscard=300;
   
    // Create the appropriate number of output files.
    
    ofstream outdataOrder, outdataAcceptance; 
    outdataOrder.open(outputOrder, ios_base :: app);
    outdataAcceptance.open(outputAcceptance, ios_base :: app);
    
    // Simulation starts here. Initialize the random generator
    
    srand(time(NULL));
    gamma=4.6;
    myHopping =(lambda-1.0)*0.7/(numG-1)+0.9;
    myCoupling =pow(N*myHopping/2.0/gamma, 0.5);
    
    // Find the range using Newton-Raphson Algorithm
    
    // Initializing the Newton's iteration
    
    myRange=0.02;
    double  acceptOld=0.0;
    double acceptNew=0.0;
    double diffAcceptance=0.0;
    double deltaRange=0.0;
    int testRun=300;
    int counter=0;
    
    t=clock();
    while ( abs(acceptOld-0.5) > 0.03){
        
        MonteCarloSimulation MC1=MonteCarloSimulation(initType, boundaryType, mySize, N, myRange, myCoupling);
        MC1.generateSetX();
        MC1.resetRejection();
        
        for(int i=0; i< testRun; ++i){
            MC1.MetropolisHastingLocalUpdate();
        }
        
        acceptOld=1.0-1.0*MC1.getRejection()/(MC1.getLatticeSize()-2.0)/testRun;
        
        cout<< acceptOld << endl;
        
        MonteCarloSimulation MC2=MonteCarloSimulation(initType, boundaryType, mySize, N, myRange+0.001, myCoupling);
        MC2.generateSetX();
        MC2.resetRejection();
        
        for(int i=0; i< testRun; ++i){
            MC2.MetropolisHastingLocalUpdate();
        }
        
        acceptNew=1.0-1.0*MC2.getRejection()/(MC2.getLatticeSize()-2.0)/testRun;
        
        diffAcceptance=(acceptNew-acceptOld)/0.001;
        deltaRange=(0.5-acceptOld)/diffAcceptance;
        myRange=myRange+deltaRange;
        counter++;
    }
    
    cout<< "Range: " << myRange<< "  " << "Acceptance rate: " << 100.0*acceptOld << " " << "# of iterations: " << counter << endl;
    
    // Initilize the Monte Carlo Simulation here
    
    MonteCarloSimulation MC=MonteCarloSimulation(initType, boundaryType, mySize, N, myRange, myCoupling);
    MC.generateSetX();
    
    // MonteCarlo sweep until equilibrium is reached
    
    for(int i=0; i< MC.EQUILIBRIUM_STEPS; ++i){
        
        MC.MetropolisHastingLocalUpdate();
    }
    
    double numDensity=0.0; 
    MC.resetRejection();
    for(int i=0; i < numMeasurements; ++i ){
        for(int j=0; j < numDiscard; ++j){  // To reduce autocorrelation, discard intermediate configurations
            MC.MetropolisHastingLocalUpdate();
        }
        numDensity=MC.computeNumberDensity();
        outdataOrder << numDensity << " ";
    }
    
    // Add new line and the parametr value to the already open files
    
    outdataOrder << myHopping << endl;
    outdataAcceptance << 100.0-100.0*MC.getRejection()/(double) (MC.getLatticeSize()-2.0)/numDiscard/numMeasurements  << " " << myHopping << endl;
    
    
    outdataOrder.close();
    outdataAcceptance.close();
    
    t=clock()-t;
    cout<< " Time elapsed was "<< t/(double) CLOCKS_PER_SEC<< endl;
    
    return 0;
}
