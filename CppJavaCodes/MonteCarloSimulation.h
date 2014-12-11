//
//  MonteCarloSimulation.h
//
//
//  Created by Omid Saremi, last edit Sept 2014. MIT License
//
//
#include <iostream>
#include <armadillo>
#include <fstream>
#include <math.h>
#include <ctime>
#include <complex>
#include <string>
#include <assert.h>

using namespace arma;
using namespace std;

class MonteCarloSimulation
{
    
public:
    
    MonteCarloSimulation(string initializationType, string boundaryType, int mylatticeSize, int myRank, double myRange, double myCoupling){
        
        assert((initializationType=="Hot" || initializationType=="Cold") && (boundaryType=="Dirichlet" || boundaryType=="Open"));
        
        // Theory parameters
        
        N=myRank;
        coupling=myCoupling;
        
        // Simulation parameters
        
        latticeSize=mylatticeSize;
        range=myRange;
        
        // Setting the boundary condition
        
        if(boundaryType=="Open"){
            startIndex=0;
            endIndex=getLatticeSize();
        }
        else{
            startIndex=1;
            endIndex=getLatticeSize()-1;
        }
        
        cx_mat iden=eye<cx_mat>(getGaugeGroup(), getGaugeGroup());
        
        // Create the configuration matrix
        
        u=cx_cube(getGaugeGroup(), getGaugeGroup(), getLatticeSize());
        
        // Boundary conditions at the end of the interval
        
        if(boundaryType=="Dirichlet"){
            u.slice(0)=iden;
            u.slice(getLatticeSize()-1)=iden;
        }
        
        // Initilization starts here. "Hot" initialization, configuration matrix is initialized to to be random
        
        if(initializationType=="Hot"){
            for(int i=getStartIndex(); i< getEndIndex(); ++i){
                u.slice(i)=generateUnitaryRandomMatrix();
            }
        } // "Cold". configuration matrix is set to be identity
        else{
            for(int i=getStartIndex(); i< getEndIndex();++i){
                u.slice(i)=iden;
            }
        }
    }
    
    cx_mat generateUnitaryRandomMatrix(){
        
        int n=getGaugeGroup();
        cx_mat randMatrix=randn<cx_mat>(n, n);
        cx_mat Hermitian=0.5*(randMatrix + trans(randMatrix))*getRange();
        cx_mat iden=eye<cx_mat>(n, n);
        
        cx_mat accum=1.0j*Hermitian;
        cx_mat unitaryGroupMember=iden+accum;
        
        for(int i=2; i < CUTOFF; ++i)
        {
            accum=accum*1.0j*Hermitian/i;
            unitaryGroupMember+=accum;
        }
        
        Extra step to improve unitarity using QR decomposition
        
        cx_mat Q(n, n);
        cx_mat R(n, n);
        qr(Q, R, unitaryGroupMember);
        return Q;
        
        return unitaryGroupMember;
    }
    
    double findDistanceToIdentity(cx_mat input){
        
        int numbRows=input.n_rows;
        cx_mat iden=eye<cx_mat>( numbRows, numbRows );
        
        return  accu(abs(input-iden))/pow(numbRows, 2.0);
    }
    
    bool isUnitary(cx_mat input){
        
        int numbRows=input.n_rows;
        cx_mat id1=input*trans(input);
        cx_mat id2=trans(input)*input;
        
        if(abs(accu(id1)-(double ) numbRows)< ERR && abs(accu(id2)-(double ) numbRows)< ERR ) return true;
        
        return false;
    }
    
    bool isUnitary(){
        
        for(int i=getStartIndex(); i< getEndIndex(); ++i){
            if(!isUnitary(u.slice(i))) return false;
        }
        
        return true;
    }
    
    void unitarizeConfiguration(){
        
        int n=getGaugeGroup();
        for(int i=getStartIndex(); i< getEndIndex(); ++i ){
            if(!isUnitary(u.slice(i))){
                cx_mat Q(n, n);
                cx_mat R(n, n);
                qr(Q, R, u.slice(i));
                u.slice(i)=Q;
            }
        }
    }
    
    double computeEnergy(){
        
        double sum=0.0;
        for(int i=getStartIndex(); i< getEndIndex()-1; ++i){
            sum +=2.0*real(trace(trans(u.slice(i+1))*u.slice(i)));
        }
        
        return sum/(double)(getEndIndex()-getStartIndex()-1)/getGaugeGroup();     // The observable is calculated in units of N and length
    }
    
    double computeNumberDensity(){
        
        double sum=0.0;
        for(int i=getStartIndex();i< getEndIndex(); ++i){
            sum +=pow(abs(trace(u.slice(i))),2.0);
        }
        
        return sum/(double)(getEndIndex()-getStartIndex())/pow(getGaugeGroup(), 2.0);     // The observable is calculated in units of N^2 and length
    }
    
    double action(){    // Action required for reweighting
        
        double sum=0.0;
        for(int i=getStartIndex();i< getEndIndex(); ++i){
            sum +=pow(abs(trace(u.slice(i))),2.0);
        }
        
        return sum;    // The observable is calculated in units of N^2 and length
    }
    
    double totalAction(){
        double a=0.1;
        double sum=0.0;
        for(int i=getStartIndex();i< getEndIndex()-1; ++i){
            sum +=-a*pow(abs(trace(u.slice(i))),2.0)-2.0*real(trace(trans(u.slice(i+1))*u.slice(i)))/a;
;
        }
        
        return sum-a*pow(abs(trace(u.slice(getEndIndex()-1))), 2.0);
        
    }
    
    double computeL1(){
        
        complex<double> sum=0.0;
        for(int i=getStartIndex();i< getEndIndex(); ++i){
            sum +=trace(u.slice(i));
        }
        
        return abs(sum)/(double)(getEndIndex()-getStartIndex())/getGaugeGroup();
    }
    
    
    double computeL2(){
        
        complex<double> sum=0.0;
        for(int i=getStartIndex();i< getEndIndex(); ++i){
            sum +=trace(u.slice(i));
        }
        
        return pow(abs(sum),2.0)/(double)pow(getEndIndex()-getStartIndex(), 2.0)/pow(getGaugeGroup(), 2.0);
    }
    
    
    double getVariationInAction(int latticeSite, cx_mat pickedXU){
        
        double a=0.1;
        double deltaS;
        // What is needed is (X-1)U
        cx_mat Y=pickedXU-u.slice(latticeSite);
        
        if(latticeSite==0){
            return 2.0*getCoupling()*real(trace(trans(u.slice(1))*Y));
        }
        
        else if(latticeSite==getLatticeSize()-1){
            return 2.0*getCoupling()*real(trace(trans(Y)*u.slice(latticeSite-1)));
        }
        
        else{
            
            deltaS=-((2.0/a)*real(accu(conj(Y) % (u.slice(latticeSite-1)+u.slice(latticeSite+1))))+a*(pow(abs(trace(pickedXU)),2.0)-pow(abs(trace(u.slice(latticeSite))),2.0)));
        }
        
        return deltaS;
    }
    
    void MetropolisHastingLocalUpdate(){
        double beta=getCoupling();
        int index;
        cx_mat pickedXU;
        double deltaS;
        for(int i=getStartIndex(); i<getEndIndex();++i){
            index=rand() % SIZE_X;
            pickedXU=X.slice(index)*u.slice(i);
            deltaS=getVariationInAction(i, pickedXU);
        
            if(deltaS < 0){
                
                u.slice(i)=pickedXU;
            }
            else if(exp(-beta*deltaS) >= (double) rand()/RAND_MAX){    
                u.slice(i)=pickedXU;
            }
            else{
                rejection++;
            }
        }
    }
    
    
    void generateSetX(){
        int n=getGaugeGroup();
        X=cx_cube(n, n, SIZE_X);
        for(int i=0; i< SIZE_X/2; ++i){
            X.slice(i)=generateUnitaryRandomMatrix();
            X.slice(i+SIZE_X/2)=trans(X.slice(i)); // Set half of the Xs to their Hermitian conjugate to respect "Detailed Balance Principle"
        }
    }
    
    
    void setLatticeSize(int mySize){
        
        latticeSize=mySize;
    }
    
    int getLatticeSize(){
        
        return latticeSize;
    }
    
    void setGaugeGroup(int rankOfGaugeGroup){
        
        N=rankOfGaugeGroup;
    }
    
    int getGaugeGroup(){
        
        return N;
    }
    
    void setRange(double rangeOfRandomMatrix){
        
        range=rangeOfRandomMatrix;
    }
    
    double getRange(){
        
        return range;
    }
    
    cx_cube getConfiguration(){
        
        return u;
    }
    
    void setCoupling(double myCoupling){
        
        coupling=myCoupling;
        
    }
    
    double getCoupling(){
        
        return coupling;
    }
    
    void resetRejection(){
        
        rejection=0;
    }
    
    int getRejection(){
        
        return rejection;
    }
    
    int getStartIndex(){
        return startIndex;
    }
    
    int getEndIndex(){
        return endIndex;
    }
    
    double static const ERR=1e-10;     // Tolerance in unitarity test
    int static const CUTOFF=35;
    int static const SIZE_X=50000;
    int static const EQUILIBRIUM_STEPS=100000;
    
private:
    
    int latticeSize;
    double coupling;
    double range;
    int N;
    int startIndex;
    int endIndex;
    
    cx_cube u;
    cx_cube X;
    int rejection;
};
