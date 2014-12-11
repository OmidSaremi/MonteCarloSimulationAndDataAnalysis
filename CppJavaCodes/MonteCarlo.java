import java.util.Random;
import java.lang.Math;
import java.util.Arrays;
/* Omid Saremi, MIT License 
 * Monte carlo Markov chain  (Metropolis-Hasting local update algorithm) implemented in Java, 
 * simulates a simple quantual mechanical path integral.   
 */

public class MonteCarlo{
    
    private int latticeSize;           // Number of sites in the latticized system 
    private double hoppingParam;       // Hopping parameter of the model
    private double range;              // Metropolis update step's range: -range < delta <= range
    private double[] u;                // Configuration vector
    private static final int MEASURED_EVERY=5;                   // Measure correlators every "MEASURED_EVERY" Monte Carlo sweeps to reduce autocorrelations 
    private static final int NUM_SWEEPS_TO_EQUILIBRIATE=500;     // To sample the distribution, a number of sweeps are needed, where no measurement is performed
    private int reject=0;                                        // Rejection count in the Metropolis local update step
    
    public MonteCarlo(int latticeSize, double hoppingParam){
        this.latticeSize=latticeSize;
        this.hoppingParam=hoppingParam;
    }
    
    public void initialize(){     // Note that u[0]=u[latticeSize-1]=0 (fixed boundary points)
        u=new double[latticeSize];   
        Random rand=new Random();
        for(int i=1;i<latticeSize-1; i++){
            u[i]=range*(rand.nextDouble()-0.5); 
        }
     return;
    }
   
    public void MetropolisLocalUpdateSweep(){
        double delta, deltaS;
        Random rand=new Random();
        for(int i=0; i< MEASURED_EVERY;i++){
            for(int j=1;j<latticeSize-1; j++){
               delta=range*(rand.nextDouble()-0.5); 
               deltaS=Math.pow(delta,2.0)+delta*(2.0*u[j]-hoppingParam*(u[j+1]+u[j-1]));
               if(deltaS < 0.0){ 
                   u[j]+=delta;
               }
               else if(Math.exp(-deltaS) >= rand.nextDouble()){
                    u[j]+=delta;
               }
               else{
                   reject++;   
                  }
            }
        }
        }
        
    public void setRange(double range){
        this.range=range;
    }
    
    public void resetRejectionCount(){
     this.reject=0;
    }
    public double findCorrelator(int power){
        double accum=0.0;
        for(int i=1; i<latticeSize-1; i++){
            accum=accum+Math.pow(u[i], power);
        }
        return accum/(latticeSize-2);
    }

    public static void main(String[] args){
        
        int latticeSize=500;
        int numMeasure=100000;
        double cor2=0.0;
        double coupling=(1.0-1.0/4.0*(0.25)*(0.25))/(1.0+1.0/4.0*(0.25)*(0.25));
        System.out.println("Hopping parameter: "+coupling);
        double acceptanceRate;
        
        MonteCarlo MC=new MonteCarlo(latticeSize, coupling);
        MC.setRange(4.2);   // Parameter range here controls the acceptance rate. The value which gives 50% reject rate is optimal
        MC.initialize();    // Perform "Hot initialization" of the configuration array 
        System.out.println(Arrays.toString(MC.u));
        
        //"Thermalize" the system to ensure, we are sampling the desired distribution function  
        
        for(int i=0; i< NUM_SWEEPS_TO_EQUILIBRIATE; i++){
        MC.MetropolisLocalUpdateSweep();
        }
        
        MC.resetRejectionCount(); // We don't care about rejection count up to this point. We are still in pre-measurement period
        
        for(int i=0;i<numMeasure; i++){
         MC.MetropolisLocalUpdateSweep(); // Note that "MEASURED_EVERY" sweeps are discarded in between measurements to have 
                                          //independent configurations sampled from the the distribution function
         cor2=cor2+MC.findCorrelator(2);  //Measuring the second moment which is related to the ground state energy of the system
        }
        
        cor2=cor2/numMeasure;
        
        System.out.println("Ground state energy is: "+ cor2);
        double exactEnergy=2.0*(1.0+1.0/64.0);
        System.out.println("Exact value: "+ exactEnergy);
        System.out.println("Rejection rate was: "+ (float) MC.reject/(latticeSize-2)/numMeasure/MEASURED_EVERY);
    }
}
    
    
