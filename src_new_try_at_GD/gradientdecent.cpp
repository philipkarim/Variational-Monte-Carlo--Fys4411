#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
//#include "mpi.h"
#include "gradientdecent.h"
#include "sampler.h"


void runGradientDecent(std::vector<double> input_variables_double, std::vector<int> input_variables_int, std::vector<bool> input_variables_bool){

    //Remember to put this in main and setting the alpha
    //before running the main metropolis. Run for about 1000 steps
    double alpha_curr=input_variables_double[0];
    double alpha_prev=alpha_curr;  //Initial guess

    double alpha_next;
    double El_deriv2, El_deriv_prev;
    double gamma_dir;
    double epsilon=1E-6;

    int iterations=10;   //Not sure how many iterations to use, 50 or 500??

    //double beta = getWaveFunction()->getParameters()[2]/getWaveFunction()->getParameters()[0];

    for(int ii=0; ii<iterations; ii++){
        //std::vector<double> parameters(2);
        //parameters[0] = alpha_curr;
        //parameters[2] = alpha_curr*beta;
        //Need to find a way to set the parameters public for each itteration
        //system->setWaveFunction             (new SimpleGaussian(system, alpha_curr, beta));

        El_deriv2=E_LDerivative(alpha_curr, input_variables_double, input_variables_int, input_variables_bool);
        El_deriv_prev=E_LDerivative(alpha_prev, input_variables_double, input_variables_int, input_variables_bool);

        gamma_dir=(alpha_curr-alpha_prev+epsilon)/(El_deriv2-El_deriv_prev);

        alpha_next=alpha_curr-gamma_dir*(El_deriv2);

        alpha_prev=alpha_curr;
        alpha_curr=alpha_next;
    }

}


    double E_LDerivative(double alpha_n, std::vector<double> in_double, std::vector<int> in_int, std::vector<bool> in_bool){
      /*
      extern double omega, beta, timeStep, stepLength, equilibration;
      extern int numberOfDimensions, numberOfParticles, seed;
      extern bool numeric, bruteforce_val;
*/
      //Based on the main function, might want to import the variables from main in a fancy way

        System* system = new System(in_int[2]);
        system->setHamiltonian              (new HarmonicOscillator(system, in_double[1]));  //Added alpha
        system->setWaveFunction             (new SimpleGaussian(system, alpha_n, in_double[2]));
        system->setInitialState             (new RandomUniform(system, in_int[0], in_int[1]));
        system->setEquilibrationFraction    (in_double[5]);
        system->setStepLength               (in_double[4]);
        system->setNumeric                  (in_bool[0]);
        system->setBruteforce               (in_bool[1]);
        system->setTimeStep                 (in_double[3]);

        system->runMetropolisSteps          (1000);

        double expectEnerg, expectE_L_deri, expectderi_dot_EL;
        std::vector<double> grad_values;
        grad_values.reserve(2);

        //Think this might cause the segmentation fault
        //m_sampler                   =Sampler();
        //Sampler testing;
        //Sampler* sample;

        //grad_values=sample->getGradientDecentValues();
        //Sampler obj;
        //grad_values =obj.getGradientDecentValues();

        grad_values[0]=1;
        grad_values[1]=1;
        grad_values[2]=1;


        expectEnerg=grad_values[0];
        expectE_L_deri=grad_values[1];
        expectderi_dot_EL=grad_values[2];
        double expectEnergy    = expectEnerg/(system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());
        double expectE_L_deriv     = expectE_L_deri/(system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());
        double expectderiv_dot_EL = expectderi_dot_EL/(system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());

        double E_L_derive=2*(expectderiv_dot_EL - expectEnergy*expectE_L_deriv)/system->getEquilibrationFraction();

        return  E_L_derive;
    }
