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

#include "gradientdecent.h"
#include "sampler.h"
#include <vector>
//#include "main.h"


GradientDecent::GradientDecent() {
}

void GradientDecent::runGradientDecent(double alpha_init, double omega, double beta,
  double timeStep, double stepLength, double equilibration, int numberOfDimensions, int numberOfParticles,
             int seed, bool numeric, bool bruteforce_val){

    //Remember to put this in main and setting the alpha
    //before running the main metropolis. Run for about 1000 steps
    double alpha_curr=alpha_init;
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

        El_deriv2=E_LDerivative(alpha_curr, omega, beta,
                     timeStep, stepLength, equilibration, numberOfDimensions, numberOfParticles,
                     seed, numeric, bruteforce_val);
        El_deriv_prev=E_LDerivative(alpha_prev, omega, beta,
                     timeStep, stepLength, equilibration, numberOfDimensions, numberOfParticles,
                     seed, numeric, bruteforce_val);

        gamma_dir=(alpha_curr-alpha_prev+epsilon)/(El_deriv2-El_deriv_prev);

        alpha_next=alpha_curr-gamma_dir*(El_deriv2);

        alpha_prev=alpha_curr;
        alpha_curr=alpha_next;
    }

}


    double GradientDecent::E_LDerivative(double alpha_n, double omega, double beta,
      double timeStep, double stepLength, double equilibration, int numberOfDimensions, int numberOfParticles,
                 int seed, bool numeric, bool bruteforce_val){
      /*
      extern double omega, beta, timeStep, stepLength, equilibration;
      extern int numberOfDimensions, numberOfParticles, seed;
      extern bool numeric, bruteforce_val;
*/
      //Based on the main function, might want to import the variables from main in a fancy way

      System* system = new System(seed);
      system->setHamiltonian              (new HarmonicOscillator(system, omega));  //Added alpha
      system->setWaveFunction             (new SimpleGaussian(system, alpha_n, beta));
      system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
      system->setEquilibrationFraction    (equilibration);
      system->setStepLength               (stepLength);
      system->setNumeric                  (numeric);
      system->setBruteforce               (bruteforce_val);
      system->setTimeStep                 (timeStep);

      system->runMetropolisSteps          (1000);

        double expectEnerg, expectE_L_deri, expectderi_dot_EL;
        std::vector<double> grad_values;
        grad_values.reserve(2);

        //Think this might cause the segmentation fault
        //m_sampler                   =Sampler();
        grad_values=getSampler()->getGradientDecentValues();
        //Sampler obj;
        //grad_values =obj.getGradientDecentValues();

        expectEnerg=grad_values[0];
        expectE_L_deri=grad_values[1];
        expectderi_dot_EL=grad_values[2];
        double expectEnergy    = expectEnerg/(m_system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());
        double expectE_L_deriv     = expectE_L_deri/(m_system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());
        double expectderiv_dot_EL = expectderi_dot_EL/(m_system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());

        double E_L_derive=2*(expectderiv_dot_EL - expectEnergy*expectE_L_deriv)/m_system->getEquilibrationFraction();
        return  E_L_derive;
    }
