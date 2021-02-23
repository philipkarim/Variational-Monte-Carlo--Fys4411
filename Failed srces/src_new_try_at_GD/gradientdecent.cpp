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

using namespace std;


void runGradientDecent(std::vector<double> input_variables_double, std::vector<int> input_variables_int, std::vector<bool> input_variables_bool){
    double alpha_curr=input_variables_double[0];
    double alpha_prev=alpha_curr;  //Initial guess
    double alpha_next;
    double El_deriv2, El_deriv_prev;
    double gamma_dir;
    double epsilon=1E-6;

    int iterations=1;   //Not sure how many iterations to use, 50 or 500??

    for(int ii=0; ii<iterations; ii++){
        El_deriv2=E_LDerivative(alpha_curr, input_variables_double, input_variables_int, input_variables_bool);
        El_deriv_prev=E_LDerivative(alpha_prev, input_variables_double, input_variables_int, input_variables_bool);

        //Division by 0, added the 5
        gamma_dir=(alpha_curr-alpha_prev+epsilon)/(5);//El_deriv2-El_deriv_prev);

        alpha_next=alpha_curr-gamma_dir*(El_deriv2);
        alpha_prev=alpha_curr;
        alpha_curr=alpha_next;
    }

}

    double E_LDerivative(double alpha_n, std::vector<double> in_double, std::vector<int> in_int, std::vector<bool> in_bool){

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
        grad_values.reserve(3);

        //Think this might cause the segmentation fault
        //m_sampler                   =Sampler();
        //Sampler testing;
        Sampler* sample2=nullptr;

        //grad_values=sample2->getGradientDecentValues();
        grad_values=sample2->return_grad();
        //Sampler obj;
        //grad_values =obj.getGradientDecentValues();

        grad_values[0]=1;
        grad_values[1]=1;
        grad_values[2]=1;


        expectEnerg=grad_values[0];
        expectE_L_deri=grad_values[1];
        expectderi_dot_EL=grad_values[2];
        double expectEnergy    = expectEnerg/1;//(system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());
        double expectE_L_deriv     = expectE_L_deri/1;//(system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());
        double expectderiv_dot_EL = expectderi_dot_EL/1;//(system->getNumberOfMetropolisSteps());//*getEquilibrationFraction());

        double E_L_derive=2*(expectderiv_dot_EL - expectEnergy*expectE_L_deriv)/1;//system->getEquilibrationFraction();

        return  E_L_derive;
    }
