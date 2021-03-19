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
#include <unistd.h>

using namespace std;


int main() {

    // Seed for the random number generator

    int seed = 2021;

    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) pow(2,19);  //16 or 17 would be nice
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z direction
    double alpha            = 0.5;          // Variational parameter.
    double timeStep         = 0.25;         // Metropolis time step (Importance sampling)
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool check_step         = false;
    bool numeric            = true;         // True->Numeric differentiation, False->Analytic
    bool bruteforce_val     = true;         // True->bruteforce, False->Importance sampling
    bool interaction        = false;
    bool GD                 = false;
    double initialAlpha     = 0.3;          //Initial alpha to start the gradient decent
    bool collectresults     =true;           //True-> aquiring large amount of results in parallel
    bool onebodydensity     =false;         //Extracting the positions to be used on the one body density
    //Writing to file
    bool GDwtf             =false;           //GD-Write to file
    bool generalwtf        =true;          //General information- write to file
    bool obdwtf            =false;          //One body density write to file

    double beta, a_length;                   //Defined under
    bool spherical;

    //Just making it easier to switch between interacting and non interacting cases
    if (interaction==true){
      //a_length=0;
      a_length=0.0043;                  //Trap length
      beta=2.82843;                     //Beta value
      spherical=false;                  //Trap symmetry
    }
    else{
      a_length         =0.0;            //Trap length
      beta             =1.0;            //Beta value
      spherical=true;                   //Trap symmetry
    }
    if (spherical==true){
      omega_z=omega;
    }


    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z, beta));  //Added alpha
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setTimeStep                 (timeStep);
    system->setStepLength               (stepLength);
    system->setEquilibrationFraction    (equilibration);
    system->setNumeric                  (numeric);
    system->setBruteforce               (bruteforce_val);
    system->setInteraction              (interaction);
    system->setTraplength               (a_length);
    system->setGD                       (GD);
    system->setGDwtf                    (GDwtf);
    system->setgeneralwtf               (generalwtf);
    
    //One body density plot
    if(onebodydensity==true){
      double bucketSize = 0.01;
      int bins = int(ceil(4 / bucketSize));

      system->setobd                    (obdwtf, bucketSize, bins);
    }
    //Gradient decent method
    if (GD==true){
      if (collectresults==true){
        int pid, pid1, pid2, pid3, pid4, pid5, pid6;

        pid=fork();       if (pid==0){system->setInitialState(new RandomUniform(system, 3, 10));
                              alpha = system->gradientDescent(0.3);}
        else{pid1=fork();if (pid1==0){system->setInitialState(new RandomUniform(system, 3, 10));
                              alpha = system->gradientDescent(0.7);}
        
        /*else{pid2=fork();if (pid2==0){system->setInitialState(new RandomUniform(system, 3, 100));
                              alpha = system->gradientDescent(0.3);}

        else{pid3=fork();if (pid3==0){system->setInitialState(new RandomUniform(system, 3, 10));
                                      system->setWaveFunction(new SimpleGaussian(system, 0.5, 2.82843));
                                      system->setTraplength(0.0043);
                                      system->setNumeric(true);
                                      interaction=true;
                                      spherical=false;
                              alpha = system->gradientDescent(0.4);}

        else{pid4=fork();if (pid4==0){system->setInitialState(new RandomUniform(system, 3, 50));
                                      system->setWaveFunction(new SimpleGaussian(system, 0.5, 2.82843));
                                      system->setTraplength(0.0043);
                                      system->setNumeric(true);
                                      interaction=true;
                                      spherical=false;
                              alpha = system->gradientDescent(0.4);}*/
                              
        /*else{pid5=fork();if (pid5==0){system->setInitialState(new RandomUniform(system, 3, 100));
                                      system->setWaveFunction(new SimpleGaussian(system, 0.5, 2.82843));
                                      system->setTraplength(0.0043);
                                      system->setNumeric(true);
                                      interaction=true;
                                      spherical=false;
                                      double bucketSize = 0.01;
                                      int bins = int(ceil(4 / bucketSize));
                                      system->setobd(obdwtf, bucketSize, bins);}
                                      
        else{pid6=fork();if (pid6==0){system->setInitialState(new RandomUniform(system, 3, 100));
                                      system->setWaveFunction(new SimpleGaussian(system, 0.5, 2.82843));
                                      system->setTraplength(0);
                                      system->setNumeric(true);
                                      interaction=true;
                                      spherical=false;
                                      double bucketSize = 0.01;
                                      int bins = int(ceil(4 / bucketSize));
                                      system->setobd(obdwtf, bucketSize, bins);}*/

        else                         {system->setInitialState(new RandomUniform(system, 3, 100));
                                      //system->setWaveFunction(new SimpleGaussian(system, 0.5, 2.82843));
                                      //system->setNumeric(true);
                                      //system->setTraplength(0.0043);
                                      //interaction=true;
                                      //spherical=false;
                              alpha = system->gradientDescent(0.3);}
        }}
      else{
      alpha = system->gradientDescent(initialAlpha);
      vector<double> parameters(2);
      parameters[0] = alpha;
      parameters[1] = beta;
      system->getWaveFunction()->setParameters(parameters);
      system->runMetropolisSteps           (numberOfSteps);
      } 
    }

    //Looking for the best step sizes by running the script in
    //parallel for each steplength and timestep
    else if(GD==false && check_step==true){
      int pid, pid1, pid2, pid3, pid4;

      pid=fork();       if (pid==0){system->checkStep(1, 0.25);}
      else{pid1=fork();if (pid1==0){system->checkStep(0.5, 0.1);}
      else{pid2=fork();if (pid2==0){system->checkStep(0.75, 0.05);}
      else{pid3=fork();if (pid3==0){system->checkStep(0.25, 0.01);}
      else{pid4=fork();if (pid4==0){system->checkStep(0.1, 0.5);}
      else                         {system->checkStep(0.05, 0.005);}
      }}}}}


    else if(collectresults==true && GD==false && check_step==false){
      int pid, pid1, pid2, pid3, pid4, pid5, pid6;

      pid=fork();       if (pid==0){system->setInitialState(new RandomUniform(system, 3, 500));}
      else{pid1=fork();if (pid1==0){system->setInitialState(new RandomUniform(system, 2, 500));}
      else{pid2=fork();if (pid2==0){system->setInitialState(new RandomUniform(system, 3, 500));}
      else{pid3=fork();if (pid3==0){system->setInitialState(new RandomUniform(system, 1, 500));
                                    system->setBruteforce                    (false);}
      else{pid4=fork();if (pid4==0){system->setInitialState(new RandomUniform(system, 2, 500));
                                    system->setBruteforce                    (false);}
      else{                         system->setInitialState(new RandomUniform(system, 3, 500));
                                    system->setBruteforce                    (false);}
      }}}}
      system->runMetropolisSteps          (numberOfSteps);
    }

    else{
      system->runMetropolisSteps          (numberOfSteps);
    }
  

    return 0;
}
