#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include <unistd.h>

using namespace std;

int main() {

    // Seed for the random number generator
    int seed = 2021;

    int numberOfDimensions  = 1;            // Set amount of dimensions
    int numberOfParticles   = 1;            // Set amount of particles
    int numberOfSteps       = (int) pow(2,19); //Amount of metropolis steps
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z direction
    double alpha            = 0.5;          // Variational parameter.
    double timeStep         = 0.25;         // Metropolis time step (Importance sampling)
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool check_step         = false;        // Simulating as function of step sizes
    bool numeric            = true;         // True->Numeric differentiation, False->Analytic
    bool bruteforce_val     = true;         // True->bruteforce, False->Importance sampling
    bool interaction        = false;        // True-> interaction, False->Not interaction
    bool GD                 = false;        // True-> gradient decent, False->Not gradietn decent
    double initialAlpha     = 0.45;         // Initial alpha to start the gradient decent
    bool collectresults     =false;         // True-> aquiring large amount of results in parallel
    bool onebodydensity     =false;         // Extracting the positions to be used on the one body density
    //Write to file
    bool GDwtf             =false;          // GD-Write to file
    bool generalwtf        =false;          // General information- write to file
    bool obdwtf            =false;          // One body density write to file

    double beta, a_length;                  // Defined under
    bool spherical;

    //Just making it easier to switch between interacting and non interacting cases
    if (interaction==true){
      //a_length=0;                     // Without the Jastrow factor remove //
      a_length=0.0043;                  // Trap length
      beta=2.82843;                     // Beta value
      spherical=false;                  // Trap symmetry
    }
    else{
      a_length         =0.0;            // Trap length
      beta             =1.0;            // Beta value
      spherical=true;                   // Trap symmetry
    }
    if (spherical==true){
      omega_z=omega;
    }

    //Setting the different values defined higher in the code
    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z, beta));
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
    
    //One body density values
    if(onebodydensity==true){
      double bucketSize = 0.01;
      int bins = int(ceil(400));
      system->setobd                    (obdwtf, bucketSize, bins);
    }
    //Gradient decent method
    if (GD==true){
      if (collectresults==true){
        int pid, pid1, pid2;
        //Parallelizing the code using 4 cores
        pid=fork();       if (pid==0){system->setInitialState(new RandomUniform(system, 3, 50));
                              alpha = system->gradientDescent(0.45);}

        else{pid1=fork();if (pid1==0){system->setInitialState(new RandomUniform(system, 3, 10));
                              alpha = system->gradientDescent(0.7);}
        
        else{pid2=fork();if (pid2==0){system->setInitialState(new RandomUniform(system, 3, 100));
                              alpha = system->gradientDescent(0.3);}

        else                         {system->setInitialState(new RandomUniform(system, 3, 100));
                              alpha = system->gradientDescent(0.45);}
        }}}
      //Running a regular gradient decent without parallel computing
      else{
      alpha = system->gradientDescent(initialAlpha);
      vector<double> parameters(2);
      parameters[0] = alpha;
      parameters[1] = beta;
      system->getWaveFunction()->setParameters(parameters);
      //system->runMetropolisSteps           (numberOfSteps);
      } 
    }

    //Looking for the best step sizes by running the script in
    //parallel for each steplength and timestep
    else if(GD==false && check_step==true){
      int pid, pid1, pid2;

      pid=fork();       if (pid==0){system->checkStep(1, 0.25);}
      else{pid1=fork();if (pid1==0){system->checkStep(0.5, 0.1);}
      else{pid2=fork();if (pid2==0){system->checkStep(0.75, 0.05);}
      else                         {system->checkStep(0.05, 0.005);}
      }}}

    //Running multiple benchmarking simulations in parallel, ensure that the computer has 6 cores in this case
    else if(collectresults==true && GD==false && check_step==false){
      int pid, pid1, pid2, pid3, pid4;

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
      //Else just run a regular simulation without GD or parallizisation
      system->runMetropolisSteps          (numberOfSteps);
    }
  

    return 0;
}
