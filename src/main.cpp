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

using namespace std;


int main() {

    //Init MPI
/*
    int numprocs, my_rank;
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
*/
    // Seed for the random number generator

    int seed = 2021;

    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) pow(2,16);  //16 or 17 would be nice
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z direction
    double alpha            = 0.5;          // Variational parameter.
    double timeStep         = 0.25;         // Metropolis time step (Importance sampling)
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool check_step         = true;
    bool numeric            = false;        // True->Numeric differentiation, False->Analytic
    bool bruteforce_val     = true;         // True->bruteforce, False->Importance sampling
    bool interaction        = false;
    bool GD                 = false;
    double initialAlpha     = 0.3;          //Initial alpha to start the gradient decent
    bool onebodydensity     =false;         //Extracting the positions to be used on the one body density
    //Writing to file
    bool GDwtf             =false;           //GD-Write to file
    bool generalwtf        =false;           //General information- write to file
    bool obdwtf            =false;           //One body density write to file

    double beta, a_length;                   //Defined under
    bool spherical;

    //Just making it easier to switch between interacting and non interacting cases
    if (interaction==true){
      a_length=0;
      //a_length=0.0043;                  //Trap length
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


    //Parallelisation
    //int my_rank, numprocs, idum;
    /*
    MPI_Init (NULL,NULL);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    */

    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z, beta));  //Added alpha
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setNumeric                  (numeric);
    system->setBruteforce               (bruteforce_val);
    system->setInteraction              (interaction);
    system->setTraplength               (a_length);
    system->setGD                       (GD);
    system->setGDwtf                    (GDwtf);
    system->setgeneralwtf               (generalwtf);
    
    if(onebodydensity==true){
      double bucketSize = 0.01;
      int bins = int(ceil(4 / bucketSize));

      system->setobd                    (obdwtf, bucketSize, bins);
    }

    if (GD==true){
      //MPI_Init(NULL,NULL);
      system->setTimeStep                 (timeStep);
      system->setStepLength               (stepLength);
      alpha = system->gradientDescent(initialAlpha);
      vector<double> parameters(2);
      parameters[0] = alpha;
      parameters[1] = beta;
      system->getWaveFunction()->setParameters(parameters);
      system->runMetropolisSteps           (numberOfSteps);
      //MPI_Finalize();
    }
    else if(GD==false && check_step==true){
      system->checkStep(stepLength, timeStep);
      }
    else{
      /*
      MPI_Init (NULL,NULL);
      MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
      */
      system->setTimeStep                 (timeStep);
      system->setStepLength               (stepLength);
      system->runMetropolisSteps          (numberOfSteps);
      //MPI_Finalize();
    }
  

    //MPI_Finalize ();

    return 0;
}
