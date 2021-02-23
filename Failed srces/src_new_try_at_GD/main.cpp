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
//#include "gradientdecent.h"


using namespace std;


int main() {

    //Init MPI
    //int numprocs, my_rank;
    //MPI_Init (&nargs, &args);
    //MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    //MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // Seed for the random number generator
    int seed = 2020;

    //Dim=1, particle=1 should give 0.5
    int numberOfDimensions  = 3;
    int numberOfParticles   = 1;
    int numberOfSteps       = (int) 1e5;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double beta             = 1;
    double timeStep         = 0.01;         // Metropolis time step (Importance sampling)
    double stepLength       = 1;            // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    bool numeric            = true;         // True->Numeric differentiation, False->Analytic
    bool bruteforce_val     = true;         // True->bruteforce, False->Importance sampling
    bool GD                 = false;         // Calculate best parameter by GD, (true or false)
    double alpha_guess_GD   =0.45;


    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega));  //Added alpha
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setBeta                     (beta);
    system->setNumeric                  (numeric);
    system->setBruteforce               (bruteforce_val);
    system->setTimeStep                 (timeStep);
    if (GD==false){

    system->runMetropolisSteps          (numberOfSteps);

    }
    else{

    //system->runGradientDecent(alpha_guess_GD);

      /*
    //std::vector<std::variant<double, int, bool>> input_variables;

    std::vector<double> input_doubles;
    std::vector<int> input_ints;
    std::vector<bool> input_bools;

    //Defining vectors with the input
    input_doubles.push_back(alpha_guess_GD);
    input_doubles.push_back(omega);
    input_doubles.push_back(beta);
    input_doubles.push_back(timeStep);
    input_doubles.push_back(stepLength);
    input_doubles.push_back(equilibration);

    input_ints.push_back(numberOfDimensions);
    input_ints.push_back(numberOfParticles);
    input_ints.push_back(seed);

    input_bools.push_back(numeric);
    input_bools.push_back(bruteforce_val);


*/
    }

    //system->gradientDecent();

    //system->setNumberOfParticles        (numeric)

    //MPI_Finalize ();

    return 0;
}
