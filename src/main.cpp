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

using namespace std;


int main() {
    // Seed for the random number generator
    int seed = 2020;

    //Dim=1, particle=1 should give 0.5
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e5;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double beta             = 1;
    double timeStep         = 0.01;         // Metropolis time step (Importance sampling)
    double stepLength       = 1;            // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    bool numeric            = false;         // True->Numeric differentiation, False->Analytic
    bool bruteforce_val     = false;         // True->bruteforce, False->Importance sampling

    System* system = new System(seed);
    system->setHamiltonian              (new HarmonicOscillator(system, omega));  //Added alpha
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setNumeric                  (numeric);
    system->setBruteforce               (bruteforce_val);
    system->setTimeStep                 (timeStep);

    system->runMetropolisSteps          (numberOfSteps);


    //system->setNumberOfParticles        (numeric)

    return 0;
}
