#include <iostream>
#include "utils.h"

using namespace std;

//#define Globalvariables søk på det
//Constants, most of them is set equal to one, for a dimensionless calculation
double h_bar=1.0;
double a, beta, alpha; //Constants to be defined later
double a_ho, r; //Diameter of bosons
double m=1, w_ho=1, w_z=1;

//Choose if the trap is speherical (S) or elliptical (E)
char trap='S';
//Choose if we have a non-interaction case (0) or an interaction (1)
int interaction=0;

int main() {

  // Seed for the random number generator

  int n_dim                 = 1;          //Number of dimensions
  int n_particle            = 1;          //Number of particles
  int n_steps         = 1000000;          //Number of steps
  double omega            = 1.0;          // Oscillator frequency.
  double alpha            = 0.5;          // Variational parameter.
  double step_length      = 0.1;          // Metropolis step length.
  double equilibration    = 0.1;          // Amount of the total steps used
  // for equilibration.

  System* system = new System(int 2021);
  system->setHamiltonian              (new HarmonicOscillator(system, omega));
  system->setWaveFunction             (new SimpleGaussian(system, alpha));
  system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
  system->setEquilibrationFraction    (equilibration);
  system->setStepLength               (stepLength);
  system->runMetropolisSteps          (numberOfSteps);

  return 0;
}
