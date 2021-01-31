#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}
/*
//These three were removed with typos, not needed for brute force metro?
void System::setPosition(const std::vector<double> &position) {
    assert(position.size() == (unsigned int) m_numberOfDimensions);
    m_position = position;
}

void System::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
*/
//_______________

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

     //_________BRUTE FORCE_____________

     int random_index=0;
     double psi_factor;
     double wfold=m_stepLength;
     std::vector<double> PositionOld=std::vector<double>();

     //Random integer generator
     std::random_device rd;
     std::mt19937_64 gen(rd());
     std::uniform_int_distribution<int> distribution(0,m_numberOfParticles-1);
     std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);

     //Random index used to choose a random particle
     random_index=distribution(gen);

     //Defining the random particle:


    PositionOld=m_particles[random_index]->getPosition();

     //Choosing a random step:
     double temp_rand=UniformNumberGenerator(gen);
     double step=m_stepLength*(temp_rand-0.5);

     //Start the step which gives movement of the particle
     for (int dim=0; dim<m_numberOfDimensions; dim++){
       m_particles[random_index]->adjustPosition(step, dim);
     }

     //Extracting the new wavefunction, and checks if it is accepted
     double wfnew=m_waveFunction->evaluate(m_particles);
     psi_factor=wfnew*wfnew/(wfold*wfold);
     //If accepted:
     if (temp_rand<=psi_factor){
        wfold=wfnew;
        return true;
     }
     else{
         m_particles[random_index]->setPosition(PositionOld);
        return false;
      }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */

         //If statement to make send the accepted steps into the sampler
         //after the system is at rest
         if (i>=numberOfMetropolisSteps*m_equilibrationFraction){
            m_sampler->sample(acceptedStep);
         }
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}
