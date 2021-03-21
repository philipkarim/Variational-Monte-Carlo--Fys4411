#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../particle.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeDoubleDerivativeNumeric (std::vector<Particle*> particles){
    //Computing the double derivative numeric, by evaluating the wavefunction 
    //using stepsize h

    //Defining some variables to be used in the calculations
    double h2, h, wfnext, wfprev, wf;
    double derivate2=0;

    //Step size
    h = 1e-4;
    h2=h*h;

    //Extracting the evaluation of the wave function in the current step 
    wf = m_system->getWaveFunction()->evaluate(particles);

    //Looping over the particles while looping over the dimensions of each particle
    //Using the adjust position function to evaluate forward and backward
    for (int i = 0; i<m_system->getNumberOfParticles(); i++){
        for (int dim = 0; dim<m_system->getNumberOfDimensions(); dim++){
            // Next position
            particles[i]->adjustPosition(h, dim);
            wfnext = m_system->getWaveFunction()->evaluate(particles);
            // Previous position
            particles[i]->adjustPosition(-2*h, dim);
            wfprev = m_system->getWaveFunction()->evaluate(particles);
            // The position now
            particles[i]->adjustPosition(h, dim);

            derivate2 -= (wfnext+wfprev-2*wf);
        }
    }
    //Returning the value of teh double derivative
    return 0.5*derivate2/h2;
    }
