#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../particle.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeDoubleDerivativeNumeric (std::vector<Particle*> particles){

double h, wfnext, wfprev, wf;
double derivate2=0;
double h2;

h = 1e-5;
h2=h*h;
//might need a smaller steplength

wf = m_system->getWaveFunction()->evaluate(particles);

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

        derivate2 -= (wfnext+wfprev-2*wf)/h2;
    }
}
return 0.5*derivate2;
}
