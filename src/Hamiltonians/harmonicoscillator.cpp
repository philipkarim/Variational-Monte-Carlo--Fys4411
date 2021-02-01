#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"


using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */
/*
    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    //double wf=m_system->getWaveFunction()->evaluate(particles);

    //Computing the kinetic energy
    kineticEnergy=(1.0/m_system->getWaveFunction()->evaluate(particles))*
                  0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);

    //Computing the potential energy
    for (int i = 0; i<m_system->getNumberOfParticles();i++){
           for(int dim = 0; dim<m_system->getNumberOfDimensions(); dim++){
               double rr = particles.at(i)->getPosition().at(dim);
               potentialEnergy+= 0.5*m_omega*m_omega*rr*rr;
           }
       }

    return kineticEnergy + potentialEnergy;
}
*/
    double r2;
    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    double wf=m_system->getWaveFunction()->evaluate(particles);
    //double wf=m_system->getWaveFunction()->evaluate(particles);

    //Computing the kinetic energy
    for (int i = 0; i<m_system->getNumberOfParticles();i++){
        kineticEnergy-=0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);

    }

    //Computing the potential energy (V_ext)
    for (int i = 0; i<m_system->getNumberOfParticles();i++){
           for(int dim = 0; dim<m_system->getNumberOfDimensions(); dim++){
               double r2 = particles.at(i)->getPosition().at(dim);
               potentialEnergy+= r2*r2;
           }
    }
//
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
      r_pos=particles[i]->getPosition();
      for (int dim=0; dim<m_system->getNumberOfDimensions(); dim++){
        r_tot+=r_pos[dim]*r_pos[dim];
        derivate2-=2*m_parameters[0]*(2*m_parameters[0]*r_tot)*wf;
      }
      derivate2+=2*m_parameters[0]*m_system->getNumberOfDimensions()*wf;
//

    potentialEnergy=0.5*m_omega*m_omega*r2;



    return (kineticEnergy + potentialEnergy)/wf;
}
