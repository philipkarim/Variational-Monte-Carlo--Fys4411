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

    double r2, r_tot, derivate2, potentialEnergyPart=0;
    double x_part, y_part, z_part;
    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    double wf=m_system->getWaveFunction()->evaluate(particles);
    std::vector<double> r_pos;

    //Computing the kinetic energy
    kineticEnergy=m_system->getWaveFunction()->computeDoubleDerivative(particles);

    //Computing the potential energy (V_ext)
    std::vector<int> dimensions_length(m_system->getNumberOfDimensions());
    std::iota(dimensions_length.begin(), dimensions_length.end(), 0);

    for(int i=0; i<m_system->getNumberOfParticles(); i++){
      r_pos=particles[i]->getPosition();

      if (dimensions_length.size()==3){
        x_part=r_pos[0]*r_pos[0];
        y_part=r_pos[1]*r_pos[1];
        z_part=r_pos[2]*r_pos[2];
      }
      else if(dimensions_length.size()==2){
        x_part=r_pos[0]*r_pos[0];
        y_part=r_pos[1]*r_pos[1];
        z_part=0;
      }
      else {
        x_part=r_pos[0]*r_pos[0];
        y_part=0;
        z_part=0;
      }

      potentialEnergyPart+=(x_part+y_part+z_part);
    }
    potentialEnergy=0.5*potentialEnergyPart*m_omega*m_omega*wf;

    cout<<m_system->getNumeric();

    return (kineticEnergy + potentialEnergy)/wf;

}
