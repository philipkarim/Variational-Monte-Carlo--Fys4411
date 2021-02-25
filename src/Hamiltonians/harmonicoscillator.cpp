#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"


using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, double omega_z) :
        Hamiltonian(system) {
    assert(omega   > 0);
    assert(omega_z > 0);
    m_omega  = omega;
    m_omega_z= omega_z;
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

    double potentialEnergyPart=0;
    double x_part, y_part, z_part;
    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    double wf=m_system->getWaveFunction()->evaluate(particles);
    std::vector<double> r_pos_CLE;

//    kineticEnergy=m_system->getWaveFunction()->computeDoubleDerivative(particles);


    //Computing the kinetic energy
    if(m_system->getNumeric()==false){
      kineticEnergy=m_system->getWaveFunction()->computeDoubleDerivative(particles);
    }
    else if(m_system->getNumeric()==true){
      kineticEnergy=m_system->getHamiltonian()->computeDoubleDerivativeNumeric(particles);
    }

    //Computing the potential energy (V_ext)
    std::vector<int> dimensions_length_CLE(m_system->getNumberOfDimensions());
    std::iota(dimensions_length_CLE.begin(), dimensions_length_CLE.end(), 0);

/*
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
      r_pos_CLE=particles[i]->getPosition();

      if (dimensions_length_CLE.size()==3){
        x_part=r_pos_CLE[0]*r_pos_CLE[0];
        y_part=r_pos_CLE[1]*r_pos_CLE[1];
        z_part=r_pos_CLE[2]*r_pos_CLE[2];
      }
      else if(dimensions_length_CLE.size()==2){
        x_part=r_pos_CLE[0]*r_pos_CLE[0];
        y_part=r_pos_CLE[1]*r_pos_CLE[1];
        z_part=0;
      }
      else {
        x_part=r_pos_CLE[0]*r_pos_CLE[0];
        y_part=0;
        z_part=0;
      }

      potentialEnergyPart+=(x_part+y_part+z_part);
    }
    potentialEnergy=0.5*potentialEnergyPart*m_omega*m_omega;
*/
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
      for(int dim=0; dim<m_system->getNumberOfDimensions(); dim++){
        if (dim==2){
          potentialEnergyPart += m_omega_z*m_omega_z*particles[i]->getPosition()[dim]*particles[i]->getPosition()[dim];
        }
        else{
        potentialEnergyPart += m_omega*m_omega*particles[i]->getPosition()[dim]*particles[i]->getPosition()[dim];
        }
      }
    }

/*
      r_pos_CLE=particles[i]->getPosition();

      if (dimensions_length_CLE.size()==3){
        x_part=r_pos_CLE[0]*r_pos_CLE[0];
        y_part=r_pos_CLE[1]*r_pos_CLE[1];
        z_part=r_pos_CLE[2]*r_pos_CLE[2];
      }
      else if(dimensions_length_CLE.size()==2){
        x_part=r_pos_CLE[0]*r_pos_CLE[0];
        y_part=r_pos_CLE[1]*r_pos_CLE[1];
        z_part=0;
      }
      else {
        x_part=r_pos_CLE[0]*r_pos_CLE[0];
        y_part=0;
        z_part=0;
      }
      potentialEnergyPart+=(x_part+y_part+z_part);
    }
    */
    potentialEnergy=0.5*potentialEnergyPart;


    return (kineticEnergy/wf + potentialEnergy);

}
