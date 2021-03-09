#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"


using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, double omega_z, double gamma) :
        Hamiltonian(system) {
    assert(omega   > 0);
    assert(omega_z > 0);
    m_omega  = omega;
    m_omega_z= omega_z;
    m_gamma  = gamma;

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

    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    double interactionEnergy=0;
    double wf=m_system->getWaveFunction()->evaluate(particles);

//    kineticEnergy=m_system->getWaveFunction()->computeDoubleDerivative(particles);

    if (m_system->getInteraction()==false){

    //Computing the kinetic energy
    if(m_system->getNumeric()==false){
      kineticEnergy=m_system->getWaveFunction()->computeDoubleDerivative(particles);
    }
    else if(m_system->getNumeric()==true){
      kineticEnergy=m_system->getHamiltonian()->computeDoubleDerivativeNumeric(particles);
    }

    //Computing the potential energy (V_ext)/trap
    std::vector<int> dimensions_length_CLE(m_system->getNumberOfDimensions());
    std::iota(dimensions_length_CLE.begin(), dimensions_length_CLE.end(), 0);

    potentialEnergy=computePotentialEnergy(particles);

    return (kineticEnergy/wf + potentialEnergy);
  }

  else{//Interacting case

    //Calculating the x,y, z parts
    potentialEnergy=computePotentialEnergyInteracting(particles);
    
    //Calculating the energy produced by interaction(V_int)
    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
        for(int j = i+1; j < m_system->getNumberOfParticles(); j++) {
            interactionEnergy+=computeEnergyInteracting(m_system->getParticles()[i],m_system->getParticles()[j]);
        }
    }

    //Calculating the double derivative part
    kineticEnergy=m_system->getHamiltonian()->computeDoubleDerivativeNumeric(particles);

    return (kineticEnergy/wf + potentialEnergy+interactionEnergy);

  }
}

//Potential energy noninteracting case
double HarmonicOscillator::computePotentialEnergy(std::vector<Particle*> particles) {
  double potentialEnergyPart=0;
  double x_part, y_part, z_part;
  std::vector<int> dimensions_length(m_system->getNumberOfDimensions());
  std::iota(dimensions_length.begin(), dimensions_length.end(), 0);
  std::vector<double> r_pos_P;
    
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
      r_pos_P=particles[i]->getPosition();

      if (dimensions_length.size()==3){
        x_part=r_pos_P[0]*r_pos_P[0];
        y_part=r_pos_P[1]*r_pos_P[1];
        z_part=r_pos_P[2]*r_pos_P[2];
      }
      else if(dimensions_length.size()==2){
        x_part=r_pos_P[0]*r_pos_P[0];
        y_part=r_pos_P[1]*r_pos_P[1];
        z_part=0;
      }
      else {
        x_part=r_pos_P[0]*r_pos_P[0];
        y_part=0;
        z_part=0;
      }

      potentialEnergyPart+=(x_part+y_part+z_part);
    }
    return 0.5*potentialEnergyPart*m_omega*m_omega;

}

//Potential energy interacting case(Technically just the xyz part og the Ep)
double HarmonicOscillator::computePotentialEnergyInteracting(std::vector<Particle*> particles) {
    double x_part, y_part, z_part;
    std::vector<int> dimensions_length(m_system->getNumberOfDimensions());
    std::iota(dimensions_length.begin(), dimensions_length.end(), 0);
    std::vector<double> r_pos_D;
    double PE=0;

   for(int i=0; i<m_system->getNumberOfParticles(); i++){
      r_pos_D=particles[i]->getPosition();

      if (dimensions_length.size()==3){
        x_part=r_pos_D[0]*r_pos_D[0];
        y_part=r_pos_D[1]*r_pos_D[1];
        z_part=m_gamma*m_gamma*r_pos_D[2]*r_pos_D[2];
      }
      else if(dimensions_length.size()==2){
        x_part=r_pos_D[0]*r_pos_D[0];
        y_part=r_pos_D[1]*r_pos_D[1];
        z_part=0;
      }
      else {
        x_part=r_pos_D[0]*r_pos_D[0];
        y_part=0;
        z_part=0;
      }

      PE+=x_part+y_part+z_part;
    }

    return 0.5*PE;

}

//Computing the energy produced by boson-boson interaction
double HarmonicOscillator::computeEnergyInteracting(Particle* particle1, Particle* particle2) {
  double r_ij = 0.0;
  double r_i, r_j;
  double trap_length=m_system->getTraplength();
  int infinity= (2^31) - 1;     //Maximum c++ value

  for(int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
    r_i=particle1->getPosition()[dim];
    r_j=particle2->getPosition()[dim];
    r_ij+= (r_i-r_j)*(r_i-r_j);
  }

  r_ij = sqrt(r_ij);
  
  if(r_ij <= trap_length){
      return infinity;
  }
  else{
      return 0.0;
  }

}
