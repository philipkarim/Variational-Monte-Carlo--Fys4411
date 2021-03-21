#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>
#include <vector>

using namespace std;

SimpleGaussian::SimpleGaussian(System* system, double alpha, double beta) :
        WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_numberOfParameters = m_parameters.size();
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
     //Implementation of the gaussian wave function. 
     //This functions evaluates the trial wave function returning the value
     //of PSI_T

    // Defining some variables to be used
    double g_func=1;
    double x_part, y_part, z_part;
    std::vector<double> r_pos_E;
    std::vector<int> dimensions_length_ev(m_system->getNumberOfDimensions());
    std::iota(dimensions_length_ev.begin(), dimensions_length_ev.end(), 0);

    //G_function of the trial wave function
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
      r_pos_E=particles[i]->getPosition();

      if (dimensions_length_ev.size()==3){
        x_part=r_pos_E[0]*r_pos_E[0];
        y_part=r_pos_E[1]*r_pos_E[1];
        z_part=m_parameters[1]*r_pos_E[2]*r_pos_E[2];
      }
      else if(dimensions_length_ev.size()==2){
        x_part=r_pos_E[0]*r_pos_E[0];
        y_part=r_pos_E[1]*r_pos_E[1];
        z_part=0;
      }
      else {
        x_part=r_pos_E[0]*r_pos_E[0];
        y_part=0;
        z_part=0;
      }

      g_func*=exp(-m_parameters[0]*(x_part+y_part+z_part));
    }

    //F_function of the trial wave function, if the interaction is off, f=1
    double f_func;
    if (m_system->getInteraction()==true){
      f_func=1.0;
      //Looping over the product og each particle correlating with eachother
      for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
          for(int j = i+1; j < m_system->getNumberOfParticles(); j++) {
              f_func*=correlation(particles[i],particles[j]);
          }
      }
    }
else{
  f_func=1.0;
}
    //Return a double value
    return g_func*f_func;
}


double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
     //Computes the value of the analytical double derivative for the non interacting case. 

    //Defining som values and variables to be used in the calculations  
     double DD_val=0;
     double wf=evaluate(particles);
     vector<double> r_pos_CDD;

    //Looping over the particles and dimensions to calculating the
    //Double derivative for each particle in each dimension
     for(int i=0; i<m_system->getNumberOfParticles(); i++){
       r_pos_CDD=particles[i]->getPosition();
       for (int dim = 0; dim<m_system->getNumberOfDimensions(); dim++){
         if (dim==2){
           DD_val+=m_parameters[1]*(2*m_parameters[0]*m_parameters[1]*r_pos_CDD[dim]*r_pos_CDD[dim]-1);
         }
         else{
          DD_val+=(2*m_parameters[0]*r_pos_CDD[dim]*r_pos_CDD[dim]-1);
         }
       }
     }

    //Return a double value
    return -(DD_val*wf*m_parameters[0]);
  }


std::vector<double> SimpleGaussian::computeQuantumForce (std::vector<double> particles){
    //Calculating the drift force used by the Metropolis-Hastings algorithm

    //Defining som values and variables to be used in the calculations  
    double x_i=particles[0];
    double y_i=particles[1];
    double z_i=particles[2];
    std::vector<double> QF_vec =  std::vector<double>();
    std::vector<int> dimensions_length_CQF(m_system->getNumberOfDimensions());
    std::iota(dimensions_length_CQF.begin(), dimensions_length_CQF.end(), 0);

    //Tried to make the loop by implementing some if statements to see if that was the reason for a
    //an error. (Kept the if statements, made the code quiet easy to read)
    if (dimensions_length_CQF.size()==3){
      QF_vec.push_back(-4*x_i*m_parameters[0]*exp(m_parameters[0]*(x_i+y_i+m_parameters[1]*z_i)));
      QF_vec.push_back(-4*y_i*m_parameters[0]*exp(m_parameters[0]*(x_i+y_i+m_parameters[1]*z_i)));
      QF_vec.push_back(-4*z_i*m_parameters[1]*m_parameters[0]*exp(m_parameters[0]*(x_i+y_i+m_parameters[1]*z_i)));
    }
    else if(dimensions_length_CQF.size()==2){
      QF_vec.push_back(-4*x_i*m_parameters[0]*exp(m_parameters[0]*(x_i+y_i+m_parameters[1]*z_i)));
      QF_vec.push_back(-4*y_i*m_parameters[0]*exp(m_parameters[0]*(x_i+y_i+m_parameters[1]*z_i)));
      z_i=0;
    }
    else {
      QF_vec.push_back(-4*x_i*m_parameters[0]*exp(m_parameters[0]*(x_i+y_i+m_parameters[1]*z_i)));
      y_i=0;
      z_i=0;
    }
    //Returns a vector consisting of up to 3 values
    return QF_vec;
}


double SimpleGaussian::correlation(Particle* particle1, Particle* particle2) {
    //Function used to evaluate the psi function.
    //Two particles are sent into the function, then the forces between them
    //are returned

    //Defining som values and variables to be used in the calculations  
    double r_ij = 0.0;
    double r_i, r_j;
    double trap_length=m_system->getTraplength();

    //Looping over the dimensions of the two particles
    for(int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
      r_i=particle1->getPosition()[dim];
      r_j=particle2->getPosition()[dim];
      r_ij+= (r_i-r_j)*(r_i-r_j);
    }

    r_ij = sqrt(r_ij);

    //Returning the forces between them
    if(r_ij <= trap_length) {
        return 0.0;
    } else {
        return 1.0-(trap_length/r_ij);
    }

}