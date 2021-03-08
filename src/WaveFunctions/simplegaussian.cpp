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

    //Add a beta value also and make the r_tot depend on beta and set beta=0
    //to get harmonic thing
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

     //Calculating a simple gaussian function, looping over each particle
     //in each dimension.

// Newest

    double r_tot=0.0;
    double g_func=1;
    std::vector<double> r_pos_E;

    double x_part, y_part, z_part;
    std::vector<int> dimensions_length_ev(m_system->getNumberOfDimensions());
    std::iota(dimensions_length_ev.begin(), dimensions_length_ev.end(), 0);


/*
    for(int i=0; i<m_system->getNumberOfParticles(); i++){
      //r_pos_E=particles[i]->getPosition();
      for (int dim=0; dim<m_system->getNumberOfDimensions(); dim++){
        if (dim==2){
          r_tot+=m_parameters[1]*particles[i]->getPosition()[dim]*particles[i]->getPosition()[dim];
        }
        else{
          r_tot+=particles[i]->getPosition()[dim]*particles[i]->getPosition()[dim];
        }
      }
      g_func*=exp(-m_parameters[0]*r_tot);
    }

*/

//G_func
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

//F_func
double f_func;
if (m_system->getInteraction()==true){
  f_func=1.0;

   for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
       for(int j = i+1; j < m_system->getNumberOfParticles(); j++) {
           f_func*=correlation(particles[i],particles[j]);
       }
   }
}
else{
  f_func=1.0;
}
/*
//interaction
    double r_jk;

    for (int j=0; j<m_system->getNumberOfParticles();j++)
    {
      for (int k=j+1; k<m_system->getNumberOfParticles();k++)
      {
        r_jk=m_system->getR_jk(j,k);
        if(r_jk>m_a){
        f=f*( 1 - m_a/m_system->getR_jk(j,k) );}
        else {return 0; }
      }
    }

    return g*f;
*/

    return g_func*f_func;
}

//Analytical differentiation
double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

     //Non interacting, probably have to add a longer bit of code and an
     //if statement when adding the interacting case.
     //This is analytically and only works on simple gaussion functions

     // For non-interacting particles
     //std::vector<double> r_pos_CDD;
     //r_pos_CDD.reserve(m_system->getNumberOfDimensions());
     double wf=evaluate(particles);
     double derivate2=0;
     double x_lap, y_lap, z_lap;
     vector<double> r_pos_CDD;
     std::vector<int> dimensions_length_CDD(m_system->getNumberOfDimensions());
     std::iota(dimensions_length_CDD.begin(), dimensions_length_CDD.end(), 0);


     //________________first term in (2)__________
     //Analytically differentiatian
     //if(m_system->getNumeric()==false){
     double DD_val=0;
     for(int i=0; i<m_system->getNumberOfParticles(); i++){

       //r_pos_CDD.clear();
       r_pos_CDD=particles[i]->getPosition();
       //r_pos_CDD.push_back(temp);

       for (int dim = 0; dim<m_system->getNumberOfDimensions(); dim++){
         if (dim==2){
           DD_val+=m_parameters[1]*(2*m_parameters[0]*m_parameters[1]*r_pos_CDD[dim]*r_pos_CDD[dim]-1);
         }
         else{
          DD_val+=(2*m_parameters[0]*r_pos_CDD[dim]*r_pos_CDD[dim]-1);
         }
       }
     }


     //interacting kinetic energy
/*   HAve a look at this
     double a=m_system->getinteractionSize();

//Interaction terms
for(int i=0; i<m_system->getNumberOfParticles(); i++){
    double r_i_square=0;
    for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
     r_i_square += particles.at(i).getPosition()[d]*
                   particles.at(i).getPosition()[d];
     }
     int d = m_system->getNumberOfDimensions()-1;
     r_i_square += particles.at(i).getPosition()[d]*
                   particles.at(i).getPosition()[d]*m_parameters[1]);

   double second=0;
   double third=0;
   double fourth=0;
   double fifth=0;
   double temp;

   for(int j=0; j < i; j++) {

       double r_ij = m_system->getDistanceMatrixij(i,j);

       temp= a / ( (r_ij-a) * r_ij );

       second += temp;

       double r_ir_j = 0;
       for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d];
       }
       int d = m_system->getNumberOfDimensions() - 1;
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d]*
                     m_parameters[2]/(m_parameters[0]);

       fourth-= temp * temp;

       fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
               ( r_ij );

   }
   for(int j = j+1; j < m_system->getNumberOfParticles(); j++){

       double r_ij = m_system->getDistanceMatrixij(i,j);

       temp = a / ( (r_ij-a) * r_ij );

       second += temp;

       double r_ir_j = 0;
       for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d];
       }
       int d = m_system->getNumberOfDimensions() - 1;
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d]*
                     m_parameters[2]/(m_parameters[0]);

       fourth-= temp * temp;

       fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
               ( r_ij );
   }

   third=second*second;

   EK_interacting+=second+third+fourth+fifth;

}
*/


    return -(DD_val*wf*m_parameters[0]);
  }


/*
       for(int i=0; i<m_system->getNumberOfParticles(); i++){

         //r_pos_CDD.clear();
         r_pos_CDD=particles[i]->getPosition();
         //r_pos_CDD.push_back(temp);


         if (dimensions_length_CDD.size()==3){
           x_lap=(2*m_parameters[0]*r_pos_CDD[0]*r_pos_CDD[0]-1);
           y_lap=(2*m_parameters[0]*r_pos_CDD[1]*r_pos_CDD[1]-1);
           z_lap=m_parameters[1]*(2*m_parameters[0]*m_parameters[1]*r_pos_CDD[2]*r_pos_CDD[2]-1);
         }
         else if(dimensions_length_CDD.size()==2){
           x_lap=(2*m_parameters[0]*r_pos_CDD[0]*r_pos_CDD[0]-1);
           y_lap=(2*m_parameters[0]*r_pos_CDD[1]*r_pos_CDD[1]-1);
           z_lap=0;
         }
         else {
           x_lap=(2*m_parameters[0]*r_pos_CDD[0]*r_pos_CDD[0]-1);
           y_lap=0;
           z_lap=0;
         }

         derivate2+=(x_lap+y_lap+z_lap);
       }

      return -derivate2*wf*m_parameters[0];
    }
*/

/*
    double pos, deriv;
int dim = m_system->getNumberOfDimensions();
int nPart = m_system->getNumberOfParticles();

double wf = evaluate(particles);//evaluate(particles);

for (int i = 0; i<nPart; i++){
    for (int j = 0; j<dim; j++){
        pos = particles[i]->getPosition()[j];
        deriv -= 2* m_parameters[0]*(2* m_parameters[0]*pos*pos)*wf;
    }
    deriv += 2* m_parameters[1]*dim*wf;
}
return deriv;

}
*/
/*
double one=0;
for(int i = 0; i < m_system->getNumberOfParticles(); i++){
    for(int d = 0; d < m_system->getNumberOfDimensions(); d++){
        one += m_parameters[d]*m_parameters[d]*
                particles[i]->getPosition()[d]*
                particles[i]->getPosition()[d];

    }
}

one*=4.0;
one-= 2 * ( (m_system->getNumberOfDimensions() - 1) * m_parameters[0] + m_parameters[2])
        * m_system->getNumberOfParticles(); //constant term

return one;
}
*/

std::vector<double> SimpleGaussian::computeQuantumForce (std::vector<double> particles){
    //Returns a vector consisting of up to 3 particles
    //Only for spherical without correlation
    //double QF=0;

    double x_i=particles[0];
    double y_i=particles[1];
    double z_i=particles[2];

    std::vector<double> QF_vec =  std::vector<double>();
    std::vector<int> dimensions_length_CQF(m_system->getNumberOfDimensions());
    std::iota(dimensions_length_CQF.begin(), dimensions_length_CQF.end(), 0);

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

    return QF_vec;
}


double SimpleGaussian::correlation(Particle* particle1, Particle* particle2) {

    double rij = 0.0;

    for(int dim = 0; dim < m_system->getNumberOfDimensions(); dim++) {
        rij+= (particle1->getPosition()[dim]-particle2->getPosition()[dim])*(particle1->getPosition()[dim]-particle2->getPosition()[dim]);
    }

    rij = sqrt(rij);

    if(rij <= m_system->getTraplength()) {
        return 0.0;
    } else {
        return 1.0-(m_system->getTraplength()/rij);
    }

}


double SimpleGaussian::computeDoubleDerivativeInteraction(std::vector<class Particle*> particles) {
    double interactionEK=0;

    const double alpha =m_parameters[0];
        const double beta  = m_parameters[1];
        const double a = m_system->getTraplength();

        std::vector<double> gradient = computeGradient(particles);


        double laplacian = 0.0;

        for(int k = 0; k < m_system->getNumberOfParticles(); k++) {

            const double xk = particles[k]->getPosition()[0];
            const double yk = particles[k]->getPosition()[1];
            const double zk = particles[k]->getPosition()[2];

            const double phi_rk = exp(-alpha*(xk*xk + yk*yk + beta*zk*zk));
            const double ddphi_k = (-4.0*alpha - 2.0*alpha*beta + (2.0*alpha*xk)*(2.0*alpha*xk)
                                    + (2.0*alpha*yk)*(2.0*alpha*yk)
                                    + (2.0*alpha*beta*zk)*(2.0*alpha*beta*zk))*phi_rk;


            std::vector<double> sum1(3);
            double sum2 = 0.0;
            double sum3 = 0.0;

            for(int j = 0; j < m_system->getNumberOfParticles(); j++) {



                if(j != k) {

                    double r_kj = 0.0;
                    const double xj = particles[j]->getPosition()[0];
                    const double yj = particles[j]->getPosition()[1];
                    const double zj = particles[j]->getPosition()[2];

                    r_kj = (xk-xj)*(xk-xj) + (yk-yj)*(yk-yj) + (zk-zj)*(zk-zj);
                    r_kj = sqrt(r_kj);

                    double du_kj = a/(r_kj*r_kj - a*r_kj);
                    double ddu_kj = -(a*(2*r_kj - a))/((r_kj*r_kj - a*r_kj)*(r_kj*r_kj - a*r_kj));

                    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
                        if(i != k) {

                            double r_ki     = 0.0;
                            const double xi = particles[i]->getPosition()[0];
                            const double yi = particles[i]->getPosition()[1];
                            const double zi = particles[i]->getPosition()[2];

                            r_ki = (xk-xi)*(xk-xi) + (yk-yi)*(yk-yi) + (zk-zi)*(zk-zi);
                            r_ki = sqrt(r_ki);

                            double du_ki = a/(r_ki*r_ki - a*r_ki);

                            double rkri_dot_rkrj = (xk-xi)*(xk-xj) + (yk-yi)*(yk-yj) + (zk-zi)*(zk-zj);

                            sum2+= ((rkri_dot_rkrj)/(r_ki*r_kj))*du_ki*du_kj;

                        }
                    }


                    sum1[0]+= (xk-xj)*(du_kj/r_kj);
                    sum1[1]+= (yk-yj)*(du_kj/r_kj);
                    sum1[2]+= (zk-zj)*(du_kj/r_kj);

                    sum3   += ddu_kj + (2.0/r_kj)*du_kj;

                }
            }

            const double grad_xk = gradient[m_system->getNumberOfDimensions()*k+0];
            const double grad_yk = gradient[m_system->getNumberOfDimensions()*k+1];
            const double grad_zk = gradient[m_system->getNumberOfDimensions()*k+2];

            double grad_rk_dot_sum1 = 2.0*(grad_xk*sum1[0] + grad_yk*sum1[1] + grad_zk*sum1[2]);
            laplacian+= (ddphi_k/phi_rk) + (grad_rk_dot_sum1/phi_rk) + sum2 + sum3;

        }

        return -0.5*laplacian;
    }

std::vector<double> SimpleGaussian::computeGradient(std::vector<Particle *> particles) {
    std::vector<double> gradient(m_system->getNumberOfDimensions()*m_system->getNumberOfParticles());
    double alpha = m_parameters[0];

    for(int p = 0; p < m_system->getNumberOfParticles(); p++) {
        for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
            gradient[m_system->getNumberOfDimensions()*p+d] = -2.0*alpha*particles[p]->getPosition()[d]*evaluate(particles);
        }
    }

    return gradient;

}
