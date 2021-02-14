#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>

#include "WaveFunctions/simplegaussian.h"


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}


bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

     //_________BRUTE FORCE_____________

     int random_index;
     double psi_factor;
     double wfold=m_waveFunction->evaluate(m_particles);
     //double wfold=m_stepLength;
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
       //Try moving this random thing in the forloop
       //double temp_rand=UniformNumberGenerator(gen);
       //double step=m_stepLength*(temp_rand-0.5);
       m_particles[random_index]->adjustPosition(step, dim);
     }

     //std::cout << wfold<<std::endl;

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


bool System::metropolisStepImportanceSampling() {
    //_________Importance sampling_____________

    //Declaring vaiables to be used:
    double part_1, part_2, green_factor, step, greenRate=0;
    //double TS=m_system->getTimeStep();
    //Defining position and quantum force vectors
    //to be used in the importance sampling
    std::vector<double> PositionOld=std::vector<double>();
    std::vector<double> QFOld=std::vector<double>();
    std::vector<double> PositionNew=std::vector<double>();
    std::vector<double> QFNew=std::vector<double>();

    //Random integer generator
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<int> distribution(0,m_numberOfParticles-1);
    std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);
    std::normal_distribution<double> Normaldistribution(0.0,1.0);

    double rand_norm  =Normaldistribution(gen);
    double random_dist=UniformNumberGenerator(gen);
    //Random index used to choose a random particle
    int random_index  =distribution(gen);

    //Defining the values of the previous position
    double wfold=m_waveFunction->evaluate(m_particles);
    PositionOld=m_particles[random_index]->getPosition();
    QFOld=m_waveFunction->computeQuantumForce(PositionOld);

    for (int dim=0; dim<m_numberOfDimensions; dim++){
        step=QFOld[dim]*m_timeStep*0.5 + sqrt(m_timeStep)*rand_norm;
        m_particles[random_index]->adjustPosition(step, dim);
    }

    // Evaluate new quantities
    double wfnew = m_waveFunction->evaluate(m_particles);
    PositionNew=m_particles[random_index]->getPosition();
    QFNew=m_waveFunction->computeQuantumForce(PositionNew);

    // Compute greens function
    for (int dim=0; dim<m_numberOfDimensions; dim++){
        part_1=PositionOld[dim]-PositionNew[dim]-0.5*m_timeStep*QFNew[dim];
        part_2=PositionNew[dim]-PositionOld[dim]-0.5*m_timeStep*QFOld[dim];
        greenRate+=(part_2*part_2)-(part_1*part_1);
    }
    greenRate = exp(greenRate/(2*m_timeStep));
    green_factor = greenRate*wfnew*wfnew/(wfold*wfold);

    // Check if the step is accepted
    if (random_dist <= green_factor) {
        wfold = wfnew;
        return true;
    }
    else {
        m_particles[random_index]->setPosition(PositionOld);
        return false;
    }


    //Extracting the new wavefunction, and checks if it is accepted
/*


    for (int cycles = 1; cycles <= NumberMCsamples; cycles++){
   // new position
   for (int i = 0; i < NumberParticles; i++) {
     for (int j = 0; j < Dimension; j++) {
       // gaussian deviate to compute new positions using a given timestep
       NewPosition(i,j) = OldPosition(i,j) + Normaldistribution(gen)*sqrt(timestep)+OldQuantumForce(i,j)*timestep*D;

     }}}

     //Choosing a random step:
     double temp_rand=UniformNumberGenerator(gen);
     double step=m_stepLength*(temp_rand-0.5);

     //Start the step which gives movement of the particle
     for (int dim=0; dim<m_numberOfDimensions; dim++){
       m_particles[random_index]->adjustPosition(step, dim);
     }

     //std::cout << wfold<<std::endl;

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

*/
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    bool acceptedStep;
    for (int i=0; i < numberOfMetropolisSteps; i++) {
      if (m_bruteforce==true){
        acceptedStep = metropolisStep();
      }
      else{
        acceptedStep = metropolisStepImportanceSampling();
      }

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */

         //If statement to send the accepted steps into the sampler
         //after the system is at rest
         if (i>=numberOfMetropolisSteps*m_equilibrationFraction){
            m_sampler->sample(acceptedStep);
         }
    }

    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
    //m_sampler->writeToFile();
}


double System::gradientDecent(){

    //Remember to put this in main and setting the alpha
    //before running the main metropolis. Run for about 1000 steps
    double alpha_curr=0.45;
    double alpha_prev=alpha_curr;  //Initial guess

    double alpha_next;
    double El_deriv2, El_deriv_prev;
    double gamma_dir;
    double epsilon;

    int iterations=10;   //Not sure how many iterations to use, 50 or 500??

    double beta = getWaveFunction()->getParameters()[2]/getWaveFunction()->getParameters()[0];

    for(int ii=0; ii<iterations; ii++){
        //std::vector<double> parameters(2);
        //parameters[0] = alpha_curr;
        //parameters[2] = alpha_curr*beta;
        //Need to find a way to set the parameters public for each itteration
        setWaveFunction(new SimpleGaussian(system,alpha_curr, beta));

        El_deriv2=E_LDerivative(alpha_curr);
        El_deriv_prev=E_LDerivative(alpha_prev);

        gamma_dir=(alpha_curr-alpha_prev+epsilon)/(El_deriv2-El_deriv_prev);

        alpha_next=alpha_curr-gamma_dir*(El_deriv2);

        alpha_prev=alpha_curr;
        alpha_curr=alpha_next;
    }

    //To compute the local energy derivative, write a function that runs the
    //runmetropolisstep function and returns the derivative function of local energy.


    return alpha_curr;

}

double System::E_LDerivative(double alpha_n){
    runMetropolisSteps(1000);
    double expectEnerg, expectE_L_deri, expectderi_dot_EL
                      =m_sampler->getGradientDecentValues();

    double expectEnergy    = expectEnerg/(m_numberOfMetropolisSteps);//*getEquilibrationFraction());
    double expectE_L_deriv     = expectE_L_deri/(m_numberOfMetropolisSteps);//*getEquilibrationFraction());
    double expectderiv_dot_EL = expectderi_dot_EL/(m_numberOfMetropolisSteps);//*getEquilibrationFraction());

    double E_L_derive=2*(expectderiv_dot_EL - expectEnergy*expectE_L_deriv)/m_equilibrationFraction;
    return  E_L_derive;
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

void System::setNumeric(bool numeric) {
    m_numeric = numeric;
}

void System::setBruteforce(bool bruteforce_val) {
    m_bruteforce = bruteforce_val;
}

void System::setTimeStep(double timeStep) {
    m_timeStep= timeStep;
}
