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

using namespace std;

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
     //double temp_rand=UniformNumberGenerator(gen);
     double step;

     //Start the step which gives movement of the particle
     for (int dim=0; dim<m_numberOfDimensions; dim++){
        step=m_stepLength*(UniformNumberGenerator(gen)-0.5);
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
     if (UniformNumberGenerator(gen)<=psi_factor){
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

}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    bool acceptedStep;

    setHistogram();

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
            m_sampler->sample(acceptedStep, i);
         }
         cout<<"The current step is "<<i<<endl;
    }

    //Chooosing what to sample
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
    if (m_general_wtf==true){m_sampler->writeToFile();}
    //if (m_GDwtf==true){m_sampler->writeToFileAlpha();}
    if (m_obdwtf==true){m_sampler->writeToFileOBD();}
    
}

double System::gradientDescent(double initialAlpha){
//Gradient descent method to find the optimal variational parameter alpha given an initial parameter initialAlpha
    int steepestDescentSteps = (int) 1e+4;
    int maxIterations=50;
    double alpha = initialAlpha;
    double beta = getWaveFunction()->getParameters()[1];
    double lambda = -0.005;
    int iterations = 0;
    double energyDerivative;
    double tol = 1e-5;
    vector<double> parameterss(2);
    m_GDalpha.push_back(alpha);

    while (iterations < maxIterations){
        parameterss[0] = alpha;
        parameterss[1] = beta;
        getWaveFunction()->setParameters(parameterss);

        runMetropolisSteps(steepestDescentSteps);
        
        energyDerivative = getSampler()->Energy_Der2();

        alpha += lambda*energyDerivative;

        cout<<energyDerivative<<endl;

        m_energyarr.push_back(m_sampler->getEnergy());
        m_GDalpha.push_back(alpha);
        iterations++;

        cout<< " New alpha = "  << alpha <<  endl;
        cout<< " Iterations = " << iterations << endl;

        if (fabs(energyDerivative) < tol){
            cout<<endl;
            cout<<"Alpha stabilized: Tol reached before iterator"<<endl;
            break;
        }
    }
    cout<<endl;

    if (m_GDwtf==true){m_sampler->writeToFileAlpha();}

    cout<<"Performing metropolisrun with best alpha"<<endl;

    //alpha = cumulativeAlpha / (maxIterations*percentAlphasToSave);
    return alpha;
}

void System::checkStep(double stepLength, double timeStep){
    setStepLength               (stepLength);
    setTimeStep                 (timeStep);
    int steps = 0;
    int maxsteps=1000000;
    vector<int> steps_list=vector<int>();
    vector<double> meanEL_list=vector<double>();

    while (steps < maxsteps){
        runMetropolisSteps(steps);

        steps_list.push_back(steps);
        meanEL_list.push_back(m_sampler->getEnergy());

        cout<< " Running: "  << steps*100/maxsteps << " %" << endl;

        steps+=10000;
    }

    m_sampler->writeToFileSteps(steps_list, meanEL_list);

}



void System::oneBodyDensity(){
    //Function to make the histrograms needed to compute the one body density
    vector<int> histogram(m_bins);
    int bucket;
    double r2 = 0;
    for (int j=0; j<getNumberOfParticles(); j++){
        r2 = 0;

        for (int d=0; d<getNumberOfDimensions(); d++){
            r2 += m_particles.at(j)->getPosition()[d]*m_particles.at(j)->getPosition()[d];
        }
        r2 = sqrt(r2);
        bucket = int(floor(r2/m_bucketSize));
        histogram[bucket] += 1;
    }

    // Update histogram
    for (int k=0; k<m_bins; k++){
       m_histogram[k] += histogram[k];
    }
}

/*
void System::setOneBodyDensity(double min, double max, int numberOfBins) {
    m_numberOfBins = numberOfBins;
    m_min = min;
    m_max = max;
    m_binWidth = (max - min) / numberOfBins;
    m_bins = (double**) calloc(m_system->getNumberOfDimensions(), sizeof(double*));
    for (int i=0; i<m_system->getNumberOfDimensions(); i++) {
        m_bins[i] = (double*) calloc(m_numberOfBins, sizeof(double));
    }
}
*/

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

void System::setInteraction(bool interaction) {
    m_interaction= interaction;
}

void System::setTraplength(double a_length) {
    m_a_length= a_length;
}

void System::setGD(bool GD) {
    m_GD = GD;
}

void System::setGDwtf(bool GDwtf) {
    m_GDwtf = GDwtf;
}

void System::setgeneralwtf(bool generalwtf) {
    m_general_wtf = generalwtf;
}

void System::setobd(bool obdwtf, double bucketSize, int bins){
    m_obdwtf = obdwtf;
    m_bucketSize=bucketSize;
    m_bins=bins;
}

void System::setHistogram(){
    vector<int> histogram(getBins());
    m_histogram = histogram;
}