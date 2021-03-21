#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../particle.h"
#include "../system.h"

using std::cout;
using std::endl;

RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    //Module to genrate the random numbers for placement of particles
    std::random_device rd;
    std::mt19937_64 gen(rd());

    // Set up the distribution for x \in [[0, 1],(can use multiple configurations)
    std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);

    //Looping over the dimensions and particles initializing each
    //particle with a random placement using a uniform distribution
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int dim=0; dim < m_numberOfDimensions; dim++) {
             // This is where the particles are placed in positions
             // using the uniform generator
             double temp_pos=(UniformNumberGenerator(gen) - 0.5);
             position.push_back(temp_pos);
        }
        //Setting the number of dimensions and appending new particle postions
        m_particles.push_back(new Particle());
        m_particles[i]->setNumberOfDimensions(m_numberOfDimensions);
        m_particles[i]->setPosition(position);
    }
}
