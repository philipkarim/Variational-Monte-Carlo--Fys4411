#pragma once
#include <vector>
#include <Math/random.h>

class System {
public:
    System();
    System(int seed);
    //Added aftertremoved with typos?
    /*
    void setPosition(const std::vector<double> &position);
    void adjustPosition(double change, int dimension);
    void setNumberOfDimensions(int numberOfDimensions);
    std::vector<double> getPosition() { return m_position; }
    */
    //
    bool metropolisStep             ();
    bool metropolisStepImportanceSampling();
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    class Random*                   getRandomEngine()   { return m_random; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }

    void setNumeric                     (bool numeric);
    bool getNumeric()                   { return m_numeric; }

    void setBruteforce                  (bool bruteforce_val);
    bool getBruteforce()                {return m_bruteforce;}
    void setTimeStep                    (double timeStep);
    double getTimeStep()                {return m_timeStep;}
    void setInteraction                    (bool interaction);
    double getInteraction()                {return m_interaction;}
    void setTraplength                  (double a_length);
    double getTraplength()                {return m_a_length;}

    double gradientDescent              (double initialAlpha);
    double findEnergyDerivative();


    //double gradientDecent();
    //double E_LDerivative                (double alpha_n);
    //bool getNumeric()          { return m_numeric; }



private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    class Random*                   m_random = nullptr;

    //true, false or nothing?
    bool m_numeric;
    bool m_bruteforce;
    double m_timeStep;
    bool m_interaction;
    double m_a_length;
};
