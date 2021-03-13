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
    double getStepLength()               {return m_stepLength;}
    void setInteraction                 (bool interaction);
    double getInteraction()             {return m_interaction;}
    void setTraplength                  (double a_length);
    double getTraplength()              {return m_a_length;}
    void checkStep                      (double stepLength, double timeStep);
    void setGD                          (bool GD);
    double getGD()                      {return m_GD;}
    void setGDwtf                       (bool GDwtf);
    void setgeneralwtf                  (bool generalwtf);
    void setobd                         (bool obdwtf, double bucketSize, int bins);
    int getBins()                       {return m_bins;}
    double getBucketSize()              {return m_bucketSize;}
    bool getObd()                       {return m_obdwtf;}
    std::vector<double>get_GDalpha()    { return m_GDalpha; }
    std::vector<double>get_energyarr()  { return m_energyarr; }

    double gradientDescent              (double initialAlpha);
    void oneBodyDensity();

    std::vector<int>m_histogram;
    void setHistogram();
    std::vector<int>getHistogram()      { return m_histogram; }




private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength=0.5;   //It said=0.1
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    class Random*                   m_random = nullptr;
    std::vector<double>             m_GDalpha;
    std::vector<double>             m_energyarr;


   //Just some variables, mostly bools
    bool m_numeric;
    bool m_bruteforce;
    double m_timeStep=0.25;
    bool m_interaction;
    double m_a_length;
    bool m_GD;
    bool m_general_wtf;
    bool m_GDwtf;
    bool m_obdwtf;
    double m_bucketSize;
    int m_bins;

};
