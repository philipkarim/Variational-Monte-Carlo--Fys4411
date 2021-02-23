#pragma once

class GradientDecent {
public:
    GradientDecent();
    //class GradientDecent*  getGradientDecent()        { return m_gradientdecent; }

    //class Sampler* getGradientDecentValues(){ return grad_list; }
    void runGradientDecent(double alpha_init, double omega, double beta,
                          double timeStep, double stepLength, double
                          equilibration, int numberOfDimensions, int
                          numberOfParticles, int seed, bool numeric, bool
                          bruteforce_val);
    double E_LDerivative(double alpha_n, double omega, double beta,
    double timeStep, double stepLength, double equilibration, int numberOfDimensions, int numberOfParticles,
                 int seed, bool numeric, bool bruteforce_val);
    class Sampler*  getSampler()        { return m_sampler; }



private:
    class Sampler*                  m_sampler = nullptr;
    class System*                   m_system = nullptr;

    //class GradientDecent*           m_gradientdecent = nullptr;
    //class Sampler*              grad_list = nullptr;
    //class Sampler*                  getSampler()        { return m_sampler; }

};

/*
extern double omega, beta, timeStep, stepLength, equilibration;
extern int numberOfDimensions, numberOfParticles, seed;
extern bool numeric, bruteforce_val;
*/
