#pragma once

class main {

public:
    class GradientDecent*  getGradientDecent()        { return m_gradientdecent; }
    //class Sampler* getGradientDecentValues(){ return grad_list; }
    //virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;

    virtual void runGradientDecent(double alpha_init, double omega, double beta,
                          double timeStep, double stepLength, double
                          equilibration, int numberOfDimensions, int
                          numberOfParticles, int seed, bool numeric, bool
                          bruteforce_val)=0;

private:
    class GradientDecent*             m_gradientdecent = nullptr;
    //class GradientDecent*           m_gradientdecent = nullptr;
    //class Sampler*              grad_list = nullptr;
    //class Sampler*                  getSampler()        { return m_sampler; }


/*
extern double omega, beta, timeStep, stepLength, equilibration;
extern int numberOfDimensions, numberOfParticles, seed;
extern bool numeric, bruteforce_val;
*/
};
