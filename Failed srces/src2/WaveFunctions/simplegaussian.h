#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double beta);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivativeNumeric(std::vector<class Particle*> particles);
    std::vector<double> computeQuantumForce (std::vector<double> particles);
    //double E_LDerivative                    (double alpha_n);

    private:

        class Sampler*                  m_sampler = nullptr;
        class System* m_system = nullptr;

};
