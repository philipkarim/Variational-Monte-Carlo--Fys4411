#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double beta);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivativeNumeric(std::vector<class Particle*> particles);
    std::vector<double> computeQuantumForce (std::vector<double> particles);
    double correlation(class Particle* particle1, class Particle* particle2);
    double computeDoubleDerivativeInteraction(std::vector<class Particle*> particles);
    std::vector<double> computeGradient(std::vector<Particle *> particles);

};
