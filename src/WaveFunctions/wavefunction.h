#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles)=0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual double correlation(class Particle* particle1, class Particle* particle2)=0;
    virtual std::vector<double> computeQuantumForce (std::vector<double> particles)=0;
    void setParameters(const std::vector<double> &parameters);

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};
