#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles)=0;
    virtual double computeDoubleDerivativeNumeric (std::vector<class Particle*> particles);
    virtual double computePotentialEnergy(std::vector<class Particle*> particles)=0;
    virtual double computePotentialEnergyInteracting(std::vector<class Particle*> particles)=0;
    virtual double computeEnergyInteracting(class Particle* particle1, class Particle* particle2)=0;

protected:
    class System* m_system = nullptr;
};
