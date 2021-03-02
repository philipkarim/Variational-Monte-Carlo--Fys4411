#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, double omega_z, double gamma);
    double computeLocalEnergy(std::vector<Particle*> particles);
    double computePotentialEnergy(std::vector<Particle*> particles);
    double computePotentialEnergyInteracting(std::vector<Particle*> particles);
    double computeEnergyInteracting(class Particle* particle1, class Particle* particle2);
private:
    double m_omega = 0;
    double m_omega_z = 0;
    double m_gamma = 0;


};
