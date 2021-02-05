#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, double alpha);
    double computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
    double alpha =0;
};
