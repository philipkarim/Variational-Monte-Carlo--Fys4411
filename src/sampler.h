#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    void writeToFile();
    double computeVariance(std::vector<double> x_sample, double x_mean);
    double getEnergy()          { return m_energy; }
    double getVariance()          { return m_variance; }
    double getAcceptRatio()          { return m_acceptRatio; }


private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergy2 = 0;
    double  m_variance = 0;
    double  m_acceptedSteps = 0;
    double  m_acceptRatio = 0;
    double time_sec;
    std::vector<double> energy_vec = std::vector<double>();
    class System* m_system = nullptr;
};
