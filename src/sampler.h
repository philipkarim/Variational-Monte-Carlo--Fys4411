#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    void writeToFile();
    void writeToFileAlpha();
    void writeToFileSteps(std::vector<int> steps_list, std::vector<double> meanEL_list);
    void writeToFileOBD();
    double computeVariance(std::vector<double> x_sample, double x_mean);
    double getEnergy()          { return m_energy; }
    double getVariance()          { return m_variance; }
    double getAcceptRatio()          { return m_acceptRatio; }
    
    //double getCumulativeEnergyDerivAvg()          { return m_E_Lderiv; }
    //double getCumulativeEnergyDerivExpectAvg()          { return m_E_Lderiv_expect; }

    double Energy_Der2()          { return 2*(m_E_Lderiv_expect-(m_E_Lderiv*m_energy)); }


    //double getGradientDecentValues()        { return m_cumulativeEnergy, m_cumulativeE_Lderiv, m_cumulativeE_Lderiv_expect; }

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
    double m_E_Lderiv=0;
    double m_E_Lderiv_expect=0;
    double m_cumulativeE_Lderiv=0;
    double m_cumulativeE_Lderiv_expect=0;

    std::vector<double> meanenergy_list = std::vector<double>();

    class System* m_system = nullptr;
};
