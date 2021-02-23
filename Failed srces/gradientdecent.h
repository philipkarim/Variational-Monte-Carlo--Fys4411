#pragma once
#include <vector>



class GD {
public:
    void runGradientDecent(std::vector<double> input_variables_double, std::vector<int> input_variables_int, std::vector<bool> input_variables_bool);
    double E_LDerivative(double alpha_n, std::vector<double> in_double, std::vector<int> in_int, std::vector<bool> in_bool);

private:

    class System* m_system = nullptr;
    class Sampler* m_sampler = nullptr;

};

















//class GradientDecent {
//public:
  //  GradientDecent();
    //class GradientDecent*  getGradientDecent()        { return m_gradientdecent; }

    //class Sampler* getGradientDecentValues(){ return grad_list; }
//void runGradientDecent(std::vector<double> input_variables_double, std::vector<int> input_variables_int, std::vector<bool> input_variables_bool);


//double E_LDerivative(double alpha_n, std::vector<double> in_double, std::vector<int> in_int, std::vector<bool> in_bool);

//class Sampler*             getGradientDecentValues()   { return grad_list; }


//public:
  //class Sampler*                  getSampler()        { return m_sampler; }
//private:
  //class Sampler*                  m_sampler = nullptr;


    //class Sampler*  getSampler()        { return m_sampler; }
//class Sampler*                  getSampler()        { return m_sampler; }
//class Sampler*                  m_sampler = nullptr;

/*
private:
    class Sampler*                  m_sampler = nullptr;
    class System*                   m_system = nullptr;

    //class GradientDecent*           m_gradientdecent = nullptr;
    //class Sampler*              grad_list = nullptr;
    //class Sampler*                  getSampler()        { return m_sampler; }

};
*/
/*
extern double omega, beta, timeStep, stepLength, equilibration;
extern int numberOfDimensions, numberOfParticles, seed;
extern bool numeric, bruteforce_val;
*/
