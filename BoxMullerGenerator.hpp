
#ifndef BOXMULLERGENERATOR_HPP_
#define BOXMULLERGENERATOR_HPP_

#include "LcgGenerator.hpp"
#include "RandomNumberGenerator.hpp"

struct RandomNumberPair {
    double Z1;
    double Z2;
};

class BoxMullerGenerator: public RandomNumberGenerator {
private:
    LcgGenerator uniformDistrubtionGenerator;
    long stdNormalRandomNumberCounter;
    bool needsNewPair;
    RandomNumberPair pair;

    RandomNumberPair generateNewPair();

public:
    BoxMullerGenerator(LcgGenerator uniformDistrubtionGenerator);
    virtual ~BoxMullerGenerator();
    long getUniformRandomNumberCounter();
    long getStdNormalRandomNumberCounter();
    double nextDouble();
};

void GetTwoCorrelatedNRV(int N,double rho,double** twodimarray);

#endif /* BOXMULLERGENERATOR_HPP_ */
