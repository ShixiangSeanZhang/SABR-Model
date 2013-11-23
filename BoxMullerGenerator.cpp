
#include "BoxMullerGenerator.hpp"
#include <iostream>
#include <cmath>

using namespace std;

BoxMullerGenerator::BoxMullerGenerator(LcgGenerator uniformDistrubtionGenerator) :
        uniformDistrubtionGenerator(uniformDistrubtionGenerator),  stdNormalRandomNumberCounter(0),
        needsNewPair(true) {
}

BoxMullerGenerator::~BoxMullerGenerator() {
}

long BoxMullerGenerator::getUniformRandomNumberCounter() {
    return uniformDistrubtionGenerator.getUniformRandomNumberCounter();
}

long BoxMullerGenerator::getStdNormalRandomNumberCounter() {
    return stdNormalRandomNumberCounter;
}

double BoxMullerGenerator::nextDouble() {
    if (needsNewPair) {
        pair = generateNewPair();
        needsNewPair = false;
        stdNormalRandomNumberCounter++;
        return pair.Z1;
    } else {
        needsNewPair = true;
        stdNormalRandomNumberCounter++;
        return pair.Z2;
    }
}

RandomNumberPair BoxMullerGenerator::generateNewPair() {
    double X = 2;
    double u1 = 0;
    double u2 = 0;

    while (X > 1) {
        u1 = uniformDistrubtionGenerator.nextDouble();
        u2 = uniformDistrubtionGenerator.nextDouble();
        u1 = 2 * u1 - 1;
        u2 = 2 * u2 - 1;
        X = pow(u1, 2) + pow(u2, 2);
    }

    double Y = sqrt(-2 * log(X) / X);

    RandomNumberPair pair;
    pair.Z1 = u1 * Y;
    pair.Z2 = u2 * Y;

    return pair;
}

void GetTwoCorrelatedNRV(int N,double rho,double** twodimarray){
    double* Z1series=new double[N];
	double* Z2series=new double[N];
	LcgGenerator lcgobj(1L, 39373, 0, pow(2,31)-1);
	BoxMullerGenerator boxmullerobj(lcgobj);
	for(int i=0;i<N;i++){
		Z1series[i]=boxmullerobj.nextDouble();
		Z2series[i]=rho*Z1series[i]+sqrt(1-pow(rho,2))*boxmullerobj.nextDouble();
	}
	twodimarray[0] = Z1series;
	twodimarray[1] = Z2series;
}

