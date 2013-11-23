#include <iostream>
#include "LcgGenerator.hpp"

using namespace std;

LcgGenerator::LcgGenerator(long x, int a, int c, int k) :
        x(x), a(a), c(c), k(k), uniformRandomNumberCounter(0) {
}

LcgGenerator::~LcgGenerator() {
}

double LcgGenerator::nextDouble() {
    x = (a * x + c) % k;
    uniformRandomNumberCounter++;
    return x / (double) k;
}

long LcgGenerator::getUniformRandomNumberCounter() const {
    return uniformRandomNumberCounter;
}
