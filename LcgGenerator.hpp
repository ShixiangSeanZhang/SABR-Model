#ifndef LCGGENERATOR_HPP_
#define LCGGENERATOR_HPP_

class LcgGenerator {
private:
    long x;
    int a;
    int c;
    int k;
    long uniformRandomNumberCounter;

public:
    LcgGenerator(long x, int a, int c, int k);
    virtual ~LcgGenerator();
    long getUniformRandomNumberCounter() const;
    double nextDouble();
};

#endif /* LCGGENERATOR_HPP_ */
