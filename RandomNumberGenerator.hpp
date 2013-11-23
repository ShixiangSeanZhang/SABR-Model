
#ifndef RANDOMNUMBERGENERATOR_HPP_
#define RANDOMNUMBERGENERATOR_HPP_

class RandomNumberGenerator {

public:
    virtual ~RandomNumberGenerator() {
    }

    virtual long getUniformRandomNumberCounter() = 0;
    virtual long getStdNormalRandomNumberCounter() = 0;
    virtual double nextDouble()=0;
};

#endif /* RANDOMNUMBERGENERATOR_HPP_ */
