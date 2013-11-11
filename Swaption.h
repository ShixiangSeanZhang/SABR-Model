#ifndef SWAPTION_H
#define SWAPTION_H

class SWAPTION{
    private:
        double T0;
        double T;
        double Sigma;
        double K;
    public:
        SWAPTION(double T0a,double Ta,double Ka,double Sigmaa);
        double A(double Tval);
        double PrecN();
        double PpayN();
        double PrecLN();
        double PpayLN();
};

#endif
