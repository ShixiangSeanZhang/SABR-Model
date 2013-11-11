#ifndef BSMODEL_H
#define BSMODEL_H

class BSMODEL{
    private:
        double T;
        double K;
        double F;
        double Sigma;
        BSMODEL(){}

    public:
        BSMODEL(double Fa,double Ka,double Ta, double Sigmaa):T(Ta),K(Ka),F(Fa),Sigma(Sigmaa){}
        double GetCDF(double x) const;
        double Getpdf(double x) const;
        double Bcall();
        double Bput();
        double Bcallnorm();
        double Bputnorm();
};




#endif
