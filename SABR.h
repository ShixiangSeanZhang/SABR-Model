#ifndef SABR_H
#define SABR_H
#include<math.h>
#include"BoxMullerGenerator.hpp"
#include<boost/numeric/ublas/io.hpp>
#include<boost/numeric/ublas/matrix.hpp>

class SABR{
    private:
        double F0;
        double K;
        double T;
        double Sigma;
        double Alpha;
        double Beta;
        double Rho;
        double epsilon;
        double Fmid;
		double m;//MC step
		double n;//MC path
        SABR(){};

    public:
		SABR(double F0a, double Ka,double Ta, double Sigmaa, double Alphaa, double Betaa, double Rhoa):
            F0(F0a),K(Ka),T(Ta),Sigma(Sigmaa),Alpha(Alphaa),Beta(Betaa),Rho(Rhoa),epsilon(Alphaa*Alphaa*Ta),Fmid(sqrt(F0*K)){}
        SABR(double F0a, double Ka,double Ta, double Sigmaa, double Alphaa, double Betaa, double Rhoa,int ma, int na):
            F0(F0a),K(Ka),T(Ta),Sigma(Sigmaa),Alpha(Alphaa),Beta(Betaa),Rho(Rhoa),epsilon(Alphaa*Alphaa*Ta),Fmid(sqrt(F0*K)),m(ma),n(na){}
        double C(double x) const;
        double C1(double x) const;
        double D(double x) const;
        double ImpCallPricen(double impvol);
        double ImpPutPricen(double impvol);
        double Impvolnormal();
		double SABRMCcallprice();
		double SABRMCputprice();


};

boost::numeric::ublas::matrix<double> SABRNormalVolCube(double Sigma, double Alpha, double Beta, double Rho,
        double maturity[],double Strike[],double ForwardSwapRate[],int m, int n);


#endif
