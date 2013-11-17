#include"SABR.h"
#include"BSModel.h"
#include<math.h>
using namespace std;
using namespace boost::numeric::ublas;

double SABR::C(double x) const{
    return pow(x,Beta);
}

double SABR::C1(double x) const{
    return Beta*pow(x,Beta-1.0);
}

double SABR::D(double x) const{
    return log(sqrt(1.0-2.0*Rho*x+x*x)+x-Rho)/(1.0-Rho);
}

double SABR::Impvolnormal(){
    double gamma1=Beta/Fmid;
    double gamma2=-Beta*(1.0-Beta)/pow(Fmid,2);
    double Delta=(Alpha/(Sigma*(1.0-Beta)))*(pow(F0,1.0-Beta)-pow(K,1.0-Beta));
    if(fabs(F0-K)>=0.00001){
        return Alpha*(F0-K)/D(Delta)*(1.0+epsilon*((2.0*gamma2-gamma1*gamma1)/24.0*pow(Sigma*C(Fmid)/Alpha,2)+
                    Rho*gamma1/4.0*Sigma*C(Fmid)/Alpha+(2.0-3.0*pow(Rho,2)/24.0)));
    }else{
        return Sigma/pow(F0,-Beta)*(1.0+epsilon*((2.0*gamma2-gamma1*gamma1)/24.0*pow(Sigma*C(F0)/Alpha,2)+
                    Rho*gamma1/4.0*Sigma*C(F0)/Alpha+(2.0-3.0*pow(Rho,2)/24.0)));
    }
}

double SABR::ImpCallPricen(double impvol){
    BSMODEL bsobj(T,K,F0,impvol);
    return bsobj.Bcallnorm();
}

double SABR::ImpPutPricen(double impvol){
    BSMODEL bsobj(T,K,F0,impvol);
    return bsobj.Bputnorm();
}

matrix<double> SABRNormalVolCube(double Sigma,double Alpha, double Beta, double Rho,
        double maturity[],double Strike[],double ForwardSwapRate[]){
    int m = sizeof(maturity)/sizeof(double);
    int n = sizeof(Strike)/sizeof(double);
    matrix<double> result(m,n);
    for(int i=0;i++;i<m){
        for(int j=0;j++;j<n){
            SABR sabrobj(ForwardSwapRate[i],ForwardSwapRate[i]+Strike[j],maturity[i],Sigma,Alpha,Beta,Rho);
            result(i,j)=sabrobj.Impvolnormal();
        }
    }
    return result;
}
