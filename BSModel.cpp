#include"BSModel.h"
#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions.hpp>
#include<math.h>
using namespace boost::math;


double BSMODEL::GetCDF(double x) const{
    normal_distribution<>myNorm;
    return cdf(myNorm,x);
}

double BSMODEL::Getpdf(double x) const{
    double PI=constants::pi<double>();
    return (1.0/sqrt(2*PI))*exp(-pow(x,2)/2.0);
}

double BSMODEL::Bcall(){
    double d1=(log((F*1.0)/K)+0.5*pow(Sigma,2)*T)*1.0/(Sigma*sqrt(T*1.0));
    double d2=d1-Sigma*sqrt(T);
    return F*GetCDF(d1)-K*GetCDF(d2);
}

double BSMODEL::Bput(){
    double d1=(log((F*1.0)/K)+0.5*pow(Sigma,2)*T)*1.0/(Sigma*sqrt(T*1.0));
    double d2=d1-Sigma*sqrt(T);
    return -F*GetCDF(-d1)+K*GetCDF(-d2);
}

double BSMODEL::Bcallnorm(){
    double d1=((F*1.0)-K)/(Sigma*sqrt(T*1.0));
    return (d1*GetCDF(d1)+Getpdf(d1))*Sigma*sqrt(T*1.0);
}

double BSMODEL::Bputnorm(){
    double d2=-((F*1.0)-K)/(Sigma*sqrt(T*1.0));
    return (d2*GetCDF(d2)+Getpdf(d2))*Sigma*sqrt(T*1.0);
}

double BSMODEL::ForSigmaLNtoN(double x){
    BSMODEL bsobj(F,K,T,x);
    return Bcall()-bsobj.Bcallnorm();
}

double BSMODEL::SecantSigmaLNtoN(){
    double xold1 = 0.01;
    double xold2 = 5.0;
    double xnew = 0.2;
    while(fabs(xnew-xold2)>=0.000001){
        double temp = xnew;
        xnew = xnew - ForSigmaLNtoN(xnew)*(xold2-xold1)/(ForSigmaLNtoN(xold2)-ForSigmaLNtoN(xold1));
        xold1 = xold2;
        xold2 = temp;
    }
    return xnew;
}

