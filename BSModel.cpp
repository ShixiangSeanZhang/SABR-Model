#include"BSModel.h"
#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions.hpp>
#include<math.h>
using namespace boost::math;


double BSMODEL::GetCDF(double x) const{
    normal_distribution<>myNorm;
    return cdf(myNorm,x);
}

double BSMODEL::Bcall(){
    double d1=(log((F*1.0)/K)+0.5*pow(Sigma,2)*T)*1.0/(Sigma*sqrt(T*1.0));
    double d2=d1-Sigma*sqrt(T);
    return F*GetCDF(d1)-K*GetCDF(d2);
}

double BSMODEL::Bput(){
    double d1=(log((F*1.0)/K)+0.5*pow(Sigma,2)*T)*1.0/(Sigma*sqrt(T*1.0));
    double d2=d1-Sigma*sqrt(T);
    return -F*GetCDF(-d1)-K*GetCDF(-d2);
}
