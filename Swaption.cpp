#include"Swaption.h"
#include"math.h"
#include"BSpline.h"
#include"BSModel.h"
#include<iostream>
extern BSpline* bs;

SWAPTION::SWAPTION(double T0a, double Ta, double Ka, double Sigmaa):
         T0(T0a),T(Ta),K(Ka),Sigma(Sigmaa){}


double SWAPTION::A(double Tval){
    int n=int((T-T0)/0.25);
    double result=0.0;
    for(int i=0;i<n;i++){
        double df=bs->DF(Tval,T0+(i+1)*0.25);
        result=0.25*df+result;
    }
    return result;
}

double SWAPTION::PrecN(){
    double F=bs->SwapRate(T0,T);
    BSMODEL bsobj(F,K,T,Sigma);
    return A(0.0)*bsobj.Bputnorm();
}

double SWAPTION::PpayN(){
    double F=bs->SwapRate(T0,T);
    BSMODEL bsobj(F,K,T,Sigma);
    return A(0.0)*bsobj.Bcallnorm();
}

double SWAPTION::PrecLN(){
    double F=bs->SwapRate(T0,T);
    BSMODEL bsobj(F,K,T,Sigma);
    return A(0.0)*bsobj.Bput();
}

double SWAPTION::PpayLN(){
    double F=bs->SwapRate(T0,T);
    BSMODEL bsobj(F,K,T,Sigma);
    return A(0.0)*bsobj.Bcall();
}
