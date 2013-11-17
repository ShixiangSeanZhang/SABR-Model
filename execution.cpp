#include<iostream>
#include"BSModel.h"
#include"SABR.h"
#include"BSpline.h"
#include"Swaption.h"
using namespace std;

extern BSpline* bs=new BSpline[1];//no construct so new

int main(){
    BSpline::N = 20;
    double knot_[] = {-0.2,-0.13,-0.06,0.0,0.02,0.06,0.15,0.3,0.8,1.6,3.1,6.1,10.1,15.1,23.1,30.1,35.0,40.0,45.0,50.0};
    BSpline::knot = knot_;
    cout.precision(12);
    double constf[]={0.000802405,0.000904085,0.000985058,0.001084941,0.001255748,0.001555918,0.001460851,0.002745368,0.013591273,0.036095721,0.033091786,0.032107426,0.031337984,0.021476918,0.015297221,0.016};
    double constl[]={0.005214085,0.005420273,0.005585051,0.005791275,0.006175834,0.007201412,0.00790712,0.005750676,0.016674008,0.033516058,0.032679726,0.034492664,0.025407151,0.032514395,0.03103284,0.033043666};
    BSpline::resultf=constf;
    BSpline::resultl=constl;
    BSMODEL BSobj(1.0,1.0,1.0,0.2);
    cout<<BSobj.GetCDF(0)<<endl;
    cout<<BSobj.Bcall()<<endl;
    cout<<BSobj.Bput()<<endl;
    cout<<BSobj.Bcallnorm()<<endl;
    cout<<BSobj.Bputnorm()<<endl;
    cout<<"Secant"<<BSobj.SecantSigmaLNtoN()<<endl;
    //F0,K,T,Sigma,Alpha,Beta,Rho
    SABR SABRobj(0.02,0.02,0.25,0.02,0.004,0.5,0.2);
    cout<<"SABR"<<endl;
    cout<<SABRobj.Impvolnormal()<<endl;
    cout<<SABRobj.ImpPutPricen(SABRobj.Impvolnormal())<<endl;
    cout<<SABRobj.ImpCallPricen(SABRobj.Impvolnormal())<<endl;
    SWAPTION Swaptionobj(0.25,1.25,0.03,0.003);
    cout<<Swaptionobj.A(0.17)<<endl;
    cout<<Swaptionobj.PrecN()<<endl;
    cout<<Swaptionobj.PpayN()<<endl;
    cout<<Swaptionobj.PrecLN()<<endl;
    cout<<Swaptionobj.PpayLN()<<endl;

}
