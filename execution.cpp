#include<iostream>
#include"BSModel.h"
#include"SABR.h"
using namespace std;

int main(){
    BSMODEL BSobj(1.0,1.0,1.0,0.2);
    cout<<BSobj.GetCDF(0)<<endl;
    cout<<BSobj.Bcall()<<endl;
    cout<<BSobj.Bput()<<endl;
    cout<<BSobj.Bcallnorm()<<endl;
    cout<<BSobj.Bputnorm()<<endl;
    SABR SABRobj(1.0,1.0,1.0,0.2,0.002,0.5,0.25);
    cout<<SABRobj.Impvolnormal()<<endl;
    //impvol=SABRobj.Impvolnormal();
    cout<<SABRobj.ImpPutPricen(SABRobj.Impvolnormal())<<endl;
    cout<<SABRobj.ImpCallPricen(SABRobj.Impvolnormal())<<endl;
    SABR SABRobj2(1.0,1.2,1.0,0.2,0.002,0.5,0.25);
    cout<<SABRobj2.Impvolnormal()<<endl;
}
