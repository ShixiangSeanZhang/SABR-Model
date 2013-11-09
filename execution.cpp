#include<iostream>
#include"BSModel.h"
using namespace std;

int main(){
    BSMODEL BSobj(1.0,1.0,1.0,0.2);
    cout<<BSobj.GetCDF(0)<<endl;
    cout<<BSobj.Bcall()<<endl;
    cout<<BSobj.Bput()<<endl;
}
