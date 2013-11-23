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
    return log((sqrt(1.0-2.0*Rho*x+x*x)+x-Rho)/(1.0-Rho));
}

double SABR::Impvolnormal(){
    double gamma1=Beta/Fmid;
    double gamma2=-Beta*(1.0-Beta)/pow(Fmid,2);
    double Delta=(Alpha/(Sigma*(1.0-Beta)))*(pow(F0,1.0-Beta)-pow(K,1.0-Beta));
    if(fabs(F0-K)>=0.00001){
        return Alpha*(F0-K)/D(Delta)*(1.0+epsilon*((2.0*gamma2-gamma1*gamma1)/24.0*pow(Sigma*C(Fmid)/Alpha,2)+
                    Rho*gamma1/4.0*Sigma*C(Fmid)/Alpha+(2.0-3.0*pow(Rho,2))/24.0));
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

double SABR::SABRMCcallprice(){
	double resultprice=0;
	double** twocorrelatednv=new double*[2];
	GetTwoCorrelatedNRV(n*m,Rho,twocorrelatednv);
	for(int j=0;j<m;j++){
		double dt=T*1.0/n;
		double resultF=F0;
		double resultsigma=Sigma;
		for(int i=0;i<n;i++){
			int temp1=i+j*n;
			resultF=resultF+resultsigma*pow(resultF,Beta)*pow(dt,0.5)*twocorrelatednv[0][temp1];
			if(resultF<0){
				resultF=0;
			}
			resultsigma=resultsigma*exp(Alpha*pow(dt,0.5)*twocorrelatednv[1][temp1]-pow(Alpha,2)*0.5*dt);
		}
		double tempresult=resultF-K;
		if(tempresult<0){
			tempresult=0;
		}
		resultprice=resultprice+tempresult;
	}
	resultprice=resultprice*1.0/m;
	
	for(int i=0;i<2;i++){
		double* tempptr=twocorrelatednv[i];
		delete[] tempptr;
	}
	delete[] twocorrelatednv;
    return resultprice;
}

double SABR::SABRMCputprice(){
	double resultprice=0;
	double** twocorrelatednv=new double*[2];
	GetTwoCorrelatedNRV(n*m,Rho,twocorrelatednv);
	for(int j=0;j<m;j++){
		double dt=T*1.0/n;
		double resultF=F0;
		double resultsigma=Sigma;
		for(int i=0;i<n;i++){
			int temp1=i+j*n;
			resultF=resultF+resultsigma*pow(resultF,Beta)*pow(dt,0.5)*twocorrelatednv[0][temp1];
			if(resultF<0){
				resultF=0;
			}
			resultsigma=resultsigma*exp(Alpha*pow(dt,0.5)*twocorrelatednv[1][temp1]-pow(Alpha,2)*0.5*dt);
		}
		double tempresult=-resultF+K;
		if(tempresult<0){
			tempresult=0;
		}
		resultprice=resultprice+tempresult;
	}
	resultprice=resultprice*1.0/m;
	
	for(int i=0;i<2;i++){
		double* tempptr=twocorrelatednv[i];
		delete[] tempptr;
	}
	delete[] twocorrelatednv;
    return resultprice;
}

matrix<double> SABRNormalVolCube(double Sigma,double Alpha, double Beta, double Rho,
        double maturity[],double Strike[],double ForwardSwapRate[],int m, int n){
    matrix<double> result(m,n);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            SABR sabrobj(ForwardSwapRate[i],ForwardSwapRate[i]+Strike[j],maturity[i],Sigma,Alpha,Beta,Rho);
            result(i,j)=sabrobj.Impvolnormal();
        }
    }
    return result;
}

