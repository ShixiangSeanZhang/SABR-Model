
#include<cppinterface.h>
#include "BSpline.h"
#include "BSModel.h"
#include "SABR.h"
#include "Swaption.h"

#pragma warning (disable : 4996)
using namespace System;


// Libor & OIS Curve
BSpline *bs = new BSpline[1];

double InitializeBSpline(MyMatrix &knots){
	int N = knots.rows();
	double *knotarray = new double[N];
	for (int i = 0; i < N; i++){
		knotarray[i] = knots(i,0);
	}
	BSpline::N = N;
	BSpline::knot =knotarray;
	return 0.0;
}

double InitializeCurves(MyMatrix &f,
						MyMatrix &l,
						MyMatrix &OIS){
							int M = f.rows();
							double *OISParam = new double[M];
							for (int i = 0; i < M; i++){
								OISParam[i] = OIS(i,0);
							}
							BSpline::resultOIS = OISParam;

							double *lParam = new double[M];
							for (int i = 0; i < M; i++){
								lParam[i] = l(i,0);
							}
							BSpline::resultl = lParam;

							double *fParam = new double[M];
							for (int i = 0; i < M; i++){
								fParam[i] = f(i,0);
							}
							BSpline::resultf = fParam;
							return 0.0;
}

double BFToExcel(int k,
				 int d,
				 double t){
					 return bs->BF(k,d,t);
}

double BF_derivativeToExcel(int k,
							int d,
							double t){
								return bs->BF_derivative(k,d,t);
}

double BF_2derivativeToExcel(int k,
							 int d,
							 double t){
								 return bs->BF_2derivative(k,d,t);
}

double BF_3derivativeToExcel(int k,
							 int d,
							 double t){
								 return bs->BF_3derivative(k,d,t);
}


double BF_integrateToExcel(int k,
						   int d,
						   double t){
							   return bs->BF_integrate(k,d,t);
}

double BF_integrate2paramToExcel(int k,
								 int d,
								 double s,
								 double t){
									 return bs->BF_integrate(k,d,s,t);
}

double BF_2integralToExcel(int k, 
						   int l,
						   int d,
						   double s,
						   double t){
							   return bs->BF_2integral(k,l,d,s,t);
}

double Tikhonov_sum_squaredToExcel(MyMatrix &fparam,
								   MyMatrix &lparam){
									   int M = fparam.rows();
									   double *f = new double[M];
									   double *l = new double[M];
									   for (int i = 0; i < M; i++){
										   f[i] = fparam(i,0);
										   l[i] = lparam(i,0);
									   }
									   return bs->Tikhonov_sum_squared(f,l);
}

double DFToExcel(double t1,
				 double t2,
				 MyMatrix &fparam){
					 int M = fparam.rows();
					 double *f = new double[M];
					 for (int i = 0; i < M; i++){
						 f[i] = fparam(i,0);
					 }
					 return bs->DF(t1, t2, f);
}

double DFToExcel_(double t1,
				  double t2){
					  return bs->DF(t1, t2, BSpline::resultf);
}


double LiborForwardToExcel(double t1, 
						   double t2,
						   MyMatrix &lparam){
							   int M = lparam.rows();
							   double *l = new double[M];
							   for (int i = 0; i < M; i++){
								   l[i] = lparam(i,0);
							   }
							   return bs->LiborForward(t1, t2, l);
}

double LiborForwardToExcel_(double t1, 
							double t2){
								return bs->LiborForward(t1, t2, BSpline::resultl);
}

double OISForwardToExcel(double t1,
						 double t2,
						 MyMatrix &OISparam){
							 int M = OISparam.rows();
							 double *OIS = new double[M];
							 for (int i = 0; i < M; i++){
								 OIS[i] = OISparam(i,0);
							 }
							 return bs->OISForward(t1, t2, OIS);
}

double OISForwardToExcel_(double t1,
						  double t2){
							  return bs->OISForward(t1, t2, BSpline::resultOIS);
}

double SwapRateToExcel(double t1,
					   double t2,
					   MyMatrix &fparam,
					   MyMatrix &lparam){
						   int M = fparam.rows();
						   double *f = new double[M];
						   double *l = new double[M];
						   for (int i = 0; i < M; i++){
							   f[i] = fparam(i,0);
							   l[i] = lparam(i,0);
						   }
						   return bs->SwapRate(t1,t2,f,l);
}

double SwapRateToExcel_(double t1,
						double t2){
							return bs->SwapRate(t1,t2,BSpline::resultf,BSpline::resultl);
}

double BasisSwapToExcel(double t1,
						double t2,
						MyMatrix &fparam,
						MyMatrix &lparam,
						MyMatrix &OISparam){
							int M = fparam.rows();
							double *f = new double[M];
							double *l = new double[M];
							double *OIS = new double[M];
							for (int i = 0; i < M; i++){
								f[i] = fparam(i,0);
								l[i] = lparam(i,0);
								OIS[i] = OISparam(i,0);
							}
							return bs->BasisSwap(t1,t2,f,l,OIS);
}

double BasisSwapToExcel_(double t1,
						 double t2){
							 return bs->BasisSwap(t1,t2,BSpline::resultf,BSpline::resultl,BSpline::resultOIS);
}


double LiborSwapValueToExcel(double t1,
							 double t2,
							 double rate,
							 MyMatrix &fparam,
							 MyMatrix &lparam){
								 int M = fparam.rows();
								 double *f = new double[M];
								 double *l = new double[M];
								 for (int i = 0; i < M; i++){
									 f[i] = fparam(i,0);
									 l[i] = lparam(i,0);
								 }
								 return bs->LiborSwapValue(t1, t2, rate, f, l);
}

double LiborSwapValueToExcel_(double t1,
							  double t2,
							  double rate){
								  return bs->LiborSwapValue(t1, t2, rate, BSpline::resultf, BSpline::resultl);
}

double BasisSwapValueToExcel(double t1,
							 double t2,
							 double rate,
							 MyMatrix &fparam,
							 MyMatrix &lparam,
							 MyMatrix &OISparam){
								 int M = fparam.rows();
								 double *f = new double[M];
								 double *l = new double[M];
								 double *OIS = new double[M];
								 for (int i = 0; i < M; i++){
									 f[i] = fparam(i,0);
									 l[i] = lparam(i,0);
									 OIS[i] = OISparam(i,0);
								 }
								 return bs->BasisSwapValue(t1, t2, rate, f, l, OIS);
}

double BasisSwapValueToExcel_(double t1,
							  double t2,
							  double rate){
								  return bs->BasisSwapValue(t1, t2, rate, BSpline::resultf, BSpline::resultl, BSpline::resultOIS);
}

//Black Model
//double F0_, double K_, double T_, double Sigma_

double BcallToExcel(MyMatrix &BlackParams){
	double F0_ = BlackParams(0,0);
	double K_ = BlackParams(1,0);
	double T_ = BlackParams(2,0);
	double Sigma_ = BlackParams(3,0);
	BSMODEL bsobj(F0_, K_, T_, Sigma_);
	return bsobj.Bcall();
}

double BputToExcel(MyMatrix &BlackParams){
	double F0_ = BlackParams(0,0);
	double K_ = BlackParams(1,0);
	double T_ = BlackParams(2,0);
	double Sigma_ = BlackParams(3,0);
	BSMODEL bsobj(F0_, K_, T_, Sigma_);
	return bsobj.Bput();
}

double BcallnormToExcel(MyMatrix &BlackParams){
	double F0_ = BlackParams(0,0);
	double K_ = BlackParams(1,0);
	double T_ = BlackParams(2,0);
	double Sigma_ = BlackParams(3,0);
	BSMODEL bsobj(F0_, K_, T_, Sigma_);
	return bsobj.Bcallnorm();
}

double BputnormToExcel(MyMatrix &BlackParams){
	double F0_ = BlackParams(0,0);
	double K_ = BlackParams(1,0);
	double T_ = BlackParams(2,0);
	double Sigma_ = BlackParams(3,0);
	BSMODEL bsobj(F0_, K_, T_, Sigma_);
	return bsobj.Bputnorm();
}

double SecantSigmaLNtoNToExcel(MyMatrix &BlackParams){
	double F0_ = BlackParams(0,0);
	double K_ = BlackParams(1,0);
	double T_ = BlackParams(2,0);
	double Sigma_ = BlackParams(3,0);
	BSMODEL bsobj(F0_, K_, T_, Sigma_);
	return bsobj.SecantSigmaLNtoN();
}


//SABR Model
//double F0_, double K_, double T_, double sigma0_, double alpha_, double beta_, double rho_
double ImpCallPricenToExcel(MyMatrix &SABRParam,double impvol){
	double F0_ = SABRParam(0,0);
	double K_ = SABRParam(1,0);
	double T_ = SABRParam(2,0);
	double sigma0_ = SABRParam(3,0);
	double volvol_ = SABRParam(4,0);
	double beta_ = SABRParam(5,0);
	double rho_ = SABRParam(6,0);
	SABR sabrobj(F0_, K_, T_, sigma0_, volvol_, beta_, rho_);
	return sabrobj.ImpCallPricen(impvol);
}

double ImpPutPricenToExcel(MyMatrix &SABRParam,double impvol){
	double F0_ = SABRParam(0,0);
	double K_ = SABRParam(1,0);
	double T_ = SABRParam(2,0);
	double sigma0_ = SABRParam(3,0);
	double volvol_ = SABRParam(4,0);
	double beta_ = SABRParam(5,0);
	double rho_ = SABRParam(6,0);
	SABR sabrobj(F0_, K_, T_, sigma0_, volvol_, beta_, rho_);
	return sabrobj.ImpPutPricen(impvol);
}

double ImpvolnormalToExcel(MyMatrix &SABRParam){
	double F0_ = SABRParam(0,0);
	double K_ = SABRParam(1,0);
	double T_ = SABRParam(2,0);
	double sigma0_ = SABRParam(3,0);
	double volvol_ = SABRParam(4,0);
	double beta_ = SABRParam(5,0);
	double rho_ = SABRParam(6,0);
	SABR sabrobj(F0_, K_, T_, sigma0_, volvol_, beta_, rho_);
	return sabrobj.Impvolnormal();
}

double SABRMCcallpriceToExcel(MyMatrix &SABRParam,int m,int n){
	double F0_ = SABRParam(0,0);
	double K_ = SABRParam(1,0);
	double T_ = SABRParam(2,0);
	double sigma0_ = SABRParam(3,0);
	double volvol_ = SABRParam(4,0);
	double beta_ = SABRParam(5,0);
	double rho_ = SABRParam(6,0);
	int m_=m;
	int n_=n;
	SABR sabrobj(F0_, K_, T_, sigma0_, volvol_, beta_, rho_,m_,n_);
	return sabrobj.SABRMCcallprice();
}

double SABRMCputpriceToExcel(MyMatrix &SABRParam,int m,int n){
	double F0_ = SABRParam(0,0);
	double K_ = SABRParam(1,0);
	double T_ = SABRParam(2,0);
	double sigma0_ = SABRParam(3,0);
	double volvol_ = SABRParam(4,0);
	double beta_ = SABRParam(5,0);
	double rho_ = SABRParam(6,0);
	int m_=m;
	int n_=n;
	SABR sabrobj(F0_, K_, T_, sigma0_, volvol_, beta_, rho_,m_,n_);
	return sabrobj.SABRMCputprice();
}

//Swaption Model
//Swaption(double T0_, double T_, double K_, double Sigma_);
double AToExcel(MyMatrix& SwaptionParams,double Tval){
	double T0_ = SwaptionParams(0,0);
	double T_ = SwaptionParams(1,0);
	double K_ = SwaptionParams(2,0);
	double Sigma_ = SwaptionParams(3,0);
	SWAPTION swtobj(T0_, T_, K_, Sigma_);
	return swtobj.A(Tval);
}

double PrecNToExcel(MyMatrix& SwaptionParams){
	double T0_ = SwaptionParams(0,0);
	double T_ = SwaptionParams(1,0);
	double K_ = SwaptionParams(2,0);
	double Sigma_ = SwaptionParams(3,0);
	SWAPTION swtobj(T0_, T_, K_, Sigma_);
	return swtobj.PrecN();
}

double PpayNToExcel(MyMatrix& SwaptionParams){
	double T0_ = SwaptionParams(0,0);
	double T_ = SwaptionParams(1,0);
	double K_ = SwaptionParams(2,0);
	double Sigma_ = SwaptionParams(3,0);
	SWAPTION swtobj(T0_, T_, K_, Sigma_);
	return swtobj.PpayN();
}

double PrecLNToExcel(MyMatrix& SwaptionParams){
	double T0_ = SwaptionParams(0,0);
	double T_ = SwaptionParams(1,0);
	double K_ = SwaptionParams(2,0);
	double Sigma_ = SwaptionParams(3,0);
	SWAPTION swtobj(T0_, T_, K_, Sigma_);
	return swtobj.PrecLN();
}

double PpayLNToExcel(MyMatrix& SwaptionParams){
	double T0_ = SwaptionParams(0,0);
	double T_ = SwaptionParams(1,0);
	double K_ = SwaptionParams(2,0);
	double Sigma_ = SwaptionParams(3,0);
	SWAPTION swtobj(T0_, T_, K_, Sigma_);
	return swtobj.PpayLN();
}

MyMatrix SABRNormalVolCubeToExcel(double Sigma,double Alpha, double Beta, double Rho, MyMatrix &maturity,MyMatrix &Strike,MyMatrix &ForwardSwapRate){
	using namespace boost::numeric::ublas;
	int m=maturity.size1();
	int n=Strike.size2();
	MyMatrix result(m,n);
	double* Strikeptr=new double[n];
	double* maturityptr=new double[m];
	double* Fsrptr=new double[m];
	for(int i=0;i<n;i++){
		Strikeptr[i]=Strike(0,i);
	}
	for(int i=0;i<m;i++){
		maturityptr[i]=maturity(i,0);
	}
	for(int i=0;i<m;i++){
		Fsrptr[i]=ForwardSwapRate(i,0);
	}
	matrix<double> ublasresult(m,n);
	ublasresult=SABRNormalVolCube(Sigma, Alpha, Beta, Rho, maturityptr,Strikeptr,Fsrptr,m,n);
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			result(i,j)=ublasresult(i,j);
		}
	}
	return result;
}