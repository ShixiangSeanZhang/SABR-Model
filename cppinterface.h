//
//
//                                                                    Test.h
//

#ifndef TEST_H
#define TEST_H

#include "xlw/MyContainers.h"
#include <xlw/CellMatrix.h>
#include <xlw/DoubleOrNothing.h>
#include <xlw/ArgList.h>
#include "BSModel.h"
#include "BSpline.h"
#include "SABR.h"
#include "Swaption.h"

using namespace xlw;

//<xlw:libraryname=MyManagedTestLibrary

//BSpline
double InitializeBSpline(MyMatrix &knots);

double InitializeCurves(MyMatrix &f,
						MyMatrix &l,
						MyMatrix &OIS);

double BFToExcel(int k,
				 int d,
				 double t);

double BF_derivativeToExcel(int k,
							int d,
							double t);

double BF_2derivativeToExcel(int k,
							 int d,
							 double t);

double BF_3derivativeToExcel(int k,
							 int d,
							 double t);

double BF_integrateToExcel(int k,
						   int d,
						   double t);

double BF_integrate2paramToExcel(int k,
								 int d,
								 double s,
								 double t);


double BF_2integralToExcel(int k, 
						   int l,
						   int d,
						   double s,
						   double t);


double Tikhonov_sum_squaredToExcel(MyMatrix &fparam,
								   MyMatrix &lparam);


double DFToExcel(double t1,
				 double t2,
				 MyMatrix &fparam);

double DFToExcel_(double t1,
				  double t2);

double LiborForwardToExcel(double t1, 
						   double t2,
						   MyMatrix &lparam);

double LiborForwardToExcel_(double t1, 
							double t2);

double OISForwardToExcel(double t1,
						 double t2,
						 MyMatrix &OISparam);

double OISForwardToExcel_(double t1,
						  double t2);

double SwapRateToExcel(double t1,
					   double t2,
					   MyMatrix &fparam,
					   MyMatrix &lparam);

double SwapRateToExcel_(double t1,
						double t2);




double BasisSwapToExcel(double t1,
						double t2,
						MyMatrix &fparam,
						MyMatrix &lparam,
						MyMatrix &OISparam);

double BasisSwapToExcel_(double t1,
						 double t2);


double LiborSwapValueToExcel(double t1,
							 double t2,
							 double rate,
							 MyMatrix &fparam,
							 MyMatrix &lparam);

double LiborSwapValueToExcel_(double t1,
							  double t2,
							  double rate);



double BasisSwapValueToExcel(double t1,
							 double t2,
							 double rate,
							 MyMatrix &fparam,
							 MyMatrix &lparam,
							 MyMatrix &OISparam);

double BasisSwapValueToExcel_(double t1,
							  double t2,
							  double rate);


//BSModel

double BcallToExcel(MyMatrix &BlackParams);
double BputToExcel(MyMatrix &BlackParams);
double BcallnormToExcel(MyMatrix &BlackParams);
double BputnormToExcel(MyMatrix &BlackParams);
double SecantSigmaLNtoNToExcel(MyMatrix &BlackParams);


//SABR Model
double ImpCallPricenToExcel(MyMatrix &SABRParam,double impvol);
double ImpPutPricenToExcel(MyMatrix &SABRParam,double impvol);
double ImpvolnormalToExcel(MyMatrix &SABRParam);
double SABRMCcallpriceToExcel(MyMatrix &SABRParam,int m,int n);
double SABRMCputpriceToExcel(MyMatrix &SABRParam,int m,int n);
MyMatrix SABRNormalVolCubeToExcel(double Sigma, double Alpha, double Beta, double Rho,
        MyMatrix &maturity,MyMatrix &Strike,MyMatrix &ForwardSwapRate);

//Swaption Model
//Swaption(double T0_, double T_, double K_, double Sigma_);
double AToExcel(MyMatrix& SwaptionParams,double Tval);
double PrecNToExcel(MyMatrix& SwaptionParams);
double PpayNToExcel(MyMatrix& SwaptionParams);
double PrecLNToExcel(MyMatrix& SwaptionParams);
double PpayLNToExcel(MyMatrix& SwaptionParams);

#endif
