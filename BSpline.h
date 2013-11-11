#ifndef BSPLINE_H_
#define BSPLINE_H_
#include<tuple>
#include<vector>

class BSpline{
public:
	static int N;
	static double *knot;
	static std::vector<std::tuple<int, int, double> >calculated;
	static std::vector<double> BFvalues;
	static double* resultf;
	static double* resultl;
	static double* resultOIS;

	BSpline();
	double BF(int k, int d, double t);  //Basis Function
	double BF_integrate(int k, int d, double t); //Integrate the basis function from negative infinity to t
	double BF_integrate(int k, int d, double S, double T);  //Integrate the basis function from S to T.
	double BF_derivative(int k, int d, double t);     //First derivative of Basis Function
	double BF_2derivative(int k, int d, double t);    //Second derivative of Basis Function
	double BF_3derivative(int k, int d, double t);    //Third derivative of Basis Function
	double BF_2integral(int k, int l, int d, double S, double T);  //The integral of cross terms of two Basis functions, k,
	double Tikhonov_sum_squared(double f[], double l[]); //The Tikhonov regularizer function
	double Tikhonov_sum_squared();

	// Functionalities
	double DF(double t1, double t2, double f[]);
	double DF(double t1, double t2);
	double LiborForward(double t1, double t2, double l[]);
	double LiborForward(double t1, double t2);
	double OISForward(double t1, double t2, double OIS[]);
	double OISForward(double t1, double t2);
	double SwapRate(double t1, double t2, double f[], double l[]);
	double SwapRate(double t1, double t2);
	double BasisSwap(double t1, double t2, double f[], double l[], double OIS[]);
	double BasisSwap(double t1, double t2);

	//Calculate swap values based on curves.
	double LiborSwapValue(double t1, double t2, double rate, double f[], double l[]);
	double LiborSwapValue(double t1, double t2, double rate);
	double BasisSwapValue(double t1, double t2, double rate, double f[], double l[], double OIS[]);
	double BasisSwapValue(double t1, double t2, double rate);
};

#endif
