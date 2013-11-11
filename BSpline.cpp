//  BSpline Library
//  September 24, 2013

#include<cmath>
#include<iostream>
#include<cstring>
#include<iomanip>
#include<algorithm>
#include"BSpline.h"

using namespace std;

int BSpline::N = 0;
double *BSpline::knot = NULL;
std::vector<std::tuple<int, int, double> > BSpline::calculated;
std::vector<double> BSpline::BFvalues;
double *BSpline::resultf = NULL;
double *BSpline::resultl = NULL;
double *BSpline::resultOIS = NULL;

/*
BSpline::~BSpline(){
delete knot;
}
*/
BSpline::BSpline(){}

double BSpline::BF(int k, int d, double t){
	// check t lies in the valid interval
	if (k < -3){
		cout << "k cannot be more negative than -3" << endl;
		exit(1);
	}
	if (k > BSpline::N){
		cout << "K cannot be more than " << BSpline::N << endl;
		exit(1);
	}
	if (t < BSpline::knot[k + 3] || t >= BSpline::knot[k+d+4]) {
		return 0.0;
	}
	// t lies in the valid interval and d == 0
	if (d==0) {
		return 1.0;
	}

	// check whether the Basis Function is already cached.
	tuple<int, int, double>tup(k,d,t);
	std::vector<tuple<int, int, double> >::iterator tupindex;
	tupindex = find(BSpline::calculated.begin(), BSpline::calculated.end(), tup);
	if (tupindex != BSpline::calculated.end()){
		return BSpline::BFvalues[tupindex - BSpline::calculated.begin()];
	}


	// t lies in the valid interval and d > 0
	double alpha, beta;
	alpha = (t - BSpline::knot[k + 3])/(BSpline::knot[k + d + 3]-BSpline::knot[k + 3]);
	beta = (BSpline::knot[k + d + 4]-t)/(BSpline::knot[k + d + 4]-BSpline::knot[k + 4]);
	double temp =  (alpha * this -> BF(k, d - 1, t) + beta * this->BF(k+1, d-1, t));
	tuple<int, int, double>tupnew(k,d,t);
	BSpline::calculated.push_back(tupnew);
	BSpline::BFvalues.push_back(temp);
	return temp;

} //debugged


double BSpline::BF_integrate(int k, int d, double t){
	double sum = 0.0;
	for (int i = k; i < BSpline::N - d - 4; i++){ 
		double temp = this->BF(i,d+1,t);
		sum += (BSpline::knot[k+d+4] - BSpline::knot[k+3]) / (d+1.0) * temp;
	}
	return sum;
} //debugged



double BSpline::BF_integrate(int k, int d, double S, double T){
	return this->BF_integrate(k, d, T) - this->BF_integrate(k, d, S);
} //debugged


double BSpline::BF_derivative(int k, int d, double t){
	// t lies in the valid interval and d > 0
	double alpha, beta;
	alpha = double(d) / (BSpline::knot[k + d + 3]-BSpline::knot[k + 3]);
	beta = double(d) / (BSpline::knot[k + d + 4]-BSpline::knot[k + 4]);
	return (alpha * this->BF(k,d-1, t) - beta * this->BF(k+1,d-1, t));
}



double BSpline::BF_2integral(int k, int l, int d, double S, double T){
	int start, end;
	for (start = 0; BSpline::knot[start] <= S; start++);
	for (end = BSpline::N - 1; BSpline::knot[end] >= T; end--);
	double* w = new double[end - start+3];
	w[0] = S;
	w[end - start + 2] = T;
	for (int i = 1; i < end - start + 2; i++){
		w[i] = BSpline::knot[start + i - 1];
	}

	double Sum = this->BF_derivative(k, d, T) * this->BF_2derivative(l, d, T) - this->BF_derivative(k, d, S) * this->BF_2derivative(l, d, S);
	double Sum1 = 0.0;
	for (int j = 1; j <= end - start + 2; j++){
		Sum1 += this->BF_3derivative(l,d,w[j-1]) * (this->BF(k, d, w[j]) - this->BF(k,d,w[j-1]));
	}
	delete w;
	return Sum - Sum1;
}


double BSpline::BF_2derivative(int k, int d, double t){
	double alpha = double(d) / (BSpline::knot[k+d+3] - BSpline::knot[k+3]);
	double beta = double(d) / (BSpline::knot[k+d+4] - BSpline::knot[k+4]);
	return alpha * this->BF_derivative(k, d-1, t) - beta * this->BF_derivative(k+1, d-1,t);
}


double BSpline::BF_3derivative(int k, int d, double t){
	double alpha = double(d) / (BSpline::knot[k + d + 3] - BSpline::knot[k + 3]);
	double beta = double(d) / (BSpline::knot[k + d + 4] - BSpline::knot[k + 4]);
	return alpha * this->BF_2derivative(k, d-1, t) - beta * this->BF_2derivative(k+1, d-1,t);
}



double BSpline::DF(double t1, double t2, double f[]){
	double Sum = 0.0;
	for (int i = -3; i <= BSpline::N - 8; i++){
		Sum += -BF_integrate(i, 3, t1, t2) * f[i+3];
	}
	return exp(Sum);
}

double BSpline::DF(double t1, double t2){
	return DF(t1, t2, resultf);
}

double BSpline::LiborForward(double t1, double t2, double l[]){
	double Sum = 0.0;
	for (int i = -3; i <= BSpline::N - 8; i++){
		Sum += BF_integrate(i, 3, t1, t2) * l[i+3];
	}
	return 1.0 / (t2-t1) * (exp(Sum)-1.0); // delta means t2-t1
}

double BSpline::LiborForward(double t1, double t2){
	return LiborForward(t1, t2, resultl);
}

double BSpline::SwapRate(double t1, double t2, double f[], double l[]){
	//start == 0 means spot starting
	int numfloatcashflow = (float(t2) - float(t1) + 0.01) / 0.25;
	int numfixedcashflow = (float(t2) - float(t1) + 0.01) / 0.5;

	//Floating Leg
	double Sum1 = 0.0;
	for (int i = 0; i < numfloatcashflow; i++){
		Sum1 += this->LiborForward(t1 + i * 0.25, t1+ (i+1) * 0.25, l) * this->DF(t1, t1+(i+1) * 0.25, f) * 0.25;
	}

	//Fixed Leg
	double Sum2 = 0.0;
	for (int i = 0; i < numfixedcashflow; i++){
		Sum2 += this->DF(t1, t1 + (i+1) * 0.5, f) * 0.5;
	}
	/*
	//Accrual at the end, last payment
	if ((t2-t1)/0.25 != floor((t2-t1)/0.25)){
	Sum1 += this->LiborForward(t1 + numfloatcashflow * 0.25, t1 + (numfloatcashflow + 1) * 0.25, l) * this->DF(t1, t1 + (numfloatcashflow + 1) * 0.25, l) * (t2 - (t1 + numfloatcashflow * 0.25));
	Sum2 += this->DF(t1, t2, f) * (t2 - (t1 + 0.5 * numfixedcashflow)); 
	//cout << setprecision(12) << Sum2 << endl;
	}
	*/
	return Sum1 / Sum2;
}


double BSpline::SwapRate(double t1, double t2){
	return SwapRate(t1, t2, resultf, resultl);
}

double BSpline::OISForward(double t1, double t2, double OIS[]){
	return LiborForward(t1, t2, OIS);
}

double BSpline::OISForward(double t1, double t2){
	return OISForward(t1, t2, resultOIS);
}

double BSpline::BasisSwap(double t1, double t2, double f[], double l[], double OIS[]){  //OIS should be the same as f[] for current standard.
	//start == 0 means spot starting
	int numLiborcashflow = (float(t2) - float(t1) + 0.01) / 0.250000; //+0.01 for handing floating number decimals
	int numOIScashflow = (float(t2) - float(t1) + 0.01) / 0.250000;   // +0.01 for handing floating number decimals

	//Floating Leg
	double Sum1 = 0.0;
	for (int i = 0; i < numLiborcashflow; i++){
		Sum1 += 0.25 *  (LiborForward(t1, t1 + (i+1) * 0.25, l) - OISForward(t1, t1 + (i+1) * 0.25, OIS)) * DF(t1, t1 + (i+1) * 0.25, f);
	}

	//Fixed Leg
	double Sum2 = 0.0;
	for (int i = 0; i < numOIScashflow; i++){
		Sum2 += 0.25 * (DF(t1, t1 + (i+1) * 0.25, f));
		//cout << setprecision(12) << Sum2 << endl;
	}
	return Sum1 / Sum2;
}

double BSpline::BasisSwap(double t1, double t2){
	return BasisSwap(t1, t2, resultf, resultl, resultOIS);
}

double BSpline::Tikhonov_sum_squared(double f[], double l[]){
	double T0 = 0.0;
	double Tmax = 30.0;
	double Sum = 0;
	for (int i = -3; i <= BSpline::N - 8; i++){
		for (int j = -3; j <= BSpline::N - 8; j++){
			double temp = BF_2integral(i,j,3,T0, Tmax);
			Sum += f[i+3] * f[j+3] * temp;
			Sum += l[i+3] * l[j+3] * temp;
		}
	}
	return Sum;
}

double BSpline::Tikhonov_sum_squared(){
	return Tikhonov_sum_squared(resultf, resultl);
}


double BSpline::LiborSwapValue(double t1, double t2, double rate, double f[], double l[]){
	int numfloatcashflow = (float(t2) - float(t1) + 0.01) / 0.25;
	int numfixedcashflow = (float(t2) - float(t1) + 0.01) / 0.5;

	//Floating Leg
	double Sum1 = 0.0;
	for (int i = 0; i < numfloatcashflow; i++){
		Sum1 += LiborForward(t1 + i * 0.25, t1+ (i+1) * 0.25, l) * DF(t1, t1+(i+1) * 0.25, f) * 0.25;
	}

	//Fixed Leg
	double Sum2 = 0.0;
	for (int i = 0; i < numfixedcashflow; i++){
		Sum2 += rate * this->DF(t1, t1 + (i+1) * 0.5, f) * 0.5;
	}
	/*
	//Accrual at the end, last payment
	if ((t2-t1)/0.25 != floor((t2-t1)/0.25)){
	Sum1 += this->LiborForward(t1 + numfloatcashflow * 0.25, t1 + (numfloatcashflow + 1) * 0.25, l) * this->DF(t1, t1 + (numfloatcashflow + 1) * 0.25, l) * (t2 - (t1 + numfloatcashflow * 0.25));
	Sum2 += this->DF(t1, t2, f) * (t2 - (t1 + 0.5 * numfixedcashflow)); 
	//cout << setprecision(12) << Sum2 << endl;
	}
	*/
	return Sum2 - Sum1;
}

double BSpline::LiborSwapValue(double t1, double t2, double rate){
	return LiborSwapValue(t1, t2, rate, resultf, resultl);
}

double BSpline::BasisSwapValue(double t1, double t2, double rate, double f[], double l[], double OIS[]){
	int numLiborcashflow = (float(t2) - float(t1) + 0.01) / 0.250000; //+0.01 for handing floating number decimals
	int numOIScashflow = (float(t2) - float(t1) + 0.01) / 0.250000;   // +0.01 for handing floating number decimals

	//Floating Leg
	double Sum1 = 0.0;
	for (int i = 0; i < numLiborcashflow; i++){
		Sum1 += 0.25 *  (this->LiborForward(t1, t1 + (i+1) * 0.25, l) - this->OISForward(t1, t1 + (i+1) * 0.25, OIS)) * DF(t1, t1 + (i+1) * 0.25, f);
	}

	//Fixed Leg
	double Sum2 = 0.0;
	for (int i = 0; i < numOIScashflow; i++){
		Sum2 += 0.25 * (DF(t1, t1 + (i+1) * 0.25, f)) * rate;
		//cout << setprecision(12) << Sum2 << endl;
	}
	return Sum2 - Sum1;
}

double BSpline::BasisSwapValue(double t1, double t2, double rate){
	return BasisSwapValue(t1, t2, rate, resultf, resultl, resultOIS);
}