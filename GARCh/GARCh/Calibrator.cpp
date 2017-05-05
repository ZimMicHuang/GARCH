#include "Calibrator.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <numeric>

using namespace std;

/////////////////////////////////////////////////////////////
/////////////// Basic Stats Functions ///////////////////////
/////////////////////////////////////////////////////////////

double Calibrator::mean(vector<double> vec) {
	double sum=0;
	for (unsigned int i=0; i<vec.size();i++) {
		sum+=vec.at(i);
	}
	return sum / vec.size();
}

double Calibrator::sd(vector<double> vec) {
	double sum = 0;
	double m = mean(vec);
	for (unsigned int i=0; i<vec.size(); i++) {
		sum+=(vec.at(i)-m)*(vec.at(i)-m);
	}
	return sqrt(sum / (vec.size()-1));
};

double Calibrator::cov(vector<double> vec, unsigned int lag) {
	double sum = 0;

	vector<double>::const_iterator first = vec.begin();
	vector<double>::const_iterator tmp1 = vec.end() - lag;
	vector<double>::const_iterator tmp2 = vec.begin() + lag;
	vector<double>::const_iterator last = vec.end();
	vector<double> vec1(first,tmp1);
	vector<double> vec2(tmp2,last);

	double m1 = mean(vec1);
	double m2 = mean(vec2);
	double s1 = sd(vec1);
	double s2 = sd(vec2);

	for (unsigned int i=0; i<vec1.size(); i++) {
		cout << vec1.at(i) <<endl;
		cout << vec2.at(i) << endl;
		sum+=(vec1.at(i)-m1)*(vec2.at(i)-m2);
	}

	return sum / (vec1.size()-1);
}

double Calibrator::autocorr(vector<double> vec, unsigned int lag) {
	double sum = 0;

	vector<double>::const_iterator first = vec.begin();
	vector<double>::const_iterator tmp1 = vec.end() - lag;
	vector<double>::const_iterator tmp2 = vec.begin() + lag;
	vector<double>::const_iterator last = vec.end();
	vector<double> vec1(first,tmp1);
	vector<double> vec2(tmp2,last);

	double m1 = mean(vec1);
	double m2 = mean(vec2);
	double s1 = sd(vec1);
	double s2 = sd(vec2);

	for (unsigned int i=0; i<vec1.size(); i++) {
		sum+=(vec1.at(i)-m1)*(vec2.at(i)-m2);
	}

	return sum / (vec1.size() - 1) / s1 / s2;
}

double Calibrator::kurt(vector<double> vec) {
	int n = vec.size();
	double m = mean(vec);
	long double sum=0;
	for (unsigned int i=0;i<vec.size();i++) {
		sum += pow((vec.at(i)-m),4.0);
	}
	return sum/(n-1)/pow( sd(vec),4.0 ) * n/(n-1);
}


/////////////////////////////////////////////////////////////
///////////// Functions for Tao Matching ////////////////////
/////////////////////////////////////////////////////////////

double Calibrator::taoEmpirical() {
	double sum = 0;
	vector<double> returnSquared = {};
	for (int i = 0; i< returns.size(); i++) {
		returnSquared.push_back(returns.at(i)*returns.at(i));
	}

	for (int i = 1; i< 120; i++) {
		sum+=autocorr(returnSquared,i);
	}

	taoEmp = sum;
	return sum;
}

double Calibrator::taoTheoretical() {
	// this function is actually called after both matching, i.e. alpha and beta, are done
	double k = kurtEmp;
	double tao = (  alpha * ( 1-(alpha+beta)*(alpha+beta)+alpha*(alpha+beta) )  ) 
		/  (  1-(alpha+beta)*(alpha+beta)+alpha*alpha  ) 
		*	(1-pow(alpha+beta,120)) / (1-alpha-beta);
	taoThr = tao;
	return tao;
}

double Calibrator::taoEquation(double x) {
	double k = kurtEmp;
	double a = (sqrt(  (k-3) * (k*(3-2*x*x)-3) )- x*(k-3))  /  3 / (k-1)  ;
	double tao = (  a * ( 1-(a+x)*(a+x)+a*(a+x) )  ) /  (  1-(a+x)*(a+x)+a*a  ) * (1-pow(a+x,120)) / (1-x-a);
	return tao-taoEmp;
}

double Calibrator::bisect(double left, double right, double tol) {
	// this function relies on a correct taoEquation function --- which 
	// requires a pre-calculated empirical kurtosis and a pre-calculated emprical tao

	unsigned int N = 1;
	while (N <= 5000) {
		double mid = (left + right) / 2;
		if (taoEquation(mid) == 0 || (right - left) / 2 < tol) {
			return mid;
		}

		if (taoEquation(mid) * taoEquation(left) > 0) {
			left = mid;
			N++;
		}
		else {
			right = mid;
			N++;
		}
	}
	return -1; // failure code -1
}



void Calibrator::solveBeta() {
	double res = bisect(0.00001, 1.00000, 0.000001);
	if (res != -1) {
		this->beta = res;
	}
	else {
		this->beta = 0;
	}
}

/////////////////////////////////////////////////////////////
////////// Functions that solves other params ///////////////
/////////////////////////////////////////////////////////////

void Calibrator::solveMuSigma() {
	mu = mean(returns) * 252;
	sigma = sd(returns) * sqrt(252);
}


void Calibrator::solveAlpha() {
	double k = kurtEmp;
	alpha = (sqrt((k - 3) * (k*(3 - 2 * beta*beta) - 3)) - beta*(k - 3)) /
		3 / (k - 1);
}

void Calibrator::solveSigma0() {
	vector<double>::const_iterator tmp1 = returns.begin();
	vector<double>::const_iterator tmp2 = returns.begin()+21;
	vector<double> tmp3(tmp1,tmp2);
	sigma0 = sd(tmp3)*sqrt(252);
}

void Calibrator::solveEpsilon0() {
	double innov = ( returns.back() - mu ) * ( returns.back() - mean(returns) );
	epsilon0 = sqrt(innov/sigma0/sigma0);
}

void Calibrator::solveKurt() {
	kurtEmp = kurt(returns);
}



/////////////////////////////////////////////////////////////
/////////////// Public Member Functions /////////////////////
/////////////////////////////////////////////////////////////

Calibrator::Calibrator(vector<double> ret, double s0) {
	this->s0 = s0;
	returns = ret;
}

void Calibrator::solve() {
	// this specific order must be followed due to implicit dependencies
	solveMuSigma(); 
	solveKurt();
	taoEmpirical();
	solveBeta();
	solveAlpha();
	taoTheoretical();
	kurtThreoretical();
	solveSigma0();
	solveEpsilon0();
}