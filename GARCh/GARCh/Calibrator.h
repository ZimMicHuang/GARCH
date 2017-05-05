#include <vector>
#include <iostream>
//#include "Simulator.h"

using namespace std;

class Calibrator{
	friend class Simulator;

private:
	//basic stats functions
	static double mean(vector<double> vec);
	static double sd(vector<double> vec);
	static double cov(vector<double> vec, unsigned int lag);
	static double autocorr(vector<double> vec,unsigned int lag);
	static double kurt(vector<double> vec);


	//helper functions that solve the clustering time equation
	double taoEquation(double alpha);
	double bisect(double left, double right, double tol);
	double taoEmpirical();

	//updates the theoretical tao
	double taoTheoretical();


	//these solve functions updates the attribute of the class and return the solved parameter
	void solveAlpha();
	void solveBeta();
	void solveMuSigma();
	void solveSigma0();
	void solveEpsilon0();
	void solveKurt();

	inline void kurtThreoretical() {
		this->kurtThr = 3 / (1 - 2 * alpha*alpha / (1 - (alpha + beta)*(alpha + beta)));
	}

	//attributes; the parameters for the GARCH(1,1)

	// from 5/2/2017, mu, sigma, and sigma0 are annualized. The factor is 1/252. 
	vector<double> returns;
	double alpha;
	double beta;
	double mu;
	double sigma;
	double sigma0;
	double epsilon0;


	// matched results; for model checking purposes
	double kurtEmp;
	double kurtThr;
	double taoEmp;
	double taoThr;


	
public:
	double s0;
	//constructors
	Calibrator(vector<double> returns,double s0);

	//getters
	inline const double getMu() { return mu; }
	inline const double getSigma() { return sigma; }
	inline const double getAlpha() { return alpha; }
	inline const double getBeta() { return beta; }
	inline const double getSigma0() { return sigma0; }


	//the main function
	void solve();//the main calibration function

	inline void print() {
		std::cout << "alpha: " << alpha << endl;
		std::cout << "beta: " << beta << endl;
		std::cout << "mu: " << mu << endl;
		std::cout << "sigma: " << sigma << endl;
		std::cout << "---------check---------" << endl;
		std::cout << "kurt of empirical data: " << kurtEmp<< endl;
		std::cout << "kurt of theoretical model: " << kurtThr<< endl;
		std::cout << "tao of empirical data: " << taoEmp<< endl;
		std::cout << "tao of theoretical model: " << taoThr<< endl;

		std::cout << "starting price is: " << s0 << endl;
	}
	
};