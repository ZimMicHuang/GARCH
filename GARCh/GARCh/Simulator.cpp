#include "Simulator.h"
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
#include <iostream>
#include <algorithm>

using namespace std;
Simulator::Simulator(Calibrator &c)
{
	Cali = &c;
	iterations = 10000;
	period = 21;
}

Simulator::Simulator(Calibrator &c, unsigned long int N)
{
	Cali = &c;
	iterations = N;
	period = 21;
}

Simulator::Simulator(Calibrator &c, unsigned int days) {
	Cali = &c;
	iterations = 10000;
	period = days;
}

Simulator::Simulator(Calibrator &c, unsigned int days, unsigned long int N) {
	Cali = &c;
	iterations = N;
	period = days;
}




vector<double> Simulator::singlePath(mt19937 &gen) {
	//Cali->print();
	double sigma = Cali->getSigma();
	double mu = Cali->getMu();
	//double mu = 0;
	double alpha = Cali->getAlpha();
	double beta = Cali->getBeta();
	double sigma0 = Cali->getSigma0();

	vector<double> epsilon = vector<double>(period+1);
	vector<double> vol = vector<double>(period+1);
	vector<double> price = vector<double>(period + 1);

	vol[0] = sigma0;
	//cout << "Sigma0: " << sigma0 << endl;
	price[0] = Cali->s0;
	//cout << "S0: " << Cali->s0 << endl;

	//create the epsilons
	//cout << "EPSILONS" << "\n" << endl;
	//normal_distribution<double> dist(0.0,1.0);
	student_t_distribution<double> dist(7.0);
	for (int i = 0; i < epsilon.size(); i++) {
		epsilon[i] = dist(gen);
		//cout << "Day " << i << " epsilon: " << epsilon[i] << endl;
	}
	//cout << "Noise generated" << endl;


	//create the sigmas
	//cout << "VOLATILITIES" << "\n" << endl;
	for (int j = 1; j < vol.size(); j++) {
		vol[j] = sqrt( (1 - alpha - beta)*sigma*sigma + 
			alpha * vol[j - 1] * epsilon[j - 1] * vol[j - 1] * epsilon[j - 1] + 
			beta * vol[j - 1] * vol[j - 1]);
		//cout << "Day " << j << " vol: " << vol[j] << endl;
	}
	//cout << "Volatilities generated" << endl;


	//prices
	//cout << "PRICES" << "\n" << endl;
	for (int k = 1; k < price.size(); k++) {
		double dS = price[k - 1] * (mu*(1.0 / 252.0) + vol[k] * sqrt(1.0 / 252.0)*epsilon[k]);
		price[k] = price[k - 1] + dS;
		//price[k] = vol[k] * epsilon[k];
		//cout << "Day " << k << " price: " << price[k] << endl;
	}
	//cout << "Price simulated" << endl;

	return price;
}

vector<vector<double>> Simulator::monteCarlo()
{
	vector<vector<double>> panel;
	mt19937 gen(time(NULL));
	for (int i = 0; i < iterations; i++) {
		//mt19937 gen(time(NULL)+rand()+i);
		panel.push_back(singlePath(gen));
		//cout << "Index " << i << " reached" << endl;
	}
	MCresults = panel;
	return panel;
}

vector<double> Simulator::getReturns()
{
	// note that these returns are returns over the period---may annualize them later
	vector<double> returns = {};
	for (int i = 0; i < MCresults.size(); i++) {
		double end = MCresults[i][MCresults[i].size() - 1];
		double beg = MCresults[i][0];
		returns.push_back(log(end / beg));
	}
	MCreturns = returns;
	sort(MCreturns.begin(), MCreturns.end());
	return returns;
}


void Simulator::report()
{
	/*
	for (unsigned int i = 0; i < MCreturns.size(); i++) {
		cout << "Path " << i + 1 << " : " << MCreturns[i] << endl;
	}
	*/
	cout << "Number of simulated paths: " << iterations << endl;
	cout << "--------------Summary Statistics (Annualized)--------------" << endl;
	cout << "Expected returns: " << Calibrator::mean(MCreturns) * 252.0 / (double)period << endl;
	cout << "Volatility: " << Calibrator::sd(MCreturns) * sqrt(252.0/(double)period) << endl;
	cout << "Kurtosis: " << Calibrator::kurt(MCreturns) << endl;
}




