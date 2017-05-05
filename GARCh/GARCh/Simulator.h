#pragma once
#include <iostream>
#include <vector>
#include <thread>
#include <random>
#include <ctime>
#include "Calibrator.h"

using namespace std;
class Simulator {
	//friend class Calibrator;
private: // may change later on
	unsigned long int iterations;
	unsigned int period;
	Calibrator *Cali;
	vector<vector<double>> MCresults;
	vector<double> MCreturns;
	vector<double> returnsSorted;

public:
	//Constructor
	Simulator(Calibrator &c);
	Simulator(Calibrator &c, unsigned long int N);
	Simulator(Calibrator &c, unsigned int days);
	Simulator(Calibrator &c, unsigned int days, unsigned long int N);

	//member functions 
private:
	//simulation
	vector<double> singlePath(mt19937 &gen);
public:
	vector<vector<double>> monteCarlo();

	//summary
	vector<double> getReturns();
	//double getVAR(int percentage);
	inline double getVAR(double percentage) {return MCreturns.at(MCreturns.size() * (1 - percentage));}
	void report();
	



};


