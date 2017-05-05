//#include "Calibrator.h"
#include "Simulator.h"
#include "dataProc.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <sstream>
#include <ctime>


using namespace std;
/*
vector<double> IStoVec(istream& is) {
	vector<double> test_rets = {};
	string line;
	getline(is, line);

	stringstream lineStream = stringstream(line);
	string cell;

	while (getline(lineStream, cell, ",")) {
		
	}

}
*/




int main() 
{
	
	
	//step1: read returns data from csv file,
	// and put the input into vector<double> format; also, store s0.
	
	string path = "G:/WORK/QF465project/GARCh/test2.csv";
	dataProc d(path);
	d.getData();


	//step2: initiate Calibrator and calibrate
	vector<double> test_rets = d.returns;
	double s0 = d.s0;

	
	Calibrator test_calibrator = Calibrator(test_rets,s0);
	test_calibrator.solve();
	test_calibrator.print();


	// step 3: initialize simulator (third constructor) and simulate
	Simulator s1 = Simulator(test_calibrator, 21, 10000);
	s1.monteCarlo();
	s1.getReturns();
	s1.report();
	cout << "95% Var:  " << s1.getVAR(0.95) << endl;
	cout << "Annualized 95% Var:" << s1.getVAR(0.95)*252.0/21.0 << endl;

	// step 4: do stuff with the returns data
	getchar();
	return 0;
}