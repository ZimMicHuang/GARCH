#include "dataProc.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;
dataProc::dataProc(string filePath)
{
	this->filePath = filePath;
	returns = {};
	s0 = 0;
}

void dataProc::getData()
{
	vector<double> raw = {};
	ifstream myFile;
	myFile.open(filePath);
	string line;
	while (getline(myFile, line))
	{
		string value;
		stringstream ss(line);
		double x = 0;
		while (getline(ss, value, ','))
		{
			//cout << value << " ";
			x = stod(value);
			raw.push_back(x);
		}
	}
	myFile.close();

	//cout << "we are here" << endl;

	this->s0 = raw[0];

	for (unsigned int i = 0; i < raw.size() - 1; i++) {
		double ret = log(raw[i] / raw[i + 1]);
		returns.push_back(ret);
	}
}
