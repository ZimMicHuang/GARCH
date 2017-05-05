#pragma once
#include <string>
#include <vector>
using namespace std;
class dataProc
{
public:
	string filePath;
	vector<double> returns;
	double s0;

	dataProc(string filePath);

	void getData();

};

