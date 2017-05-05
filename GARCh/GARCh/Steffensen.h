#include "baseSolver.h"
#include <math.h>
#include <iostream>
#include <assert.h>
#include <vector>


//tested 4/11/2017 works fine with the following test code:
/*
#include "Steffensen.h"
#include <cmath>

double PlusOneCube(double x) { return pow((x+1),3); }

int main() 
{
	SteffensenSolver SteffSolver(-0.5,PlusOneCube);
	SteffSolver.tol = 0.000005;
	SteffSolver.solve();
	SteffSolver.print();
	getchar();
	return 0;
}
*/

/*
class SteffensenSolver:public NonlinearSolver
{
	friend class Calibrator;
private:
	double Xcurrent, Xprevious, X0;
	long int iteration;
	Calibrator *Cali;
	double(*f) (double, Calibrator);

public:
	inline SteffensenSolver(double X0guess, double (*Fsteff) (double , Calibrator), Calibrator c) 
	{
		X0 = X0guess;
		f = Fsteff;
		Cali = &c;
		Xcurrent = Xprevious = X0;
		iteration = 0;
		MAX_ITERATION = 50000;
	}

	inline double Gsteff( double x , double (*f)(double, Calibrator)) 
	{
		double gx;
		gx = (f(x+f(x,*Cali), *Cali) - f(x, *Cali))/f(x, *Cali);
		return gx;
	}
	inline double solve() 
	{
		assert (iteration < MAX_ITERATION );
		if (fabs(f(Xcurrent,*Cali)) <= tol) { return Xcurrent; }
		else 
		{	using namespace std;
			Xcurrent = Xprevious - f(Xprevious, *Cali)/Gsteff(X0,f);
			cout << Xcurrent << endl;
			Xprevious = Xcurrent; 
			iteration++; 
			cout<<iteration<<endl;
			return solve(); }
	}
	inline void print() 
	{
		using namespace std;
		cout << "Data on Steffenson's Method requested" << endl;
		cout << "The root is: " << Xcurrent << endl;
		cout << "The numnber of iteration is: " << iteration << endl;
	}
};

*/

#include "baseSolver.h"
#include <math.h>
#include <iostream>
#include <assert.h>

class SteffensenSolver :public NonlinearSolver
{
private:
	double Xcurrent, Xprevious, X0;
	long int iteration;
public:
	inline SteffensenSolver(double X0guess, double(*Fsteff)(double))
	{
		X0 = X0guess;
		f = Fsteff;
		Xcurrent = Xprevious = X0;
		iteration = 0;
		MAX_ITERATION = 1000;
	}
	inline double Gsteff(double x, double(*f)(double))
	{
		double gx;
		gx = (f(x + f(x)) - f(x)) / f(x);
		return gx;
	}
	inline double solve()
	{
		assert(iteration < MAX_ITERATION);
		if (fabs(f(Xcurrent)) <= tol) { return Xcurrent; }
		else
		{
			using namespace std;
			Xcurrent = Xprevious - f(Xprevious) / Gsteff(X0, f);
			cout << Xcurrent << endl;
			Xprevious = Xcurrent;
			iteration++;
			cout << iteration << endl;
			return solve();
		}
	}
	inline void print()
	{
		using namespace std;
		cout << "Data on Steffenson's Method requested" << endl;
		cout << "The root is: " << Xcurrent << endl;
		cout << "The numnber of iteration is: " << iteration << endl;
	}
};