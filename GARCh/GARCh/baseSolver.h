#pragma once

class NonlinearSolver
{
public: 
	double (*f) (double);
	double tol;
	long int MAX_ITERATION;
	inline NonlinearSolver () { }
	inline NonlinearSolver (double (*f) (double)) { }
	virtual double solve()=0;
};

