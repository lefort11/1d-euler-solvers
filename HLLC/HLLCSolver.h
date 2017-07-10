#ifndef REOSOLVER_HLLCSOLVER_H
#define REOSOLVER_HLLCSOLVER_H

#include "../Common/Solver.h"

/*********************************************************/
/******* HLLC gas dynamics one dimensional solver*********/
/********************************************************/

class HLLCSolver: public NumericalSolver
{

public:

	HLLCSolver(std::vector<riemann::Vec3> const& w0, double left = 0.0,
			   double right = 1.0, unsigned Nx = 100, double gamma = 5.0/3.0):
			NumericalSolver(w0, left, right, Nx, gamma)
	{}

private:

	//!@brief Calulates flux for the i+1/2 cell's boundary
	//!@param i Cell's number
	virtual riemann::Vec3 CalculateFlux(int i) const;


};


#endif //REOSOLVER_HLLCSOLVER_H
