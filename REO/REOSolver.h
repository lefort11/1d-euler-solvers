#ifndef REOSOLVER_REOSOLVER_H
#define REOSOLVER_REOSOLVER_H

#include "../Common/Solver.h"


/**************************************************************************/
/******     Gas dynamics one dimensional solver              **************/
/******     based on Roe-Einfeldt-Osher Scheme               **************/
/**************************************************************************/

class REOSolver: public NumericalSolver
{


private:

	//!@brief calculates third-order Osher's flux
	//!@param i - cell's number
	//!@param plus = true if calculating F^{+} flux according to formula, = false => F^{-}
	std::array<riemann::Vec3, 3> CalculateFPlus(int i, bool plus) const; //Calculates F^{m+}_{i+1/2}


	//!@brief Calculates Roe averages, fills parameter variables
	//!@param i cell's number
	//!@param ro_l left density @param ro_r right density @param u_l left velocity @param u_r right velocity
	//!@param P_l left pressure @param P_r right pressure
	//!@param delta_s Difference of characteristic vector (Riemann invariants)( == l * (q_i+1 - q_i), l - left eigenvector)
	//!@param r_star right eigenvector at the discontinuity, @param lambdas eigenvalues
	void CalculateStarComponents(int i,
								 double ro_l,
								 double ro_r,
								 double u_l,
								 double u_r,
								 double P_l,
								 double P_r,
								 riemann::Vec3& delta_s, std::array<riemann::Vec3, 3>& r_star, riemann::Vec3& lambdas) const;


	//!@brief Calculates Roe-Einfeldt-Osher flux
	//!@param i cell's number
	virtual riemann::Vec3 CalculateFlux(int i) const;


public:

	REOSolver(std::vector<riemann::Vec3> const& w0, double left = 0.0, double right = 1.0,
			  unsigned Nx = 100, double gamma = 5.0/3.0): NumericalSolver(w0, left, right, Nx, gamma)
	{}



};


#endif //REOSOLVER_REOSOLVER_H
