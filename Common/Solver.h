#ifndef REOSOLVER_SOLVER_H
#define REOSOLVER_SOLVER_H

#include <vector>
#include <assert.h>
#include "../Maths/Algebra.h"
#include "SpaceMesh.h"

using namespace riemann;



/************************************************************************************************/
/******* Abstract class for all following solvers (including exact Riemann solver) *************/
/************************************************************************************************/
class Solver
{

protected:

	double m_gamma;

	unsigned m_Nx;

	double m_hx;

	double m_delta_t = 0.0;



public:

	Solver(	double leftX, double rightX, unsigned Nx, double gamma = 5.0/3.0):
			m_Nx(Nx),
			m_hx((rightX - leftX) / Nx),
			m_gamma(gamma)
	{}


	virtual void Calculate(double time) = 0;

	virtual SpaceMesh GetSolution() const = 0;




};


/********************************************************************************************/
/************ Abstract class as a base for all following approximate solvers ***************/
/*******************************************************************************************/
class NumericalSolver: public Solver
{

protected:


	const std::vector<riemann::Vec3> m_W_Initial; // W = {ro, u, P}

//	double m_leftX, m_rightX;


	SpaceMesh m_currentState;


private:

	virtual riemann::Vec3 CalculateFlux(int i) const = 0; // Calculates F_{i + 1/2}



protected:


	template<typename T>
	Vec3 RungeKuttaTVDStep(Vec3 const& current_q, T const& f) const
	{

		auto const q_1 = current_q + m_delta_t * f(current_q);

		//	auto const q_2 = 3.0/4.0 * current_q + 1.0/4.0 * q_1 + 1.0/4.0 * m_delta_t * f(q_1);

		//	auto const next_q = 1.0/3.0 * current_q + 2.0/3.0 * q_2 + 2.0/3.0 * m_delta_t * f(q_2);

		//	return next_q;

		return q_1;

	}


	//!@brief Gets left and right values
	//!@param ro_l left density @param ro_l right density @param u_l left velocity @param u_r right velocity
	//!@param P_l left pressure P_r right pressure
	//!@param i Cell's number
	virtual void GetWLR(int i, double& ro_l, double& ro_r, double& u_l, double& u_r, double& P_l, double& P_r) const
	{
		ro_l = m_currentState[i].density; // == ro_i
		ro_r = m_currentState[i + 1].density; // == ro_i+1

		u_l = m_currentState[i].velocity; // == u_i
		u_r = m_currentState[i + 1].velocity; // == U_i+1

		P_l = m_currentState[i].pressure; // == p_i
		P_r = m_currentState[i + 1].pressure; // == p_i+1
	}



	//!@brief Calculates next time step according to Courant criteria
	void CalculateTimeStep()
	{
		auto lambda_max = 0.0;
		for(int i = 0; i < m_Nx; ++i)
		{
			auto c0 = std::sqrt(m_gamma * m_currentState[i].pressure / m_currentState[i].density);

			if(fabs(m_currentState[i].velocity + c0) > lambda_max)
				lambda_max = fabs(m_currentState[i].velocity + c0);
			if(fabs(m_currentState[i].velocity - c0) > lambda_max)
				lambda_max = fabs(m_currentState[i].velocity - c0);

		}

		double const sigma = 0.1;

		m_delta_t = sigma * m_hx / lambda_max;

//		double const phi = 1.0/3.0;

//		double const beta = 4.0;

//		return 2 * m_hx / (lambda_max * (5 - phi + (1 - phi) * beta));
	}


public:

	NumericalSolver(std::vector<riemann::Vec3> const& w0, double left = 0.0, double right = 1.0,
	unsigned Nx = 100, double gamma = 5.0/3.0): 	Solver(left, right, Nx, gamma),
												   	m_W_Initial(w0),
													m_currentState(Nx)
	{

	}



	//!@brief Does all calculations according to conservation law scheme
	//!@param time At this time we search the solution
	virtual void Calculate(double time);


	virtual SpaceMesh GetSolution() const
	{
		return m_currentState;
	}

};



/*
class ExactRiemannSolver: public Solver
{

private:

	double const Uright;
	double const Rright;
	double const Pright;
//	double const E01 = 0.0;

	double const Rleft;
	double const Pleft;
	double const Uleft;
//	double const E02 = 0;


//	double const Gam = 5.0/3.0;
	double const C2 = 2.04124145231932E+0000;
	double const V3  = 7.87760379795509E-0001;
	double const P3  = 2.51163357568276E+0000;
	double const DUV = 1.91890023213048E+0000;
	double const R3p = 1.69643057679328E+0000;
	double const R3m = 1.32318863224284E+0000;

private:

	double Gam1, Gam2, Gam3, Gam4;


	double m_leftX, m_rightX;
	SpaceMesh m_Solution;

public:

	ExactRiemannSolver(riemann::Vec3 const& w0Left,
					   riemann::Vec3 const& w0Right,
					   double leftX,
					   double rightX,
					   unsigned Nx,
					   double gamma): Solver(leftX, rightX, Nx, gamma),
									  m_leftX(leftX),
									  m_rightX(rightX),
									  Rright(w0Right[0]),
									  Uright(w0Right[1]),
									  Pright(w0Right[2]),
									  Rleft(w0Left[0]),
									  Uleft(w0Left[1]),
									  Pleft(w0Left[2]),
									  m_Solution(Nx)
	{
		Gam1 = (m_gamma + 1) / 2;
		Gam2 = (m_gamma-1)/(m_gamma+1);
		Gam3 = 2*m_gamma/(m_gamma-1);
		Gam4 = (m_gamma-1)/2;

		//C2 = std::sqrt(m_gamma * Pleft/ Rleft);
	}


	double C_(double X, double T);
	double Vt_ (double X, double T);
	double Rt_ (double X, double T);
	double Pt_ (double X, double T);
	double Et_ (double X, double T);


	virtual void Calculate(double time);

	virtual SpaceMesh GetSolution() const
	{
		return m_Solution;
	}
};
 */

#endif //REOSOLVER_SOLVER_H
