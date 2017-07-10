#ifndef REOSOLVER_EXACTRIEMANNSOLVER_H
#define REOSOLVER_EXACTRIEMANNSOLVER_H


#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "../Common/Solver.h"


/*****************************************************************************/
/*********** Exact Riemann Solver by E.F.Toro ********************************/
/*****************************************************************************/



class ExactRiemannSolver: public Solver
{

	double g1, g2, g3, g4, g5, g6, g7, g8;

	double                // density, velocity, pressure, speed of sound
			rl, ul, pl, cl,   // in left region
			rr, ur, pr, cr;   // in right region

	double m_x0;


	SpaceMesh m_solution;



public:
	ExactRiemannSolver(riemann::Vec3 const& w0left,
					   riemann::Vec3 const& w0right,
					   double leftX,
					   double rightX,
					   unsigned Nx,
					   double gamma,
					   double x0): Solver(leftX, rightX, Nx, gamma),
								   m_solution(Nx),
								   m_x0(x0)
	{
		// compute gamma related constants
		g1 = (m_gamma - 1.0)/(2.0*m_gamma);
		g2 = (m_gamma + 1.0)/(2.0*m_gamma);
		g3 = 2.0*m_gamma/(m_gamma - 1.0);
		g4 = 2.0/(m_gamma - 1.0);
		g5 = 2.0/(m_gamma + 1.0);
		g6 = (m_gamma - 1.0)/(m_gamma + 1.0);
		g7 = (m_gamma - 1.0)/2.0;
		g8 = m_gamma - 1.0;



		rl = w0left[0];
		rr = w0right[0];

		ul = w0left[1];
		ur = w0right[1];

		pl = w0left[2];
		pr = w0right[2];

	}


	virtual void Calculate(double time)
	{

		// compute sound speeds
		cl = sqrt(m_gamma*pl/rl);
		cr = sqrt(m_gamma*pr/rr);

		// the pressure positivity condition is tested for
		if (g4*(cl+cr) <= (ur-ul)) {

			std::cerr << "the initial data is such that vacuum is generated"
				 << "\nstopping program" << std::endl;
			exit(1);
		}


		double pm, um;

		// exact solution for pressure and velocity in star region is found
		starpu(pm, um, 1.0);


		for(int i = 0; i < m_Nx; ++i)
		{
			auto xpos = (i + 0.5) * m_hx;
			auto s = (xpos - m_x0)/time;


			double rs, us, ps;
			sample(pm, um, s, rs, us, ps);

			m_solution[i].density = rs;
			m_solution[i].pressure = ps;
			m_solution[i].velocity = us;
		}
	}

	SpaceMesh GetSolution() const
	{
		return m_solution;
	}


private:


	void starpu(
			double &p,
			double &u,
			const double pscale)
	{
		// purpose: to compute the solution for pressure and
		//          velocity in the Star Region

		const int nriter = 20;
		const double tolpre = 1.0e-6;
		double change, fl, fld, fr, frd, pold, pstart, udiff;

		// guessed value pstart is computed
		guessp(pstart);
		pold = pstart;
		udiff = ur - ul;

		std::cout << "----------------------------------------\n"
			 << "   Iteration number     Change\n"
			 << "----------------------------------------" << std::endl;

		int i = 1;
		for ( ; i <= nriter; i++) {
			prefun(fl, fld, pold, rl, pl, cl);
			prefun(fr, frd, pold, rr, pr, cr);
			p = pold - (fl + fr + udiff)/(fld + frd);
			change = 2.0*fabs((p - pold)/(p + pold));
			std::cout << '\t' << i <<  "\t\t" << change << std::endl;
			if (change <= tolpre)
				break;
			if (p < 0.0)
				p = tolpre;
			pold = p;
		}
		if (i > nriter)
			std::cout << "divergence in Newton-Raphson iteration" << std::endl;

		// compute velocity in star region
		u = 0.5*(ul + ur + fr - fl);
		std::cout << "----------------------------------------\n"
			 << "     Pressure           Velocity\n"
			 << "----------------------------------------\n"
			 << "     " << p/pscale << "\t\t" <<  u << '\n'
			 << "----------------------------------------" << std::endl;
	}

	void sample(
			const double pm,
			const double um,
			const double s,
			double &d,
			double &u,
			double &p)
	{
		// purpose: to sample the solution throughout the wave
		//          pattern. Pressure pm and velocity um in the
		//          star region are known. Sampling is performed
		//          in terms of the 'speed' s = x/t. Sampled
		//          values are d, u, p

		double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

		if (s <= um) {
			// sampling point lies to the left of the contact discontinuity
			if (pm <= pl) {
				// left rarefaction
				shl = ul - cl;
				if (s <= shl) {
					// sampled point is left data state
					d = rl;
					u = ul;
					p = pl;
				} else {
					cml = cl*pow(pm/pl, g1);
					stl = um - cml;
					if (s > stl) {
						// sampled point is star left state
						d = rl*pow(pm/pl, 1.0/m_gamma);
						u = um;
						p = pm;
					} else {
						// sampled point is inside left fan
						u = g5*(cl + g7*ul + s);
						c = g5*(cl + g7*(ul - s));
						d = rl*pow(c/cl, g4);
						p = pl*pow(c/cl, g3);
					}
				}
			} else {
				// left shock
				pml = pm/pl;
				sl = ul - cl*sqrt(g2*pml + g1);
				if (s <= sl) {
					// sampled point is left data state
					d = rl;
					u = ul;
					p = pl;
				} else {
					// sampled point is star left state
					d = rl*(pml + g6)/(pml*g6 + 1.0);
					u = um;
					p = pm;
				}
			}
		} else {
			// sampling point lies to the right of the contact discontinuity
			if (pm > pr) {
				// right shock
				pmr = pm/pr;
				sr  = ur + cr*sqrt(g2*pmr + g1);
				if (s >= sr) {
					// sampled point is right data state
					d = rr;
					u = ur;
					p = pr;
				} else {
					// sampled point is star right state
					d = rr*(pmr + g6)/(pmr*g6 + 1.0);
					u = um;
					p = pm;
				}
			} else {
				// right rarefaction
				shr = ur + cr;
				if (s >= shr) {
					// sampled point is right data state
					d = rr;
					u = ur;
					p = pr;
				} else {
					cmr = cr*pow(pm/pr, g1);
					str = um + cmr;
					if (s <= str) {
						// sampled point is star right state
						d = rr*pow(pm/pr, 1.0/m_gamma);
						u = um;
						p = pm;
					} else {
						// sampled point is inside left fan
						u = g5*(-cr + g7*ur + s);
						c = g5*(cr - g7*(ur - s));
						d = rr*pow(c/cr, g4);
						p = pr*pow(c/cr, g3);
					}
				}
			}
		}
	}


	void guessp(double &pm)
	{
		// purpose: to provide a guessed value for pressure
		//          pm in the Star Region. The choice is made
		//          according to adaptive Riemann solver using
		//          the PVRS, TRRS and TSRS approximate
		//          Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1

		double cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr,
				qmax, quser, um;

		quser = 2.0;

		// compute guess pressure from PVRS Riemann solver
		cup = 0.25*(rl + rr)*(cl + cr);
		ppv = 0.5*(pl + pr) + 0.5*(ul - ur)*cup;
		ppv = std::max(0.0, ppv);
		pmin = std::min(pl,  pr);
		pmax = std::max(pl,  pr);
		qmax = pmax/pmin;

		if (qmax <= quser && (pmin <= ppv && ppv <= pmax))
			pm = ppv;     // select PVRS Riemann solver
		else {
			if (ppv < pmin) {
				// select Two-Rarefaction Riemann solver
				pq = pow(pl/pr, g1);
				um = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
				ptl = 1.0 + g7*(ul - um)/cl;
				ptr = 1.0 + g7*(um - ur)/cr;
				pm = 0.5*(pow(pl*ptl, g3) + pow(pr*ptr, g3));
			} else {
				// select Two-Shock Riemann solver with PVRS as estimate
				gel = sqrt((g5/rl)/(g6*pl + ppv));
				ger = sqrt((g5/rr)/(g6*pr + ppv));
				pm = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
			}
		}
	}

	void prefun(
			double &f,
			double &fd,
			double &p,
			double &dk,
			double &pk,
			double &ck)
	{
		// purpose: to evaluate the pressure functions
		//          fl and fr in exact Riemann solver
		//          and their first derivatives

		double ak, bk, pratio, qrt;

		if (p <= pk) {
			// rarefaction wave
			pratio = p/pk;
			f = g4*ck*(pow(pratio, g1) - 1.0);
			fd = (1.0/(dk*ck))*pow(pratio, -g2);
		} else {
			//  shock wave
			ak = g5/dk;
			bk = g6*pk;
			qrt = sqrt(ak/(bk + p));
			f = (p - pk)*qrt;
			fd = (1.0 - 0.5*(p - pk)/(bk + p))*qrt;
		}
	}



};


#endif //REOSOLVER_EXACTRIEMANNSOLVER_H
