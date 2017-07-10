#include "REOSolver.h"




//!@brief Calculates F_{i+1/2} - flux in the i-th cell;
riemann::Vec3 REOSolver::CalculateFlux(int i) const
{

	//Getting i-th and (i+1)-th node data
	// -----------------------------------//
	double ro_l, ro_r, u_l, u_r, P_l, P_r;

	GetWLR(i, ro_l, ro_r, u_l, u_r, P_l, P_r); // Getting cell data

	auto const eps_l = P_l / ( ro_l * (m_gamma - 1.0) );
	auto const eps_r = P_r / ( ro_r * (m_gamma - 1.0) );

	auto const h_l = eps_l + P_l / ro_l + u_l * u_l / 2.0;
	auto const h_r = eps_r + P_r / ro_r + u_r * u_r / 2.0;

	auto const c_l = sqrt( m_gamma * P_l / ro_l );
	auto const c_r = sqrt( m_gamma * P_r / ro_r );
	//--------------------------------------//



	//Calulating first order approximation term (R-E)
	//------------------------------------------------------------------------//
	riemann::Vec3 F_i(ro_l * u_l, ro_l * u_l * u_l + P_l, ro_l * u_l * h_l);

	riemann::Vec3 F_i_1(ro_r * u_r, ro_r * u_r * u_r + P_r, ro_r * u_r * h_r);

	auto const lambda_i_1 = u_l - c_l;
	auto const lambda_i_plus_1_3 = u_r + c_r;


	std::array<riemann::Vec3, 3> r_star;
	riemann::Vec3 delta_s;
	riemann::Vec3 lambda_star;

	CalculateStarComponents(i, ro_l, ro_r, u_l, u_r, P_l, P_r, delta_s, r_star, lambda_star);

	std::array<double, 3> lambdas = {
			std::fmin(lambda_star[0], lambda_i_1),
			lambda_star[1],
			std::fmax(lambda_star[2], lambda_i_plus_1_3)
	}; // enthropy fix

	/*std::array<double, 3> lambdas = {
			lambda_star[0],
			lambda_star[1],
			lambda_star[2]
	}; */


	riemann::Vec3 flux = 0.5 * (F_i + F_i_1);
	for(int m = 0; m < 3; ++m)
	{
		flux = flux - 0.5 * std::fabs(lambdas[m]) * delta_s[m] * r_star[m];
	}
	//---------------------------------------------------------------------------------//


	//High order approximation term (Osher scheme) calculation
	//----------------------------------------------------------------------------//
	riemann::Vec3 OsherTerm = {0.0, 0.0, 0.0};

	auto const F_plus_i_minus_1div2 = CalculateFPlus(i-1, true);
	auto const F_plus_i_plus_1div2 = CalculateFPlus(i, true);


	auto const F_minus_i_plus_1div2 = CalculateFPlus(i, false);
	auto const F_minus_i_plus_3div2 = CalculateFPlus(i + 1, false);


	static double const phi = 1.0/3.0;
	static double const beta = 4.0;

	for(int m = 0; m < 3; ++m)
	{
		OsherTerm = OsherTerm + (1 + phi) / 4 * riemann::minmod(F_plus_i_plus_1div2[m], beta * F_plus_i_minus_1div2[m])
					+ (1 - phi) / 4 * riemann::minmod(beta * F_plus_i_plus_1div2[m], F_plus_i_minus_1div2[m])
					- (1 + phi) / 4 * riemann::minmod(F_minus_i_plus_1div2[m],
												 beta * F_minus_i_plus_3div2[m])
					- (1 - phi) / 4 * riemann::minmod(beta * F_minus_i_plus_1div2[m],
												 F_minus_i_plus_3div2[m]);

	}

	//-------------------------------------------------------------------------------------------------//


	return flux + OsherTerm;
//	return flux;

}


//!@brief calculates F^{m+}_{i+1/2} if plus and F^{m-}_{i+1/2} if !plus - 3rd order approximation term
std::array<riemann::Vec3, 3> REOSolver::CalculateFPlus(int i, bool plus) const
{

	//Getting i-th and (i+1)-th node data
	double ro_l, ro_r, u_l, u_r, P_l, P_r;
	GetWLR(i, ro_l, ro_r, u_l, u_r, P_l, P_r);


	std::array<riemann::Vec3, 3> r_star;
	riemann::Vec3 delta_s;
	riemann::Vec3 lambda_star;
	CalculateStarComponents(i, ro_l, ro_r, u_l, u_r, P_l, P_r, delta_s, r_star, lambda_star);

	riemann::Vec3 lambdas;

	if(plus)
	{
		lambdas = {
				std::fmax(lambda_star[0], 0.0),
				std::fmax(lambda_star[1], 0.0),
				std::fmax(lambda_star[2], 0.0)
		};
	}
	else
	{
		lambdas = {
				std::fmin(lambda_star[0], 0.0),
				std::fmin(lambda_star[1], 0.0),
				std::fmin(lambda_star[2], 0.0)
		};

	}

	std::array<riemann::Vec3, 3> F;
	for(int m = 0; m < 3; ++m)
	{
		F[m] = lambdas[m] * delta_s[m] * r_star[m];
	}

	return F;

}






//!@brief Calculates r_star, delta_s, and lambdas for q* values
void REOSolver::CalculateStarComponents(int i,
									 double ro_l, double ro_r,
									 double u_l, double u_r,
									 double P_l, double P_r,
									 riemann::Vec3 &delta_s, std::array<riemann::Vec3, 3> &r_star,
									 riemann::Vec3 &lambda_star) const
{

	auto const eps_l = P_l / ( ro_l * (m_gamma - 1.0) );
	auto const eps_r = P_r / ( ro_r * (m_gamma - 1.0) );

	auto const h_l = eps_l + P_l / ro_l + u_l * u_l / 2.0;
	auto const h_r = eps_r + P_r / ro_r + u_r * u_r / 2.0;


	//calculating q* components
	auto const ro_star = sqrt(ro_l * ro_r);
	auto const u_star = ( sqrt(ro_l) * u_l + sqrt(ro_r) * u_r ) / ( sqrt(ro_l) + sqrt(ro_r) );
	auto const h_star = ( sqrt(ro_l) * h_l + sqrt(ro_r) * h_r ) / ( sqrt(ro_l) + sqrt(ro_r) );
	auto const c_star = sqrt( (m_gamma - 1.0) * (h_star - u_star * u_star / 2.0) );


	lambda_star = {
			u_star - c_star,
			u_star,
			u_star + c_star
	};


	r_star = {
			riemann::Vec3(1.0, u_star - c_star, h_star - u_star * c_star),

			riemann::Vec3(2.0, 2.0 * u_star, u_star * u_star),

			riemann::Vec3(1.0, u_star + c_star, h_star + u_star*c_star)
	};

	delta_s = {
			0.5/(c_star * c_star) * ( (P_r - P_l) - ro_star * c_star * (u_r - u_l) ),
			0.5/(c_star * c_star) * ( c_star * c_star * (ro_r - ro_l) - (P_r - P_l) ),
			0.5/(c_star * c_star) * ( (P_r - P_l) + ro_star * c_star * (u_r - u_l) ),


	};
}


