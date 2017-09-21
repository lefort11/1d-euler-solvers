#include "HLLCSolver.h"

riemann::Vec3 HLLCSolver::CalculateFlux(int i) const
{

	// ********* pressure estimate **************/
	double ro_l, ro_r, u_l, u_r, P_l, P_r;

	GetWLR(i, ro_l, ro_r, u_l, u_r, P_l, P_r);

	auto const c_l = sqrt(m_gamma * P_l/ro_l);
	auto const c_r = sqrt(m_gamma * P_r/ro_r);

	auto const ro_bar = 0.5 * (ro_l + ro_r);
	auto const c_bar = 0.5 * (c_l + c_r);

	auto const P_pvrs = 0.5 * (P_l + P_r) - 0.5 * (u_r - u_l) * ro_bar * c_bar;

	auto const P_star = std::fmax(0, P_pvrs);
	//********************************************//

	// ********** wave speed estimates ************** //
	auto const q_l = P_star <= P_l ? 1 : sqrt( 1 + (m_gamma + 1) / (2 * m_gamma) * (P_star / P_l - 1));
	auto const q_r = P_star <= P_r ? 1 : sqrt( 1 + (m_gamma + 1) / (2 * m_gamma) * (P_star / P_r - 1));

	auto const S_l = u_l - c_l * q_l;
	auto const S_r = u_r + c_r * q_r;

	auto const S_star = ( P_r - P_l + ro_l * u_l * (S_l - u_l) - ro_r * u_r * (S_r - u_r) ) /
			( ro_l*(S_l - u_l) - ro_r * (S_r - u_r) );

	// ********************************************* //

	// *************** HLLC flux ******************** //


	auto const eps_l = P_l / ((m_gamma - 1) * ro_l);
	auto const eps_r = P_r / ((m_gamma - 1) * ro_r);

	auto const E_l = (eps_l + u_l * u_l / 2);
	auto const E_r = (eps_r + u_r * u_r / 2);

	auto const H_l = eps_l + P_l / ro_l + u_l * u_l / 2;
	auto const H_r = eps_r + P_r / ro_r + u_r * u_r / 2;

	riemann::Vec3 const F_l = {
			ro_l * u_l,
			ro_l * u_l * u_l + P_l,
			ro_l * u_l * H_l
	};

	riemann::Vec3 const F_r = {
			ro_r * u_r,
			ro_r * u_r * u_r + P_r,
			ro_r * u_r * H_r
	};


	if(0 <= S_l)
	{
		return F_l;
	}
	else if((S_l <= 0) && (0 <= S_star))
	{
		riemann::Vec3 U_star_l = {
				1.0,
				S_star,
				E_l + (S_star - u_l) * (S_star + P_l / (ro_l * (S_l - u_l)))
		};

		U_star_l = ro_l * (S_l - u_l) / (S_l - S_star) * U_star_l;

		riemann::Vec3 U_l = {
				ro_l,
				ro_l * u_l,
				ro_l * E_l
		};

		riemann::Vec3 HLLCflux = F_l + S_l * (U_star_l - U_l);

		return HLLCflux;


	}
	else if( (S_star <= 0) && (0 <= S_r))
	{
		riemann::Vec3 U_star_r = {
				1.0,
				S_star,
				E_r + (S_star - u_r) * (S_star + P_r / (ro_r * (S_r - u_r)))
		};


		U_star_r = ro_r * (S_r - u_r) / (S_r - S_star) * U_star_r;

		riemann::Vec3 U_r = {
				ro_r,
				ro_r * u_r,
				ro_r * E_r
		};


		riemann::Vec3 HLLCflux = F_r + S_r * (U_star_r - U_r);

		return HLLCflux;

	}
	else if(0 >= S_r)
	{
		return F_r;
	}

	return {0.0,0.0,0.0};
}
