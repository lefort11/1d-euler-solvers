#ifndef REOSOLVER_WENOSOLVER_H
#define REOSOLVER_WENOSOLVER_H

#include "../Common/Solver.h"
#include "../REO/REOSolver.h"
#include "../HLLC/HLLCSolver.h"

#include <type_traits>


/**************************************************************************************************************/
/******* Gas dynamics WENO Solver based on WENO reconstruction and one of approximate Riemann solvers ********/
/************************************************************************************************************/


// T is a Riemann Approximate Solver class
template <class T>
class WENOSolver: public T
{

private:

	//!@brief Makes the WENO component-wise reconstruction for q = {ro, ro*u, ro*E}
	//!@param i Cell's number @param bias =true if stencil is biased to right
	riemann::Vec3 Reconstruct(int i, bool biasRight) const //component-wise reconstruction
	{
		static riemann::Vec3 const gammas = {
				1.0/10.0,
				3.0/5.0,
				3.0/10.0
		};

		static auto const b_1 = 13.0/12.0;
		static auto const b_2 = 1.0/4.0;

		static auto const eps = 1e-6;

		static riemann::Vec3 const a_1 = {
				1.0/3.0,
				-7.0/6.0,
				11.0/6.0
		};

		static riemann::Vec3 const a_2 = {
				-1.0/6.0,
				5.0/6.0,
				1.0/3.0
		};

		static riemann::Vec3 const a_3 = {
				1.0/3.0,
				5.0/6.0,
				-1.0/6.0
		};

		std::vector<riemann::Vec3> q_s(5);

		riemann::Vec3 result;

		if(!biasRight)
		{
			for(int k = 0; k < 5; ++k)
			{

				auto const ro = T::m_currentState[i + k - 2].density;
				auto const u = T::m_currentState[i + k - 2].velocity;
				auto const P = T::m_currentState[i + k - 2].pressure;

				auto const enrg = P / ((T::m_gamma - 1) * ro);
				auto const E = enrg + 0.5 * u * u;



				q_s[k] = {

						ro,

						ro * u,

						ro * E
				};
			}

			riemann::Vec3 const q_1 = a_1[0] * q_s[0] + a_1[1] * q_s[1] + a_1[2] * q_s[2];

			riemann::Vec3 const q_2 = a_2[0] * q_s[1] + a_2[1] * q_s[2] + a_2[2] * q_s[3];

			riemann::Vec3 const q_3 = a_3[0] * q_s[2] + a_3[1] * q_s[3] + a_3[2] * q_s[4];


			for( int m = 0; m < 3; ++m )
			{
				auto const beta_1 = b_1 * riemann::sqr(q_s[0][m] - 2 * q_s[1][m] + q_s[2][m]) +
									b_2 * riemann::sqr(q_s[0][m] - 4 * q_s[1][m] + 3 * q_s[2][m]);

				auto const beta_2 = b_1 * riemann::sqr(q_s[1][m] - 2 * q_s[2][m] + q_s[3][m]) +
									b_2 * riemann::sqr(q_s[1][m] - q_s[3][m]);

				auto const beta_3 = b_1 * riemann::sqr(q_s[2][m] - 2 * q_s[3][m] + q_s[4][m]) +
									b_2 * riemann::sqr(3 * q_s[2][m] - 4 * q_s[3][m] + q_s[4][m]);



				auto const w_t_1 = gammas[0] / riemann::sqr(eps + beta_1);
				auto const w_t_2 = gammas[1] / riemann::sqr(eps + beta_2);
				auto const w_t_3 = gammas[2] / riemann::sqr(eps + beta_3);

				auto const w_summ = w_t_1 + w_t_2 + w_t_3;

				auto const w_1 = w_t_1 / w_summ;
				auto const w_2 = w_t_2 / w_summ;
				auto const w_3 = w_t_3 / w_summ;

				result[m] = w_1 * q_1[m] + w_2 * q_2[m] + w_3 * q_3[m];

			}
		}
		else
		{
			for(int k = 0; k < 5; ++k)
			{

				auto const ro = T::m_currentState[i + k - 1].density;
				auto const u = T::m_currentState[i + k - 1].velocity;
				auto const P = T::m_currentState[i + k - 1].pressure;

				auto const enrg = P / ((T::m_gamma - 1) * ro);
				auto const E = enrg + 0.5 * u * u;
				q_s[k] = {

						ro,

						ro * u,

						ro * E

				};
			}

			riemann::Vec3 const q_1 = a_2[0] * q_s[0] + a_2[1] * q_s[1] + a_2[2] * q_s[2];

			riemann::Vec3 const q_2 = a_3[0] * q_s[1] + a_3[1] * q_s[2] + a_3[2] * q_s[3];

			riemann::Vec3 const q_3 = a_1[2] * q_s[2] + a_1[1] * q_s[3] + a_1[0] * q_s[4];


			for( int m = 0; m < 3; ++m )
			{
				auto const beta_1 = b_1 * riemann::sqr(q_s[0][m] - 2 * q_s[1][m] + q_s[2][m]) +
									b_2 * riemann::sqr(q_s[0][m] - 4 * q_s[1][m] + 3 * q_s[2][m]);

				auto const beta_2 = b_1 * riemann::sqr(q_s[1][m] - 2 * q_s[2][m] + q_s[3][m]) +
									b_2 * riemann::sqr(q_s[1][m] - q_s[3][m]);

				auto const beta_3 = b_1 * riemann::sqr(q_s[2][m] - 2 * q_s[3][m] + q_s[4][m]) +
									b_2 * riemann::sqr(3 * q_s[2][m] - 4 * q_s[3][m] + q_s[4][m]);



				auto const w_t_1 = gammas[2] / riemann::sqr(eps + beta_1);
				auto const w_t_2 = gammas[1] / riemann::sqr(eps + beta_2);
				auto const w_t_3 = gammas[0] / riemann::sqr(eps + beta_3);

				auto const w_summ = w_t_1 + w_t_2 + w_t_3;

				auto const w_1 = w_t_1 / w_summ;
				auto const w_2 = w_t_2 / w_summ;
				auto const w_3 = w_t_3 / w_summ;

				result[m] = w_1 * q_1[m] + w_2 * q_2[m] + w_3 * q_3[m];

			}
		}


		return result;
	}

	//!@brief Makes the WENO characteristic-wise reconstruction for q = {ro, ro*u, ro*E}
	//!@param i Cell's number @param bias =true if stencil is biased to right
	riemann::Vec3 ReconstructChar(int i, bool biasRight) const //characteristic-wise reconstructuion
	{
		riemann::Matrix3x3 R, Rinv;


		FillRandRinv(i, R, Rinv);

		static riemann::Vec3 const gammas = {
				1.0/10.0,
				3.0/5.0,
				3.0/10.0
		};

		static auto const b_1 = 13.0/12.0;
		static auto const b_2 = 1.0/4.0;

		static auto const eps = 1e-6;

		static riemann::Vec3 const a_1 = {
				1.0/3.0,
				-7.0/6.0,
				11.0/6.0
		};

		static riemann::Vec3 const a_2 = {
				-1.0/6.0,
				5.0/6.0,
				1.0/3.0
		};

		static riemann::Vec3 const a_3 = {
				1.0/3.0,
				5.0/6.0,
				-1.0/6.0
		};

		std::vector<riemann::Vec3> q_s(5);

		riemann::Vec3 result;

		if(!biasRight)
		{
			for(int k = 0; k < 5; ++k)
			{

				auto const ro = T::m_currentState[i + k - 2].density;
				auto const u = T::m_currentState[i + k - 2].velocity;
				auto const P = T::m_currentState[i + k - 2].pressure;

				auto const enrg = P / ((T::m_gamma - 1) * ro);
				auto const E = enrg + 0.5 * u * u;



				q_s[k] = {

						ro,

						ro * u,

						ro * E
				};

				q_s[k] = Rinv * q_s[k];
			}

			riemann::Vec3 const q_1 = a_1[0] * q_s[0] + a_1[1] * q_s[1] + a_1[2] * q_s[2];

			riemann::Vec3 const q_2 = a_2[0] * q_s[1] + a_2[1] * q_s[2] + a_2[2] * q_s[3];

			riemann::Vec3 const q_3 = a_3[0] * q_s[2] + a_3[1] * q_s[3] + a_3[2] * q_s[4];


			for( int m = 0; m < 3; ++m )
			{
				auto const beta_1 = b_1 * riemann::sqr(q_s[0][m] - 2 * q_s[1][m] + q_s[2][m]) +
									b_2 * riemann::sqr(q_s[0][m] - 4 * q_s[1][m] + 3 * q_s[2][m]);

				auto const beta_2 = b_1 * riemann::sqr(q_s[1][m] - 2 * q_s[2][m] + q_s[3][m]) +
									b_2 * riemann::sqr(q_s[1][m] - q_s[3][m]);

				auto const beta_3 = b_1 * riemann::sqr(q_s[2][m] - 2 * q_s[3][m] + q_s[4][m]) +
									b_2 * riemann::sqr(3 * q_s[2][m] - 4 * q_s[3][m] + q_s[4][m]);



				//auto const w_t_1 = gammas[0] / riemann::sqr(eps + beta_1);
				//auto const w_t_2 = gammas[1] / riemann::sqr(eps + beta_2);
				//auto const w_t_3 = gammas[2] / riemann::sqr(eps + beta_3);

				auto const w_t_1 = gammas[0] * (1 + beta_2) * (1 + beta_3);
				auto const w_t_2 = gammas[1] * (1 + beta_1) * (1 + beta_3);
				auto const w_t_3 = gammas[2] * (1 + beta_1) * (1 + beta_2);



				auto const w_summ = w_t_1 + w_t_2 + w_t_3;

				auto const w_1 = w_t_1 / w_summ;
				auto const w_2 = w_t_2 / w_summ;
				auto const w_3 = w_t_3 / w_summ;

				result[m] = w_1 * q_1[m] + w_2 * q_2[m] + w_3 * q_3[m];

			}
		}
		else
		{
			for(int k = 0; k < 5; ++k)
			{

				auto const ro = T::m_currentState[i + k - 1].density;
				auto const u = T::m_currentState[i + k - 1].velocity;
				auto const P = T::m_currentState[i + k - 1].pressure;

				auto const enrg = P / ((T::m_gamma - 1) * ro);
				auto const E = enrg + 0.5 * u * u;
				q_s[k] = {

						ro,

						ro * u,

						ro * E

				};

				q_s[k] = Rinv * q_s[k];
			}

			riemann::Vec3 const q_1 = a_2[0] * q_s[0] + a_2[1] * q_s[1] + a_2[2] * q_s[2];

			riemann::Vec3 const q_2 = a_3[0] * q_s[1] + a_3[1] * q_s[2] + a_3[2] * q_s[3];

			riemann::Vec3 const q_3 = a_1[2] * q_s[2] + a_1[1] * q_s[3] + a_1[0] * q_s[4];


			for( int m = 0; m < 3; ++m )
			{
				auto const beta_1 = b_1 * riemann::sqr(q_s[0][m] - 2 * q_s[1][m] + q_s[2][m]) +
									b_2 * riemann::sqr(q_s[0][m] - 4 * q_s[1][m] + 3 * q_s[2][m]);

				auto const beta_2 = b_1 * riemann::sqr(q_s[1][m] - 2 * q_s[2][m] + q_s[3][m]) +
									b_2 * riemann::sqr(q_s[1][m] - q_s[3][m]);

				auto const beta_3 = b_1 * riemann::sqr(q_s[2][m] - 2 * q_s[3][m] + q_s[4][m]) +
									b_2 * riemann::sqr(3 * q_s[2][m] - 4 * q_s[3][m] + q_s[4][m]);



//				auto const w_t_1 = gammas[2] / riemann::sqr(eps + beta_1);
//				auto const w_t_2 = gammas[1] / riemann::sqr(eps + beta_2);
//				auto const w_t_3 = gammas[0] / riemann::sqr(eps + beta_3);

				auto const w_t_1 = gammas[0] * (1 + beta_2) * (1 + beta_3);
				auto const w_t_2 = gammas[1] * (1 + beta_1) * (1 + beta_3);
				auto const w_t_3 = gammas[2] * (1 + beta_1) * (1 + beta_2);

				auto const w_summ = w_t_1 + w_t_2 + w_t_3;

				auto const w_1 = w_t_1 / w_summ;
				auto const w_2 = w_t_2 / w_summ;
				auto const w_3 = w_t_3 / w_summ;

				result[m] = w_1 * q_1[m] + w_2 * q_2[m] + w_3 * q_3[m];

			}
		}

		result = R * result;

		return result;
	}


	//!@brief Redefined GetWLR method that is used to calculate Fluxes.
	//!@brief Now q_l and q_r vectors are composed of reconstructed vectors
	virtual void GetWLR(int i, double& ro_l, double& ro_r, double& u_l, double& u_r, double& P_l, double& P_r) const
	{
		auto const q_l = ReconstructChar(i, false);
		auto const q_r = ReconstructChar(i, true);


		ro_l = q_l[0];
		ro_r = q_r[0];

		u_l = q_l[1] / ro_l;
		u_r = q_r[1] / ro_r;

		auto const E_l = q_l[2] / ro_l;
		auto E_r = q_r[2] / ro_r;

		auto const eps_l = E_l - 0.5 * u_l * u_l;
		auto const eps_r = E_r - 0.5 * u_r * u_r;

		P_l = (T::m_gamma - 1) * eps_l * ro_l;
		P_r = (T::m_gamma - 1) * eps_r * ro_r;
	}

	//!@brief Fills R and R^(-1) matrices, where R is the matrix composed of eigenvectors of Jacobian
	void FillRandRinv(int i, riemann::Matrix3x3& R, riemann::Matrix3x3& Rinv) const
	{
		double ro_l = 0.0, ro_r = 0.0, u_l = 0.0, u_r = 0.0, P_l = 0.0, P_r = 0.0;

		T::GetWLR(i, ro_l, ro_r, u_l, u_r, P_l, P_r);

		auto const eps_l = P_l / ((T::m_gamma - 1) * ro_l);
		auto const eps_r = P_r / ((T::m_gamma - 1) * ro_r);

		auto const E_l = eps_l + u_l * u_l * 0.5;
		auto const E_r = eps_r + u_r * u_r * 0.5;

		auto const h_l = E_l + P_l / ro_l;
		auto const h_r = E_r + P_r / ro_r;

/*	riemann::Vec3 q_star = {
			(ro_l + ro_r) * 0.5,
			(ro_l * u_l + ro_r * u_r) * 0.5,
			(ro_l * E_l + ro_r * E_r) * 0.5
	};

	auto ro_star = q_star[0];
	auto u_star = q_star[1] / ro_star;
	auto E_star = q_star[2] / ro_star;

	auto eps_star = E_star - 0.5 * u_star * u_star;

	auto P_star = (m_gamma - 1) * eps_star * ro_star;

	auto c_star = std::sqrt(m_gamma * P_star / ro_star);

	auto h_star = E_star + P_star/ro_star; */

//	auto ro_star = sqrt(ro_l * ro_r);
		auto const u_star = ( sqrt(ro_l) * u_l + sqrt(ro_r) * u_r ) / ( sqrt(ro_l) + sqrt(ro_r) );
		auto const h_star = ( sqrt(ro_l) * h_l + sqrt(ro_r) * h_r ) / ( sqrt(ro_l) + sqrt(ro_r) );
		auto const c_star = sqrt( (T::m_gamma - 1.0) * (h_star - u_star * u_star / 2.0) );

		//filling R here
		R[0][0] = 1;
		R[0][1] = 2;
		R[0][2] = 1;
		R[1][0] = u_star - c_star;
		R[1][1] = 2 * u_star;
		R[1][2] = u_star + c_star;
		R[2][0] = h_star - u_star * c_star;
		R[2][1] = u_star * u_star;
		R[2][2] = h_star + u_star * c_star;

		//filling Rinv here

		Rinv[0][0] = 0.5 * u_star * u_star + u_star * c_star / (T::m_gamma - 1);
		Rinv[0][1] = -u_star - c_star / (T::m_gamma - 1);
		Rinv[0][2] = 1;

		Rinv[1][0] = c_star * c_star / (T::m_gamma - 1) - 0.5 * u_star * u_star;
		Rinv[1][1] = u_star;
		Rinv[1][2] = -1;

		Rinv[2][0] = 0.5 * u_star * u_star - u_star * c_star / (T::m_gamma - 1);
		Rinv[2][1] = -u_star + c_star / (T::m_gamma - 1);
		Rinv[2][2] = 1;

		Rinv *= (T::m_gamma - 1) / (2 * c_star * c_star);


	}


public:

	WENOSolver(std::vector<riemann::Vec3> const& w0, double left = 0.0, double right = 1.0,
			   unsigned Nx = 100, double gamma = 5.0/3.0): T(w0, left, right, Nx, gamma)
	{
		static_assert(std::is_base_of<NumericalSolver, T>::value, "T must be derived from NumericalSolver class");
	}



};



#endif //REOSOLVER_WENOSOLVER_H
