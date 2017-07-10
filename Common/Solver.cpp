#include "Solver.h"





//!@brief does all stuff
void NumericalSolver::Calculate(double time)
{


	//Initializing mesh
	for(unsigned i = 0; i < m_Nx; ++i)
	{
		m_currentState[i].density = m_W_Initial[i][0];
		m_currentState[i].velocity = m_W_Initial[i][1];
		m_currentState[i].pressure = m_W_Initial[i][2];

	}



	//initializing q^0 vector
	std::vector<riemann::Vec3> currentQ(m_Nx);

	for(unsigned i = 0; i < m_Nx; ++i)
	{
		//internal energy
		double epsilon = m_currentState[i].pressure / (m_currentState[i].density * (m_gamma - 1.0));

		//full energy
		double E = epsilon + m_currentState[i].velocity * m_currentState[i].velocity / 2.0;

		currentQ[i] = riemann::Vec3(m_currentState[i].density, m_currentState[i].density * m_currentState[i].velocity,
								m_currentState[i].density * E);
	}

	// ~ q^{n+1}
	std::vector<riemann::Vec3> nextQ(m_Nx, {0.0, 0.0, 0.0});


	auto currentTime = 0.0;
	auto tau = CalculateTimeStep();

	if(tau > time)
		tau = time;
	currentTime += tau;

	while( (tau != 0) )
	{

		//calculating q^{n+1} vectors for each x_i
#pragma omp parallel for
		for(int i = 0; i < m_Nx; ++i)
		{

			riemann::Vec3 leftFlux = CalculateFlux(i-1);
			riemann::Vec3 rightFlux = CalculateFlux(i);

			nextQ[i] = currentQ[i] - tau / m_hx * (rightFlux - leftFlux);
		}

		//filling solution here
		for(unsigned i = 0; i < m_Nx; ++i)
		{

			m_currentState[i].density = nextQ[i][0];
			m_currentState[i].velocity = nextQ[i][1] / m_currentState[i].density;

			double E = nextQ[i][2] / m_currentState[i].density;
			double epsilon = E - m_currentState[i].velocity * m_currentState[i].velocity / 2.0;


			m_currentState[i].pressure = (m_gamma - 1.0) * epsilon * m_currentState[i].density;
		}

		tau = CalculateTimeStep();
		if(currentTime + tau> time)
			tau = time - currentTime;
		currentTime += tau;


		currentQ = nextQ;
	}
}


/*
double ExactRiemannSolver::C_(double X, double T)
{
	return C2-Gam2*(X/T+C2);
}

double ExactRiemannSolver::Et_(double X, double T)
{
	return Pt_(X,T)/((m_gamma-1)*Rt_(X,T))+0.5*riemann::sqr(Vt_(X,T));
}


double ExactRiemannSolver::Pt_(double X, double T)
{
	if (X<-C2*T)
		return Pleft;
	else if (X< (Gam1*V3-C2)*T)
		return Pleft*std::exp(std::log(C_(X,T)/C2)*Gam3);
	else if (X< DUV*T)
		return P3;
	else  return Pright;
}

double ExactRiemannSolver::Rt_(double X, double T)
{
	if (X<-C2*T)
		return Rleft;
	else if (X < (Gam1*V3-C2)*T)
		return Rleft*std::exp(std::log(C_(X,T)/C2)/Gam4);
	else if (X < V3*T)
		return R3m;
	else if (X< DUV*T)
		return R3p;
	else
		return Rright;
}

double ExactRiemannSolver::Vt_(double X, double T)
{
	if (X<-C2*T)
		return Uleft;
	else if (X < (Gam1*V3-C2)*T)
		return (X+C2*T)/(Gam1*T);
	else if (X < DUV*T)
		return V3;
	else
		return Uright;
}

void ExactRiemannSolver::Calculate(double time)
{
	for(int i = 0; i < m_Nx; ++i)
	{
		m_Solution[i].density = Rt_(m_leftX + (i + 0.5) * m_hx, time);
		m_Solution[i].velocity = Vt_(m_leftX + (i + 0.5) * m_hx, time);
		m_Solution[i].pressure = Pt_(m_leftX + (i + 0.5) * m_hx, time);
	}

}
 */










