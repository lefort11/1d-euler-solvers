#include <iostream>
#include <fstream>
#include "REO/REOSolver.h"

#include "ExactRiemannSolver/ExactRiemannSolver.h"

#include "HLLC/HLLCSolver.h"

#include "WENO/WENOSolver.h"

namespace riemann
{
	Vec3 NormL2(SpaceMesh func, double h)
	{
		riemann::Vec3 norm = {0,0,0};
		for(int i = 0; i < func.size(); ++i)
		{
			norm[0] += sqr(func[i].density);
			norm[1] += sqr(func[i].velocity);
			norm[2] += sqr(func[i].pressure);
		}

		norm = norm * h;

		norm[0] = std::sqrt(norm[0]);
		norm[1] = std::sqrt(norm[1]);
		norm[2] = std::sqrt(norm[2]);

		return norm;

	}


	Vec3 NormC(SpaceMesh func, double h, double to)
	{
		riemann::Vec3 normC = {0.0, 0.0, 0.0};
		for(int i = 0; (i + 0.5) * h <= to; ++i)
		{
			if(std::fabs(func[i].density) > normC[0])
				normC[0] = std::fabs(func[i].density);
			if(std::fabs(func[i].velocity) > normC[1])
				normC[1] = std::fabs(func[i].velocity);
			if(std::fabs(func[i].pressure) > normC[2])
				normC[2] = std::fabs(func[i].pressure);
		}
		return normC;
	}

}


class RiemannTest
{

public:
	riemann::Vec3 const left;
	riemann::Vec3 const right;
	double const x0;

	double const leftX, rightX;

	unsigned const Nx = 200;
	double const hx = (rightX - leftX)/Nx;

	double const gamma;

public:
	RiemannTest(riemann::Vec3 const& left, riemann::Vec3 const& right, double leftX, double rightX,
				double const x0, double const gamma):
			left(left),
			right(right),
			leftX(leftX),
			rightX(rightX),
			x0(x0),
			gamma(gamma)
	{}

	std::vector<riemann::Vec3> GetW0() const
	{
		std::vector<riemann::Vec3> w0(Nx);

		for(int i = 0; i < Nx; ++i)
		{
			if(leftX + (i + 0.5) * hx < x0)
				w0[i] = left;
			else
				w0[i] = right;
		}
		return w0;

	}

};


int main(int argc, char** argv)
{
	RiemannTest test2({1.0, 0.75, 1.0}, {0.125, 0.0, 0.1}, 0.0, 1.0, 0.3, 1.4);
	RiemannTest test1({2.0, 0.0, 5.0}, {1.0, 0.0, 1.0}, 0.0, 1.0, 0.5, 5.0/3.0);

	RiemannTest test3({1.0, -2.0, 0.4}, {1.0, 2.0, 0.4}, 0.0, 1.0, 0.5, 5.0/3.0);

	RiemannTest test4({8.0, 0.0, 480.0}, {1.0, 0.0, 1.0}, 0.0, 1.0, 0.5, 5.0/3.0);

	double time = 0.177;


	RiemannTest* pTest = &test1;

	REOSolver reoSolver(pTest->GetW0(), pTest->leftX, pTest->rightX, pTest->Nx, pTest->gamma);
//	reoSolver.Calculate(time);

	WENOSolver<HLLCSolver> wenoSolver(pTest->GetW0(), pTest->leftX, pTest->rightX, pTest->Nx, pTest->gamma);
	wenoSolver.Calculate(time);

	auto solution = wenoSolver.GetSolution();

	std::ofstream densityResultsFile, pressureResultsFile, velocityResultsFile;

	densityResultsFile.open("results/density_numerical_1.txt");
	pressureResultsFile.open("results/pressure_numerical_1.txt");
	velocityResultsFile.open("results/velocity_numerical_1.txt");

	for(int i = 0; i < pTest->Nx; ++i)
	{
		densityResultsFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << solution[i].density << std::endl;
		velocityResultsFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << solution[i].velocity << std::endl;
		pressureResultsFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << solution[i].pressure << std::endl;

	}

	densityResultsFile.close();


	ExactRiemannSolver exactRiemannSolver(pTest->left, pTest->right, pTest->leftX, pTest->rightX,
			pTest->Nx, pTest->gamma, pTest->x0);
	exactRiemannSolver.Calculate(time);
	auto exactSolution = exactRiemannSolver.GetSolution();

	std::ofstream exactDensityResultsFile, exactVelocityResultsFile, exactPressureResultsFile;


	exactDensityResultsFile.open("results/density_exact_1.txt");
	exactVelocityResultsFile.open("results/velocity_exact_1.txt");
	exactPressureResultsFile.open("results/pressure_exact_1.txt");

	for(int i = 0; i < pTest->Nx; ++i)
	{
		exactDensityResultsFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << exactSolution[i].density << std::endl;
		exactVelocityResultsFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << exactSolution[i].velocity << std::endl;
		exactPressureResultsFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << exactSolution[i].pressure << std::endl;
	}


	exactDensityResultsFile.close();
	riemann::Vec3 normL2 = riemann::NormL2(solution - exactSolution, pTest->hx);

	auto normC = riemann::NormC(solution - exactSolution, pTest->hx, 0.5);


	std::cout << "Density norm: L2 " << normL2[0] << " C " << normC[0] << std::endl;
	std::cout << "Velocity norm: L2 " << normL2[1] << " C " << normC[1] << std::endl;
	std::cout << "Pressure norm: L2 " << normL2[2] << " C " << normC[2] << std::endl;

	riemann::Vec3 amplitudeContact = {1.75 - 1.35, 10000, 10000};
	riemann::Vec3 amplitudeShock = { 1.75, 0.79, 2.49 };

	riemann::Vec3 pointContact = {}, pointShock = {};
	for(int k = 100; k < 150; ++k)
	{
		if(std::fabs(solution[k].density - exactSolution[k].density) > (amplitudeContact[0] / 140))
			pointContact[0]++;


	}
	for(int k = 150; k < 200; ++k)
	{
		if(std::fabs(solution[k].density - exactSolution[k].density) > (amplitudeShock[0] / 400))
			pointShock[0]++;
		if(std::fabs(solution[k].velocity - exactSolution[k].velocity) > (amplitudeShock[1] / 400))
			pointShock[1]++;
		if(std::fabs(solution[k].pressure - exactSolution[k].pressure) > (amplitudeShock[2] / 400))
			pointShock[2]++;
	}

	std::cout << pointContact[0] << " " << pointContact[1] << " " << pointContact[2] << std::endl;
	std::cout << pointShock[0] << " " << pointShock[1] << " " << pointShock[2] << std::endl;


/*
	HLLCSolver hllcSolver(pTest->GetW0(), pTest->leftX, pTest->rightX, pTest->Nx, pTest->gamma);
	hllcSolver.Calculate(time);

	auto hllcSolution = hllcSolver.GetSolution();

	std::ofstream hllcDensityResFile, hllcVelocityResFile, hllcPressureResFile;
	hllcDensityResFile.open("results/density_hllc.txt");
	hllcVelocityResFile.open("results/velocity_hllc.txt");
	hllcPressureResFile.open("results/pressure_hllc.txt");


	for(int i = 0; i < pTest->Nx; ++i)
	{
		hllcDensityResFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << hllcSolution[i].density << std::endl;
		hllcVelocityResFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << hllcSolution[i].velocity << std::endl;
		hllcPressureResFile << pTest->leftX + (i + 0.5) * pTest->hx << " " << hllcSolution[i].pressure << std::endl;

	}
*/

	return 0;

}