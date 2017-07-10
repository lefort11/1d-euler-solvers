#include "Algebra.h"
#include "assert.h"

double riemann::sgn(double x)
{
	if( x > 0 )
		return 1;
	if( x < 0 )
		return -1;
	return 0;
}

double riemann::minmod(double x, double y)
{
	return 0.5 * (sgn(x) + sgn(y)) * std::min(fabs(x), fabs(y));
}


riemann::Vec3 riemann::minmod(riemann::Vec3 const &x, riemann::Vec3 const &y)
{
	return riemann::Vec3(riemann::minmod(x[0], y[0]), riemann::minmod(x[1], y[1]), riemann::minmod(x[2], y[2]));
}


double riemann::dot(Vec3 const &first, Vec3 const &second)
{
	return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}


const riemann::Vec3 riemann::operator+(Vec3 const &left, Vec3 const &right)
{
	return Vec3(left[0] + right[0], left[1] + right[1], left[2] + right[2]);
}

const riemann::Vec3 riemann::operator-(Vec3 const &left, Vec3 const &right)
{
	return Vec3(left[0] - right[0], left[1] - right[1], left[2] - right[2]);
}

const riemann::Vec3 riemann::operator*(Vec3 const &vec, double alpha)
{
	return Vec3(alpha * vec[0], alpha * vec[1], alpha * vec[2]);
}


const riemann::Vec3 riemann::operator*(double alpha, Vec3 const &vec)
{
	return Vec3(alpha * vec[0], alpha * vec[1], alpha * vec[2]);
}

riemann::Vec3& riemann::operator+=(riemann::Vec3 &left, riemann::Vec3 const &right)
{
	left[0] += right[0];
	left[1] += right[1];
	left[2] += right[2];
	return left;
}

riemann::Matrix3x3& riemann::operator*=(Matrix3x3 &mat, double alpha)
{
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			mat[i][j] *= alpha;
	return mat;
}

double riemann::sqr(double x)
{
	return x*x;
}


riemann::Vec3 riemann::sqr(riemann::Vec3 const& vec)
{
	riemann::Vec3 result;
	for(int i = 0; i < 3; ++i)
	{
		result[i] = riemann::sqr(vec[i]);
	}
	return result;
}


const riemann::Vec3 riemann::operator*(riemann::Matrix3x3 const &mat, riemann::Vec3 const &vec3)
{
	riemann::Vec3 result;
	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < 3; ++j)
			result[i] += mat[i][j] * vec3[j];
	}
	return result;
}