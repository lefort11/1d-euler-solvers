#ifndef REOSOLVER_ALGEBRA_H
#define REOSOLVER_ALGEBRA_H

#pragma once

#include <array>
#include <cmath>

#include "../Common/SpaceMesh.h"

namespace riemann
{
	class Vec3: public std::array<double, 3>
	{



	public:

		Vec3(double x = 0.0, double y = 0.0, double z = 0.0): std::array<double, 3>{{x, y, z}}
		{}


		friend const Vec3 operator+(Vec3 const& left, Vec3 const& right);

		friend const Vec3 operator-(Vec3 const& left, Vec3 const& right);

		friend const Vec3 operator*(double alpha, Vec3 const& vec);

		friend const Vec3 operator*(Vec3 const& vec, double alpha);

		friend riemann::Vec3& operator+=(riemann::Vec3& first, riemann::Vec3 const& right);



	};


	class Matrix3x3
	{
		class Row
		{
			std::array<double, 3> elements;
		public:
			double& operator[] (const int j)
			{
				if(j >= 3)
					throw 1;
				return elements[j];
			}
			double const& operator[] (const int j) const
			{
				if(j >= 3)
					throw 1;
				return elements[j];
			}
		};
		std::array<Row, 3> rows;

	public:
		Row& operator[] (const int i)
		{
			if(i >= 3)
				throw 1;
			return rows[i];
		}
		Row const& operator[] (const int i) const
		{
			if(i >= 3)
				throw 1;
			return rows[i];
		}


		friend const riemann::Vec3 operator*(Matrix3x3 const& mat, riemann::Vec3 const& vec3);

		friend Matrix3x3& operator*=(Matrix3x3& mat, double alpha);

	};


	Matrix3x3& operator*=(Matrix3x3& mat, double alpha);

	const riemann::Vec3 operator*(Matrix3x3 const& mat, riemann::Vec3 const& vec3);


	riemann::Vec3& operator+=(riemann::Vec3& first, riemann::Vec3 const& right);


	const Vec3 operator+(Vec3 const& left, Vec3 const& right);

	const Vec3 operator-(Vec3 const& left, Vec3 const& right);


	const Vec3 operator*(double alpha, Vec3 const& vec);

	const Vec3 operator*(Vec3 const& vec, double alpha);

	double dot(Vec3 const &first, Vec3 const &second);

	double sgn(double x);

	double minmod(double x, double y);

	Vec3 minmod(riemann::Vec3 const& x, riemann::Vec3 const& y);

	double sqr(double x);

	Vec3 sqr(riemann::Vec3 const& vec);



}









#endif //REOSOLVER_ALGEBRA_H
