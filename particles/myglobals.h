#pragma once

#include <vector>

#include <Eigen/Core>
typedef Eigen::Vector3d Vector3Type;

//    SPHERE_FIBONACCI_POINTS computes sphere points on a Fibonacci spiral.
//
//    John Burkardt - 21 October 2013 / This code is distributed under the GNU LGPL license.
//
//  Reference:
//
//    Richard Swinbank, James Purser,
//    Fibonacci grids: A novel approach to global modelling. July 2006
//
inline std::vector< Vector3Type > sphere_fibonacci_points ( int n = 100 )
{
	double cphi;
	double i_r8,n_r8;
	double phi,sphi,theta;
	const double pi = 3.141592653589793;

	phi = ( 1.0 + std::sqrt ( 5.0 ) ) / 2.0;
	n_r8 = ( double ) ( n );

	std::vector< Vector3Type > points;

	for (int j = 0; j < n; j++ )
	{
		i_r8 = ( double ) ( - n + 1 + 2 * j );
		theta = 2.0 * pi * i_r8 / phi;
		sphi = i_r8 / n_r8;
		cphi = std::sqrt ( ( n_r8 + i_r8 ) * ( n_r8 - i_r8 ) ) / n_r8;

		points.push_back( Vector3Type(cphi * std::sin ( theta ), cphi * std::cos ( theta ), sphi) );
	}

	return points;
}
