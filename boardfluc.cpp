/*
 * board_fluc.cpp
 *
 *  Created on: 2016/12/23
 *      Author: kklab
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include "myMath.h"
using namespace std;

Functions F;

/* integration interval */
const double x0 = 1e-10;
const double x_max = 500.0;

/* the density function of operating time */
inline double operating_time(double *x) {
	/*
	 * (r/x) * exp(-(lambda + mu)*x) * pow(mu/lambda, r/2) * besseli(2*x*sqrt(lambda*mu), r)
	 * x[0] :x, x[1] :r, x[2] :lambda, x[3] :mu
	 */
	return (x[1]/x[0]) * exp(-(x[2] + x[3])*x[0]) * pow(x[3]/x[2], x[1]/2) * F.besseli(2*x[0]*sqrt(x[2]*x[3]), (int)x[1]);
}

inline double f_u(double *x) {
	/*
	 * f_A(t, r_A, l_A, m_A) * (1 - integrate(f_B, 1e-10, t, r_B, l_B, m_B))
	 * x[0] :x
	 * x[1] :r_A, x[2] :lambda_A, x[3] :mu_A
	 * x[4] :r_B, x[5] :lambda_B, x[6] :mu_B
	 */
	double args_A[4];
	for (int i=0; i<4; i++) {
		args_A[i] = x[i];
	}
	return operating_time(args_A) * F.gauss_legendre(operating_time, x0, x[0], 32, 3, x[4], x[5], x[6]);
}

inline double f_d(double *x) {
	/*
	 * f_B(t, r_B, l_B, m_B) * (1 - integrate(f_A, 1e-10, t, r_A, l_A, m_A))
	 * x[0] :x
	 * x[1] :r_A, x[2] :lambda_A, x[3] :mu_A
	 * x[4] :r_B, x[5] :lambda_B, x[6] :mu_B
	 */
	double args_B[4];
	args_B[0] = x[0];
	for (int i=1; i<4; i++) {
		args_B[i] = x[i+3];
	}
	return operating_time(args_B) * F.gauss_legendre(operating_time, x0, x[0], 32, 3, x[4], x[5], x[6]);
}

int main() {
    cout << "Hello C++ World" << endl;
    cout << F.simpson_rule(F.gamma_inner, 0, 1000, 1, 10.0) << endl;
    cout << F.gauss_legendre(F.gamma_inner, x0, 1000, 80, 1, 10.0) << endl;
    cout << F.gauss_legendre(f_u, x0, x_max, 32, 6, 586.0, 2.682564, 3.104103, 58.0, 2.699103, 2.887436) << endl;
    return 0;
}

