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
const double x_max = 1000.0;

/* density function of operating time */
inline double operating_time(double *x) {
	/*
	 * (r/x) * exp(-(lambda + mu)*x) * pow(mu/lambda, r/2) * besselI(2*x*sqrt(lambda*mu), r)
	 * x[0] :x
	 * x[1] :r
	 * x[2] :lambda
	 * x[3] :mu
	 */
	return (x[1]/x[0]) * exp(-(x[2] + x[3])*x[0]) * pow(x[3]/x[2], x[1]/2) * F.besselI(2*x[0]*sqrt(x[2]*x[3]), x[1]);
}

inline double f_u(double *x) {
	/*
	 * f_A(t, r_A, l_A, m_A) * (1 - simpson_rule(f_B, 1e-10, t, r_B, l_B, m_B)$value)
	 * x[0] :x
	 * x[1] :r_A
	 * x[2] :lambda_A
	 * x[3] :mu_A
	 * x[4] :r_B
	 * x[5] :lambda_B
	 * x[6] :mu_B
	 * */
	double args_A[4];
	args_A[0] = x[0];
	for (int i=1; i<4; i++) {
		args_A[i] = x[i];
	}
	return operating_time(args_A) * (1 - F.simpson_rule(operating_time, x0, x[0], 3, x[4], x[5], x[6]));
}
double ot_trans(double x, int r, double lambda, double mu) {
	double rhoinv = mu/lambda;
	double besarg = 2*log(x)*sqrt(lambda*mu)/(lambda+mu);
	return(pow(-1, r+1) * (r/log(x)) * pow(rhoinv, r/2) * F.besselI(besarg, r));
}

double expinv(double *x) {
	return exp(-x[0]);
}
int main() {
    cout << "Hello C++ World" << endl;
    cout << F.simpson_rule(F.gamma_inner, 0, 1000, 1, 10.0) << endl;
    cout << F.simpson_rule(expinv, 0, 100, 0) << endl;
    //cout << F.simpson_rule(operating_time, x0, x_max, 3, 10.0, 0.3, 0.4) << endl;
    //cout << F.simpson_rule(f_u, x0, x_max, 6, 10.0, 0.3, 0.4, 10.0, 0.3, 0.4) << endl;
    double arg[7] = {1000, 586, 2.682564, 3.104103, 58, 2.699103, 2.887436};
    cout << f_u(arg) << endl;
    /*
    double res;
    double tmp;
    double diff = 10;
    int N = 128;
    tmp = simpson(N, x0, x_max, 1e-08, 5, 0.3, 0.4);
    while (diff >= 1e-08) {
    	N *= 2;
    	res = simpson(N, x0, x_max, 1e-08, 5, 0.3, 0.4);
    	diff = res - tmp;
    	tmp = res;
    }
    cout << res << endl;
    */
    return 0;
}

