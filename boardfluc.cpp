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
const double x_max = 100000.0;
const int n = 500;

/* the density function of operating time */
inline double operating_time(double *x) {
	/*
	 * (r/x) * exp(-(lambda + mu)*x) * pow(mu/lambda, r/2) * besseli(2*x*sqrt(lambda*mu), r)
	 * x[0] :x, x[1] :r, x[2] :lambda, x[3] :mu
	 */
	return (x[1]/x[0]) * pow(x[3]/x[2], x[1]/2) * F.besseli(2*x[0]*sqrt(x[2]*x[3]), (int)x[1], (x[2] + x[3])*x[0]);
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
	//return operating_time(args_A) * (1 - F.gauss_legendre(operating_time, x0, x[0], n, 3, x[4], x[5], x[6]));
	//return operating_time(args_A) * F.romberg(operating_time, x0, x[0], 3, x[4], x[5], x[6]);
	return operating_time(args_A) * (1 - F.gauss_chebyshev(operating_time, x0, x[0], n, 3, x[4], x[5], x[6]));
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
	//return operating_time(args_B) * (1 - F.gauss_legendre(operating_time, x0, x[0], n, 3, x[1], x[2], x[3]));
	return operating_time(args_B) * F.romberg(operating_time, x0, x[0], 3, x[1], x[2], x[3]);
}

double gamma_dist(double x, double alpha, double beta) {
	double a[2] = {x, alpha};
	return (1/(F.gamma_inner(a)*pow(beta, alpha))) * pow(x,alpha-1.0) * exp(-x/beta);
}

double ope_prob(double lambda, double mu, double r) {
	if (lambda > mu) {
		return pow(mu/lambda, r);
	} else {
		return 1.0;
	}
}

int main() {
    cout << "Hello C++ World" << endl;
    string rootpath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source";
    string datayear = "\\2007";
    string filepath = "\\parameters_1pieces.csv";
    double r_U_A, r_D_A, r_U_B, r_D_B, lambda_A, lambda_B, mu_A, mu_B;
    vector<string> strdata = fio.readcsv(rootpath+datayear+filepath, true);
    for (vector<string>::iterator itr = strdata.begin(); itr!=strdata.end(); ++itr) {
    	r_U_B = im.stod(fio.split(*itr, ',').at(5));
    	r_U_A = im.stod(fio.split(*itr, ',').at(6));
    	r_D_B = im.stod(fio.split(*itr, ',').at(7));
    	r_D_A = im.stod(fio.split(*itr, ',').at(8));
    	lambda_B = im.stod(fio.split(*itr, ',').at(9));
    	lambda_A = im.stod(fio.split(*itr, ',').at(10));
    	mu_A = im.stod(fio.split(*itr, ',').at(11));
    	mu_B = im.stod(fio.split(*itr, ',').at(12));

    	//cout << F.gauss_legendre(f_u, x0, x_max, n, 6, r_U_A, lambda_A, mu_A, r_U_B, lambda_B, mu_B) << endl;
    	//cout << F.gauss_legendre(f_d, x0, x_max, n, 6, r_U_A, lambda_A, mu_A, r_U_B, lambda_B, mu_B) << endl;
    	//cout << F.gauss_legendre(f_u, x0, x_max, n, 6, r_D_A, lambda_A, mu_A, r_D_B, lambda_B, mu_B) << endl;
    	//cout << F.gauss_legendre(f_d, x0, x_max, n, 6, r_D_A, lambda_A, mu_A, r_D_B, lambda_B, mu_B) << endl;

    	cout << F.gauss_chebyshev(f_u, x0, x_max, n, 6, r_U_A, lambda_A, mu_A, r_U_B, lambda_B, mu_B) << endl;
    	//cout << "p_UU, " << ope_prob(lambda_A, mu_A, r_U_A) - F.romberg(f_u, x0, x_max, 6, r_U_A, lambda_A, mu_A, r_U_B, lambda_B, mu_B) << endl;
    	//cout << "p_UD, " << ope_prob(lambda_B, mu_B, r_U_B) - F.romberg(f_d, x0, x_max, 6, r_U_A, lambda_A, mu_A, r_U_B, lambda_B, mu_B) << endl;
    	//cout << "p_DU, " << ope_prob(lambda_A, mu_A, r_D_A) - F.romberg(f_u, x0, x_max, 6, r_D_A, lambda_A, mu_A, r_D_B, lambda_B, mu_B) << endl;
    	//cout << "p_DD, " << ope_prob(lambda_B, mu_B, r_D_B) - F.romberg(f_d, x0, x_max, 6, r_D_A, lambda_A, mu_A, r_D_B, lambda_B, mu_B) << endl;
    }
    return 0;
}

