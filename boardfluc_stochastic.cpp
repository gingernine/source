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
const double x_max = 10000.0;
const int n = 8;

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
	return operating_time(args_A) * (1 - F.gauss_legendre(operating_time, x0, x[0], n, 3, x[4], x[5], x[6]));
	//return operating_time(args_A) * F.romberg(operating_time, x0, x[0], 3, x[4], x[5], x[6]);
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
	return operating_time(args_B) * (1 - F.gauss_legendre(operating_time, x0, x[0], n, 3, x[1], x[2], x[3]));
	//return operating_time(args_B) * F.romberg(operating_time, x0, x[0], 3, x[1], x[2], x[3]);
}

double gamma_dist(double x, double alpha, double beta) {
	double gamma_alpha;
	if (alpha == 3.805877022) {
		gamma_alpha = 4.727357;
	} else if (alpha == 0.673427365) {
		gamma_alpha = 1.342197;
	} else if (alpha == 0.680112458) {
		gamma_alpha = 1.330693;
	} else if (alpha == 3.849663935) {
		gamma_alpha = 4.983704;
	}
	return (1/(gamma_alpha*pow(beta, alpha))) * pow(x,alpha-1.0) * exp(-x/beta);
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
    string dirpath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\gauss_roots";
    string filename = "\\legrootsweights.csv";
    vector<string> strleg1 = fio.readcsv(dirpath+filename, true);
    vector<string> strleg2 = fio.readcsv(dirpath+filename, true);
    vector<string> strleg3 = fio.readcsv(dirpath+filename, true);

    string rootpath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source";
    string datayear = "\\2007";
    string filepath = "\\parameters_1pieces.csv";
    double r_U_A, r_D_A, r_U_B, r_D_B, lambda_A, lambda_B, mu_A, mu_B;
    vector<string> strdata = fio.readcsv(rootpath+datayear+filepath, true);
    for (vector<string>::iterator itr = strdata.begin(); itr!=strdata.end(); ++itr) {
    	lambda_B = im.stod(fio.split(*itr, ',').at(9));
    	lambda_A = im.stod(fio.split(*itr, ',').at(10));
    	mu_A = im.stod(fio.split(*itr, ',').at(11));
    	mu_B = im.stod(fio.split(*itr, ',').at(12));

    	double r1,r2,r3, w1,w2,w3, sum_UU = 0.0, sum_UD = 0.0, sum_DU = 0.0, sum_DD = 0.0;
    	for (vector<string>::iterator itr1 = strleg1.begin(); itr1!=strleg1.end(); ++itr1) {
    		if(n == im.stoi(fio.split(*itr1, ',').at(0))) {
    			r1 = im.stod(fio.split(*itr1, ',').at(1));
    			w1 = im.stod(fio.split(*itr1, ',').at(2));
    			for (vector<string>::iterator itr2 = strleg2.begin(); itr2!=strleg2.end(); ++itr2) {
    				if(n == im.stoi(fio.split(*itr2, ',').at(0))) {
    					r2 = im.stod(fio.split(*itr2, ',').at(1));
    					w2 = im.stod(fio.split(*itr2, ',').at(2));
    					for (vector<string>::iterator itr3 = strleg3.begin(); itr3!=strleg3.end(); ++itr3) {
    						if(n == im.stoi(fio.split(*itr3, ',').at(0))) {
    							r3 = im.stod(fio.split(*itr3, ',').at(1));
    							w3 = im.stod(fio.split(*itr3, ',').at(2));
    							r1 = (x_max - x0)*0.5*r1 + (x_max + x0)*0.5;
    							r2 = (x_max - x0)*0.5*r2 + (x_max + x0)*0.5;
    							r3 = (x_max - x0)*0.5*r3 + (x_max + x0)*0.5;
    							double x[7] = {r1, r2, lambda_A, mu_A, r3, lambda_B, mu_B};
    							sum_UU += w1*w2*w3*f_u(x) * gamma_dist(r2, 3.805877022, 110.4830751) * gamma_dist(r3, 0.680112458, 59.24521486);
    							sum_UD += w1*w2*w3*f_d(x) * gamma_dist(r2, 3.805877022, 110.4830751) * gamma_dist(r3, 0.680112458, 59.24521486);
    							sum_DU += w1*w2*w3*f_u(x) * gamma_dist(r2, 0.673427365, 59.86672222) * gamma_dist(r3, 3.849663935, 109.5604804);
    							sum_DD += w1*w2*w3*f_d(x) * gamma_dist(r2, 0.673427365, 59.86672222) * gamma_dist(r3, 3.849663935, 109.5604804);
    						}
    					}
    				}
    		    }
    		}
    	}

    	cout << (x_max - x0)*0.5*sum_UU << endl;
    	cout << (x_max - x0)*0.5*sum_UD << endl;
    	cout << (x_max - x0)*0.5*sum_DU << endl;
    	cout << (x_max - x0)*0.5*sum_DD << endl;
    }
    return 0;
}

