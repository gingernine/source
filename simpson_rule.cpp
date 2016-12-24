/*
 * simpson_rule.cpp
 *
 *  Created on: 2016/12/23
 *      Author: kklab
 */

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
using namespace std;

/* bessel function of the first kind */
double besselI(double x, int r){
	double sum = 1;
	double diff = 1;
	int i = 1;
	while (diff >=1e-08) {
		diff = 1;
		for (int j=1; j<=i; j++){
			diff *= (x*x/4) / (j*(r+j));
		}
		sum += diff;
		i++;
	}
	for (int j=0; j<r; j++){
		sum /= r-j;
	}
	return(pow(x/2, r) * sum);
}

/* gamma function */
const double alpha = 9.0;
double f(double x, double alpha){
    return (exp(-x)) * pow(x, alpha-1);
}

/* integration interval */
const double x0 = 1e-10;
const double x_max = 600.0;

/* density function of operating time */
double operating_time(double x, int r, double lambda, double mu){
	return((r/x) * exp(-(lambda + mu)*x) * pow(mu/lambda, r/2) * besselI(2*x*sqrt(lambda*mu), r));
}

double ot_trans(double x, int r, double lambda, double mu){
	double rhoinv = mu/lambda;
	double besarg = 2*log(x)*sqrt(lambda*mu)/(lambda+mu);
	return(pow(-1, r+1) * (r/log(x)) * pow(rhoinv, r/2) * besselI(besarg, r));
}

/* simpson rule */
double simpson(int n, double lower,double upper, double error, int r, double lambda, double mu){
    double delta = (upper - lower) / (2*n);
    double sum_o = 0; // odd
    double sum_e = 0; // even
    double f_lower = operating_time(lower, r, lambda, mu);
    double f_upper = operating_time(upper, r, lambda, mu);
    for(int i=2; i<=2*n-2; i+=2){
        sum_e += operating_time(lower + delta * i, r, lambda, mu);
    }
    for(int i=1; i<=2*n-1; i+=2){
    	sum_o += operating_time(lower + delta * i, r, lambda, mu);
    }
    return (f_lower + 2*sum_e + 4*sum_o + f_upper) * delta / 3.0;
}

int main() {
    cout << "Hello C++ World" << endl;
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
    return 0;
}
