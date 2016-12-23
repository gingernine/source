/*
 * simpson_rule.cpp
 *
 *  Created on: 2016/12/23
 *      Author: kklab
 */

#include <iostream>
#include <cmath>
using namespace std;

/* gamma function */
const double alpha = 9.0;
double f(double x, double alpha){
    return (exp(-x)) * pow(x, alpha-1);
}

/* integration interval */
const double x0 = 0.0000001; // x0を０としないのは、xが０の時に発散するため
const double x_max = 100000.0;

/* simpson rule */
double simpson(double lower,double upper, double error){
	int n = (int)pow( 1/error * pow(upper-lower, 5) / 2880, 0.25);
    double delta = (upper - lower) / (2*n);
    double sum;
    for(int i=1; i<=n; i++){
        sum += f(lower + delta * (2*i-2), alpha)
        		+ 4 * f(lower + delta * (2*i-1), alpha)
				+ f(lower + delta * 2*i, alpha);
    }
    return sum*delta/3.0;
}

int main() {
    cout << "Hello C++ World" << endl;
    cout << simpson(x0, x_max, 1e-08) << endl;
    return 0;
}
