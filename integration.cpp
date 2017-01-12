/*
 * integration.cpp
 *
 *  Created on: 2017/01/10
 *      Author: kklab
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include "myMath.h"
using namespace std;

Functions F;

int main() {
    cout << "Hello C++ World" << endl;
    cout << F.trapezoidal_rule(F.gamma_inner, 0, 1000, 1, 10.0) << endl;
    cout << F.simpson_rule(F.gamma_inner, 0, 1000, 1, 10.0) << endl;
    for (int x=0; x<=1000; x++) {
    	cout << F.besseli(1.0*x, 10)
    			<< ", "
    			<< 1.0/sqrt(2.0*M_PI*x) * exp(1.0*x)
				<< endl;
    }
    return 0;
}


