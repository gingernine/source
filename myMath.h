/*
 * myMath.h
 *
 *  Created on: 2016/12/31
 *      Author: kklab
 */

#ifndef MYMATH_H_
#define MYMATH_H_

#include <cmath>
#include <stdarg.h>
#include <iostream>

using namespace std;

class Functions {
public:
	Functions(void){} // constructor
	virtual ~Functions(void){} // destructor

	inline static double besselI(double x, double r){
		/* bessel function of the first kind */
		double sum = 1;
		double diff = 1;
		int i = 1;
		while (diff >=1e-08) {
			diff = 1;
			for (int j=1; j<=i; j++){
				diff *= (x*x*0.25) / (j*(r+j));
			}
			sum += diff;
			i++;
		}
		for (int j=0; j<r; j++){
			sum /= r-j;
		}
		return(pow(x*0.5, r) * sum);
	}

	static double gamma_inner(double *x){
		/* gamma function */
	    return (exp(-x[0])) * pow(x[0], x[1]-1);
	}

	inline static double simpson_rule(double (*func)(double *x), double lower, double upper, int len, ...) {
		/*
		 * simpson rule
		 * 可変長引数 ... は，関数が複数パラメータをもつ場合の対処
		 * 関数が変数の外にパラメータを持たない場合は len = 0 とする．
		 * n     : 2 * n が積分区間分割数
		 * args  : 関数に渡す引数の配列
		 * delta : 積分区間分割幅
		 * sum   : 積分近似値
		 * sum_o : 分点が奇数番目の箇所の積分値総和
		 * sum_e : 分点が偶数番目の箇所の積分値総和
		 * tmp   : 積分値比較用の格納場所
		 */
		int n = 32;
		double args[len+1];
		va_list ARGS;
		va_start(ARGS, len);
		if (len >= 1) {
			for (int i=1; i<=len; i++) {
				args[i] = va_arg(ARGS, double);
			}
		}
		double delta, sum=0, sum_o=0, sum_e=0, tmp=1;
		while (abs(sum - tmp) >= 1e-4) {
			cout << tmp << endl;
			delta = (upper - lower) / (2 * n);
			tmp = sum; // update
			sum_e = 0;
			sum_o = 0;
			args[0] = lower;
			sum = func(args);
			args[0] = upper;
			sum += func(args);
			for(int i=2; i<=2*n-2; i+=2){
				args[0] = lower + delta * i;
				sum_e += func(args);
			}
			for(int i=1; i<=2*n-1; i+=2){
				args[0] = lower + delta * i;
				sum_o += func(args);
			}
			sum += 2 * sum_e + 4 * sum_o;
			sum *= delta / 3.0;
			n *= 2;
		}
		return sum;
	    va_end(ARGS);
	}
};



#endif /* MYMATH_H_ */
