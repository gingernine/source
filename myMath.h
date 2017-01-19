/*
 * myMath.h
 *
 *  Created on: 2016/12/31
 *      Author: kklab
 */

#ifndef MYMATH_H_
#define MYMATH_H_

#include <cmath>
#include <limits>
#include <stdarg.h>
#include <iostream>
#include <vector>
#include "fileio.h"

using namespace std;
InfoManager im;
Fileio fio;

class Functions {
public:
	Functions(void){} // constructor
	virtual ~Functions(void){} // destructor

	double besseli0(double x, double e) {
		/* Return the modified Bessel functioin of the first order
		 * for any real x and n = 0 */

		double ax, ans, y;
		if((ax=fabs(x)) < 3.75) {
			y = x / 3.75;
			y *= y;
			ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
						+ y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
			ans *= exp(-e);
		} else {
			y = 3.75 / ax;
			ans = (exp(ax-e) / sqrt(ax)) * (0.39894228 + y*(0.1328592e-1
						+ y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2
						+ y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1
						+ y*0.392377e-2))))))));
		}
		return ans;
	}

	double besseli1(double x, double e) {
		/* Return the modified Bessel functioin of the first order
			 * for any real x and n = 1 */

		double ax, ans, y, t;
		if((ax=fabs(x)) < 3.75) {
			y = x / 3.75;
			y *= y;
			ans = ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
					+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
			ans *= exp(-e);
		} else {
			y = 3.75 / ax;
			ans = 0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
					-y*0.420059e-2));
			ans = 0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
					+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
			ans *= (exp(ax-e)/sqrt(ax));
		}
		return x < 0.0 ? -ans : ans;
	}

	double besseli(double x, int n, double e) {
		/* Return the modified Bessel functioin of the first order In(x)
		 * for any real x and any integer n */

		const double ACC = 200.0;
		const int IEXP = numeric_limits<double>::max_exponent/2;
		int j, k;
		double bi, bim, bip, dum, tox, ans;

		if (n == 0) {
			return besseli0(x, e);
		}
		if (n == 1) {
			return besseli1(x, e);
		}
		if (x*x <= 8.0*numeric_limits<double>::min()) {
			return 0.0;
		} else {
			tox = 2.0 / fabs(x);
			bip = ans = 0.0;
			bi = 1.0;
			for (j=2*(n+int(sqrt(ACC*n))); j>0; j--) {
				bim = bip + j * tox * bi;
				bip = bi;
				bi = bim;
				dum = frexp(bi, &k);
				if (k > IEXP) {
					ans = ldexp(ans, -IEXP);
					bi = ldexp(bi, -IEXP);
					bip = ldexp(bip, -IEXP);
				}
				if (j == n) {
					ans = bip;
				}
			}
			ans *= besseli0(x, e) / bi;
			return x < 0.0 ? -ans : ans;
		}
	}

	static double gamma_inner(double *x){
		/* gamma function */
	    return exp(-x[0]) * pow(x[0], x[1]-1);
	}

	inline static double term_by_term() {
		double sum;
		return sum;
	}

	inline static double gauss_legendre(double (*func)(double *x), double lower, double upper, int n, int len, ...) {
		/*
		 * trapezoidal rule
		 * 可変長引数 ... は，関数が複数パラメータをもつ場合の対処
		 * n     : legendre多項式の次数
		 * len   : 被積分関数のパラメータの個数．なければ 0 を渡す．
		 * args  : 関数に渡す引数の配列
		 * sum   : 積分近似値
		 */
		double args[len+1];
		va_list ARGS;
		va_start(ARGS, len);
		if (len >= 1) {
			for (int i=1; i<=len; i++) {
				args[i] = va_arg(ARGS, double);
			}
		}
		va_end(ARGS);
		string dirpath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\gauss_roots";
		string filename = "\\legrootsweights.csv";
		vector<string> strdata = fio.readcsv(dirpath+filename, true);
		double root, weight, sum = 0.0;
		for (vector<string>::iterator itr = strdata.begin(); itr!=strdata.end(); ++itr) {
			if(n == im.stoi(fio.split(*itr, ',').at(0))) {
				root = im.stod(fio.split(*itr, ',').at(1));
				weight = im.stod(fio.split(*itr, ',').at(2));
				args[0] = (upper - lower)*0.5*root + (upper + lower)*0.5;
				sum += weight * func(args);
			}
			if (n < im.stoi(fio.split(*itr, ',').at(0))) {
				break;
			}
		}
		return sum * (upper - lower) * 0.5;
	}

	inline static double gauss_chebyshev(double (*func)(double *x), double lower, double upper, int n, int len, ...) {
			/*
			 * trapezoidal rule
			 * 可変長引数 ... は，関数が複数パラメータをもつ場合の対処
			 * n     : chebyshev多項式の次数
			 * len   : 被積分関数のパラメータの個数．なければ 0 を渡す．
			 * args  : 関数に渡す引数の配列
			 * sum   : 積分近似値
			 */
			double args[len+1];
			va_list ARGS;
			va_start(ARGS, len);
			if (len >= 1) {
				for (int i=1; i<=len; i++) {
					args[i] = va_arg(ARGS, double);
				}
			}
			va_end(ARGS);
			double root, weight = M_PI/n, sum = 0.0;
			for (int i=1; i<=n; i++) {
				root = cos((2.0*i-1.0)*M_PI/(2.0*n));
				args[0] = (upper - lower)*0.5*root + (upper + lower)*0.5;
				sum += sqrt(1.0-root*root) * func(args);
			}
			return weight * (upper - lower) * 0.5 * sum;
		}

	inline static double trapezoidal_rule(double (*func)(double *x), double lower, double upper, int len, ...) {
		/*
		 * trapezoidal rule
		 * 可変長引数 ... は，関数が複数パラメータをもつ場合の対処
		 * n     : 積分区間分割数
		 * len   : 被積分関数のパラメータの個数．なければ 0 を渡す．
		 * args  : 関数に渡す引数の配列
		 * delta : 積分区間分割幅
		 * sum   : 積分近似値
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
		double delta, sum=0.0, tmp=1.0;
		while (abs(sum - tmp) >= 1e-4) {
			cout << tmp << " with error " << abs(sum-tmp) << endl;
			delta = 1.0 * (upper - lower) / n;
			tmp = sum; // update
			args[0] = lower;
			sum = func(args) / 2.0;
			args[0] = upper;
			sum += func(args) / 2.0;
			for(int i=1; i<n; i++){
				args[0] = lower + delta * i;
				sum += func(args);
			}
			sum *= delta;
			n *= 2;
		}
		return sum;
	    va_end(ARGS);
	}

	inline static double simpson_rule(double (*func)(double *x), double lower, double upper, int len, ...) {
		/*
		 * Simpson rule
		 * 可変長引数 ... は，関数が複数パラメータをもつ場合の対処
		 * n     : 2 * n が積分区間分割数
		 * len   : 被積分関数のパラメータの個数．なければ 0 を渡す．
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
		double delta, sum=0.0, sum_o=0.0, sum_e=0.0, tmp=1.0;
		while (abs(sum - tmp) >= 1e-4) {
			delta = 1.0 * (upper - lower) / (2 * n);
			tmp = sum; // update
			sum_e = 0.0;
			sum_o = 0.0;
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
			sum += 2.0 * sum_e + 4.0 * sum_o;
			sum *= delta / 3.0;
			n *= 2;
		}
		return sum;
	    va_end(ARGS);
	}

	inline static double romberg(double (*func)(double *x), double lower, double upper, int len, ...) {
		/*
		 * Romberg method
		 * 可変長引数 ... は，関数が複数パラメータをもつ場合の対処
		 * m     : 漸化式の m (=1,2,3,...)
		 * 2^m   : 積分区間分割数
		 * len   : 被積分関数のパラメータの個数．なければ 0 を渡す．
		 * args  : 関数に渡す引数の配列
		 * delta : 積分区間分割幅
		 * sum   : 積分近似値
		 * temp1 : 漸化式の S_{m-1,k} の値を格納するベクトル．
		 * temp2 : 漸化式の S_{m,k} の値を格納するベクトル．
		 */

		double args[len+1];
		va_list ARGS;
		va_start(ARGS, len);
		if (len >= 1) {
			for (int i=1; i<=len; i++) {
				args[i] = va_arg(ARGS, double);
			}
		}

		int m=1;
		double delta, sum=0.0, diff=0.0;
		vector<double> temp1, temp2;

		/* 台形則で m=0 の分計算する． */
		delta = 1.0 * (upper - lower) / pow(2,m-1);
		args[0] = lower;
		sum = func(args) / 2.0;
		args[0] = upper;
		sum += func(args) / 2.0;
		for(int i=1; i<pow(2,m-1); i++){
			args[0] = lower + delta * i;
			sum += func(args);
		}
		sum *= delta;
		temp1.push_back(sum);

		while (true) {
			/* 台形則で積分値を計算する． */
			delta = 1.0 * (upper - lower) / pow(2,m);
			args[0] = lower;
			sum = func(args) / 2.0;
			args[0] = upper;
			sum += func(args) / 2.0;
			for(int i=1; i<pow(2,m); i++){
				args[0] = lower + delta * i;
				sum += func(args);
			}
			sum *= delta;

			/* 漸化式を利用する． */
			temp2.push_back(sum);
			diff = abs(sum - temp1[m-1]);
			for (int k=1; k<=m; k++) {
				if (diff < 1e-10 && m > 5) {
					return temp2[k-1];
				}
				temp2.push_back((pow(4,k)*temp2[k-1] - temp1[k-1]) / (pow(4,k) - 1));
				diff = abs(temp2[k] - temp2[k-1]);
			}
			if (diff < 1e-10 && m > 5) {
				return temp2[m];
			}
			temp1 = temp2;
			temp2.clear();
			m += 1;
		}
		va_end(ARGS);
	}
};


#endif /* MYMATH_H_ */
