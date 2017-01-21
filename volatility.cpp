/*
 * volatility.cpp
 *
 *  Created on: 2017/01/19
 *      Author: kklab
 */

#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "myMath.h"

using namespace std;

Functions F;

double calc_prob(int u, int d, double p_UU, double p_UD, double p_DU, double p_DD) {
	/*
	 *  変動の系列が与えられた下で変動の確率を計算する．
	 */
	double p_U, p_D, ans;
	int k;
	p_U = p_DU / (p_UD + p_DU); // 初期分布
	p_D = p_UD / (p_UD + p_DU); // 初期分布
	if (d==0) {
		return p_U * pow(p_UU, (u-1));
	}
	if (u==0) {
		return p_D * pow(p_DD, d-1);
	}
	if (u > d) {
		ans = p_D * F.binom(u-1, d-1) * pow(p_UD, d-1) * pow(p_DU, d) * pow(p_UU, u-d);
		for (k=0; k<=d-1; k++) {
			ans += p_U * F.binom(d-1, k) * pow(p_UD, k+1) * pow(p_DU, k) * pow(p_UU, u-2-k) * pow(p_DD, d-1-k) * (F.binom(u-1, k)*p_UU + F.binom(u-1, k+1)*p_DU);
		}
		for (k=0; k<=d-2; k++) {
			ans += p_D * F.binom(u-1, k) * pow(p_UD, k) * pow(p_DU, k+1) * pow(p_UU, u-1-k) * pow(p_DD, d-2-k) * (F.binom(d-1, k)*p_DD + F.binom(d-1, k+1)*p_UD);
		}
	} else if (u < d) {
		ans = p_U * F.binom(d-1, u-1) * pow(p_UD, u) * pow(p_DU, u-1) * pow(p_DD, d-u);
		for (k=0; k<=u-2; k++) {
			ans += p_U * F.binom(d-1, k) * pow(p_UD, k+1) * pow(p_DU, k) * pow(p_UU, u-2-k) * pow(p_DD, d-1-k) * (F.binom(u-1, k)*p_UU + F.binom(u-1, k+1)*p_DU);
		}
		for (k=0; k<=u-1; k++) {
			ans += p_D * F.binom(u-1, k) * pow(p_UD, k) * pow(p_DU, k+1) * pow(p_UU, u-1-k) * pow(p_DD, d-2-k) * (F.binom(d-1, k)*p_DD + F.binom(d-1, k+1)*p_UD);
		}
	} else {
		ans = p_U * pow(p_UD, u) * pow(p_DU, u-1) + p_D * pow(p_UD, u-1) * pow(p_DU, u);
		for (k=0; k<=u-2; k++) {
			ans += p_U * F.binom(d-1, k) * pow(p_UD, k+1) * pow(p_DU, k) * pow(p_UU, u-2-k) * pow(p_DD, d-1-k) * (F.binom(u-1, k)*p_UU + F.binom(u-1, k+1)*p_DU);
			ans += p_D * F.binom(u-1, k) * pow(p_UD, k) * pow(p_DU, k+1) * pow(p_UU, u-1-k) * pow(p_DD, d-2-k) * (F.binom(d-1, k)*p_DD + F.binom(d-1, k+1)*p_UD);
		}
	}
	return ans;
}

int main() {

	string rootpath = "C:\\Users\\kklab\\Desktop\\yurispace\\board_fluctuation\\src\\nikkei_needs_output";
	string subdir = "\\statistics_of_the_limit_order_book\\move_frequency";
	string datayear = "\\2007";
	string sessions[2] = { "\\morning", "\\afternoon" };
	string dirpath = rootpath + subdir + datayear;
	string filenames[4] = { "1min_int.csv", "3min_int.csv", "5min_int.csv", "10min_int.csv" };

	string filepath, session, name, date, rowdate;
	vector<string> data, probmat;
	double p_UU, p_UD, p_DU, p_DD, prob, EX, VX;
	int T, up, down, N;

	for (int s=0; s<sizeof(sessions)/sizeof(sessions[0]); s++) {
		session = sessions[s];
		for (int n=0; n<sizeof(filenames)/sizeof(filenames[0]); n++) {
			name = filenames[n];
			filepath = dirpath + session + name;
			data = fio.readcsv(filepath, false);
			vector<double> table(2, 0.0); //期待値と分散を計算する．

			for (vector<string>::iterator itr=data.begin(); itr!=data.end(); ++itr) {
				vector<double> vec(2, 0.0);
				date = fio.split(*itr, ',').at(0);
				cout << date << session << name <<endl;
				filepath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source\\2007\\probability_10pieces.csv"; //確率を取得する．
				probmat = fio.readcsv(filepath, true);
				rowdate = date.substr(1,4) + "/" + date.substr(5,2) + "/" +date.substr(7,2) + "/" + session.substr(1);

				for (vector<string>::iterator itr2=probmat.begin(); itr2!=probmat.end(); ++itr2) {
					if (fio.split(*itr2, ',').at(0) == rowdate) {
						p_UU = im.stod(fio.split(*itr2, ',').at(1));
						p_UD = im.stod(fio.split(*itr2, ',').at(2));
						p_DU = im.stod(fio.split(*itr2, ',').at(3));
						p_DD = im.stod(fio.split(*itr2, ',').at(4));
						break;
					}
				}
				T=0;
				for (int p=1; p<fio.split(*itr, ',').size(); p++) {
					T += im.stoi(fio.split(*itr, ',').at(p)); // 変動回数
				}
				EX=0.0;
				VX=0.0;
				double pp = 0.0;
				for (N=-T; N<=T; N+=2) {
					up = (T+N)/2;
					down = (T-N)/2;
					prob = calc_prob(up, down, p_UU, p_UD, p_DU, p_DD);
					pp += prob;
					EX += N * prob;
					VX += N * N * prob;
				}
				cout << pp << endl;
				VX = VX - EX * EX;
				//cout << EX << ", " << VX << endl;
			}
		}
	}

	return 0;
}
