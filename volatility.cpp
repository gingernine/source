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
#include "fileio.h"

using namespace std;

InfoManager im;
Fileio fio;

vector<int> dec2bin(int n, int m) {
	/* 10進数 n を2進数ベクトルに変換する．返すベクトルの長さを m に揃える． */
	vector<int> bin(m, 0);
	if (n == 0) {
		return bin;
	}
	int a, i=0;
	while (n != 0) {
		a = n % 2;
		bin.at(i) = a;
		n = (n-a)/2;
		i++;
	}
	return bin;
}

double calc_prob(vector<int> pattern, double p_UU, double p_UD, double p_DU, double p_DD) {
	/*
	 *  変動の系列が与えられた下で変動の確率を計算する．
	 *  pattern の要素は0:down, 1:up と対応させる．
	 */
	double ans;
	if (pattern.at(0) == 1) {
		ans = p_DU/(p_DU + p_UD);
	} else {
		ans = p_UD/(p_DU + p_UD);
	}
	for (int i=0; i<=pattern.size()-2; ++i) {
		if (pattern.at(i)==1 && pattern.at(i+1)==1) {
			ans *= p_UU;
		} else if (pattern.at(i)==1 && pattern.at(i+1)==0) {
			ans *= p_UD;
		} else if (pattern.at(i)==0 && pattern.at(i+1)==1) {
			ans *= p_DU;
		} else if (pattern.at(i)==0 && pattern.at(i+1)==0) {
			ans *= p_DD;
		}
	}
	return ans;
}

int main() {

	int m=0;
	for (int i=1; i<1; i++){
		m++;
	}
	cout << m << endl;
	/*
	string rootpath = "C:\\Users\\kklab\\Desktop\\yurispace\\board_fluctuation\\src\\nikkei_needs_output";
	string subdir = "\\statistics_of_the_limit_order_book\\move_frequency";
	string datayear = "\\2007";
	string sessions[2] = { "\\morning", "\\afternoon" };
	string dirpath = rootpath + subdir + datayear;
	string filenames[4] = { "1min_int.csv", "3min_int.csv", "5min_int.csv", "10min_int.csv" };

	string filepath, session, name, date, rowdate;
	vector<string> data, probmat;
	vector<int> pattern;
	double p_UU, p_UD, p_DU, p_DD, prob, X, EX, VX;
	int m, u;

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
				filepath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source\\2007\\probability_1pieces.csv"; //確率を取得する．
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
				m=0;
				for (int p=1; p<fio.split(*itr, ',').size(); p++) {
					m += im.stoi(fio.split(*itr, ',').at(p)); // 変動回数
				}
				EX=0.0;
				VX=0.0;
				pattern.clear();
				for(int i=0; i<pow(2,20); i++) {
					pattern = dec2bin(i, m); // paths
					u = std::accumulate(pattern.begin(), pattern.end(), 0);
					X = 10.0 * u - 10.0 * (m-u);
					prob = calc_prob(pattern, p_UU, p_UD, p_DU, p_DD);
					EX += X * prob;
					VX += X * X * prob;
				}
				VX = VX - EX * EX;
				cout << EX << ", " << VX << endl;
			}
		}
	}
	*/
	return 0;
}
