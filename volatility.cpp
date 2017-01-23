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
#include <float.h>
#include <fstream>
#include "myMath.h"

using namespace std;

Functions F;

int convert_to_seconds(string time) {
	/* "hhmmss"(時分秒)表示される時間を始点0時0分0秒からの秒で表示する． */
	int h = im.stoi(time.substr(0, 2)) * 3600;
	int m = im.stoi(time.substr(2, 2)) * 60;
	int s = im.stoi(time.substr(4, 2));
	return h + m + s;
}

double calc_prob(const int u, const int d, const double p_UU, const double p_UD, const double p_DU, const double p_DD) {
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

double realized_volatility(double timeint, string filepath, bool usequote) {
	vector<string> series = fio.readcsv(filepath, false);
	vector<string>::iterator itr=series.begin();
	double valln, vallntmp, rv=0.0;
	int inttime, timestamp, back, b, k=0;
	string tq, tqjudge="Quote";

	timestamp = convert_to_seconds(fio.split(*itr, ',').at(1).substr(0,2)
			+ fio.split(*itr, ',').at(1).substr(3,2)
			+ fio.split(*itr, ',').at(1).substr(6,2));
	timestamp += timeint;
	if (usequote) {
		while (fio.split(*itr, ',').at(2) != tqjudge) {
			++itr;
			k++;
		}
		vallntmp = log((im.stod(fio.split(*itr, ',').at(5)) + im.stod(fio.split(*itr, ',').at(7))) * 0.5);
	} else {
		tqjudge = "Trade";
		vallntmp = log(im.stod(fio.split(*itr, ',').at(3)));
	}

	for (int i=0; i<series.size()-k; i++) {
		tq = fio.split(*itr, ',').at(2);
		inttime = convert_to_seconds(fio.split(*itr, ',').at(1).substr(0,2)
				+ fio.split(*itr, ',').at(1).substr(3,2)
				+ fio.split(*itr, ',').at(1).substr(6,2));
		if (i == series.size()-k-1) {
			--itr;
			// 場の終わりは時間間隔を狭めてボラティリティを足す．
			while(fio.split(*itr, ',').at(2) != tqjudge) {
				--itr;
			}
			if (usequote) {
				valln = log((im.stod(fio.split(*itr, ',').at(5)) + im.stod(fio.split(*itr, ',').at(7))) * 0.5);
			} else {
				valln = log(im.stod(fio.split(*itr, ',').at(3)));
			}
			rv += (valln - vallntmp) * (valln - vallntmp);
			break;
		}
		if (inttime >= timestamp && tq == tqjudge) {
			--itr;
			back = 1;
			while(fio.split(*itr, ',').at(2) != tqjudge) {
				--itr;
				back++;
			}
			if (usequote) {
				valln = log((im.stod(fio.split(*itr, ',').at(5)) + im.stod(fio.split(*itr, ',').at(7))) * 0.5);
			} else {
				valln = log(im.stod(fio.split(*itr, ',').at(3)));
			}
			for (b=0; b<back; b++) {
				++itr;
			}
			rv += (valln - vallntmp) * (valln - vallntmp);
			timestamp += timeint;
			vallntmp = valln;
		}
		++itr;
	}
	return rv;
}

double startval(string filepath, bool usequote) {
	vector<string> series = fio.readcsv(filepath, false);
	vector<string>::iterator itr=series.begin();

	if (usequote) {
		while (fio.split(*itr, ',').at(2) != "Quote") {
			++itr;
		}
		return (im.stod(fio.split(*itr, ',').at(5)) + im.stod(fio.split(*itr, ',').at(7))) * 0.5;
	} else {
		return im.stod(fio.split(*itr, ',').at(3));
	}
}

int main() {

	string unit = "30";
	bool observed = true;
	string rootpath = "C:\\Users\\kklab\\Desktop\\yurispace\\board_fluctuation\\src\\nikkei_needs_output";
	string subdir = "\\statistics_of_the_limit_order_book\\move_frequency";
	string datayear = "\\2007";
	string sessions[2] = { "\\morning", "\\afternoon" };
	string dirpath = rootpath + subdir + datayear;

	string filepath, session, date, rowdate, wfilepath;
	vector<string> data, probmat, expects;
	double p_UU, p_UD, p_DU, p_DD, prob, EX, VX;
	int T, up, down, N;

	for (int s=0; s<sizeof(sessions)/sizeof(sessions[0]); s++) {
		session = sessions[s];
		filepath = dirpath + session + "1min_int.csv";
		data = fio.readcsv(filepath, false);
		vector<double> table(2, 0.0); //期待値と分散を計算する．

		ofstream ofs((dirpath + session + "_" + unit + "unit_volatility_observed.csv").c_str());

		for (vector<string>::iterator itr=data.begin(); itr!=data.end(); ++itr) {
			vector<double> vec(2, 0.0);
			date = fio.split(*itr, ',').at(0);
			cout << date << session << endl;
			if (observed) {
				filepath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source" + datayear + "\\probability_observed.csv"; //確率を取得する．
			} else {
				filepath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source" + datayear + "\\probability2_" + unit + "pieces.csv"; //確率を取得する．
			}
			probmat = fio.readcsv(filepath, true);
			filepath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source" + datayear + "\\Expectations_" + unit + "pieces.csv"; //変動回数を取得する．
			expects = fio.readcsv(filepath, true);
			rowdate = date.substr(1,4) + "/" + date.substr(5,2) + "/" +date.substr(7,2) + "/" + session.substr(1);
			vector<string>::iterator itr3 = expects.begin();
			for (vector<string>::iterator itr2=probmat.begin(); itr2!=probmat.end(); ++itr2) {
				if (fio.split(*itr2, ',').at(0) == rowdate) {
					p_UU = im.stod(fio.split(*itr2, ',').at(1));
					p_UD = im.stod(fio.split(*itr2, ',').at(2));
					p_DU = im.stod(fio.split(*itr2, ',').at(3));
					p_DD = im.stod(fio.split(*itr2, ',').at(4));
					T = im.stoi(fio.split(*itr3, ',').at(3));
					break;
				}
				++itr3;
			}
			EX=0.0;
			VX=0.0;
			double start = startval(rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", false);
			for (N=-T; N<=T; N+=2) {
				up = (T+N)/2;
				down = (T-N)/2;
				prob = calc_prob(up, down, p_UU, p_UD, p_DU, p_DD);
				EX += (log(start + 10.0*N) - log(start)) * prob;
				VX += (log(start + 10.0*N) - log(start)) * (log(start + 10.0*N) - log(start)) * prob;
			}
			VX = VX - EX * EX;

			ofs << rowdate
				<< ","
				<< VX
				<< ","
				<< realized_volatility(60.0,
					rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", false)
				<< ","
				<< realized_volatility(300.0,
					rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", false)
				<< ","
				<< realized_volatility(600.0,
					rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", false)
				<< ","
				<< realized_volatility(60.0,
					rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", true)
				<< ","
				<< realized_volatility(300.0,
					rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", true)
				<< ","
				<< realized_volatility(600.0,
					rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", true)
				<<
				endl;
			cout << rowdate
							<< ","
							<< VX
							<< ","
							<< realized_volatility(60.0,
								rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", false)
							<< ","
							<< realized_volatility(300.0,
								rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", false)
							<< ","
							<< realized_volatility(600.0,
								rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", false)
							<< ","
							<< realized_volatility(60.0,
								rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", true)
							<< ","
							<< realized_volatility(300.0,
								rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", true)
							<< ","
							<< realized_volatility(600.0,
								rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", true)
							<<
							endl;
		}
	}
	return 0;
}
