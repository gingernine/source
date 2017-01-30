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
		return p_U * pow(p_UU, u-1);
	}
	if (u==0) {
		return p_D * pow(p_DD, d-1);
	}
	if (u > d) {
		ans = p_U*p_UD*pow(p_UU, u-1)*pow(p_DD, d-1)
				+ p_U*F.binom(u-1, d)*pow(p_UD, d)*pow(p_DU, d)*pow(p_UU, u-1-d)
				+ p_D*p_DU*pow(p_UU, u-1)*pow(p_DD, d-1);
		for (k=1; k<=d-1; k++) {
			ans += pow(p_UD, k)*pow(p_DU, k)*pow(p_UU, u-1-k)*pow(p_DD, d-1-k)*(
					F.binom(u-1, k)*F.binom(d-1, k)*p_U*p_UD + F.binom(u-1, k)*F.binom(d-1, k-1)*p_U*p_DD
					+ F.binom(u-1, k)*F.binom(d-1, k)*p_D*p_DU + F.binom(u-1, k-1)*F.binom(d-1, k)*p_D*p_UU);
		}
	} else if (u < d) {
		ans = p_U*p_UD*pow(p_UU, u-1)*pow(p_DD, d-1)
				+ p_D*F.binom(d-1, u)*pow(p_UD, u)*pow(p_DU, u)*pow(p_DD, d-1-u)
				+ p_D*p_DU*pow(p_UU, u-1)*pow(p_DD, d-1);
		for (k=1; k<=u-1; k++) {
			ans += pow(p_UD, k)*pow(p_DU, k)*pow(p_UU, u-1-k)*pow(p_DD, d-1-k)*(
					F.binom(u-1, k)*F.binom(d-1, k)*p_U*p_UD + F.binom(u-1, k)*F.binom(d-1, k-1)*p_U*p_DD
					+ F.binom(u-1, k)*F.binom(d-1, k)*p_D*p_DU + F.binom(u-1, k-1)*F.binom(d-1, k)*p_D*p_UU);
		}
	} else {
		ans = p_U*p_UD*pow(p_UU, u-1)*pow(p_DD, u-1)
				+ p_D*p_DU*pow(p_UU, u-1)*pow(p_DD, u-1);
		for (k=1; k<=u-1; k++) {
			ans += pow(p_UD, k)*pow(p_DU, k)*pow(p_UU, u-1-k)*pow(p_DD, u-1-k)*(
					F.binom(u-1, k)*F.binom(u-1, k)*p_U*p_UD + F.binom(u-1, k)*F.binom(u-1, k-1)*p_U*p_DD
					+ F.binom(u-1, k)*F.binom(u-1, k)*p_D*p_DU + F.binom(u-1, k-1)*F.binom(u-1, k)*p_D*p_UU);
		}
	}
	return ans;
}

double calc_cor(const double p_UU, const double p_UD, const double p_DU, const double p_DD) {
	/*
	 *  推移確率が与えられた下で変動の系列相関を計算する．
	 */
	double p_U, p_D, ans;
	p_U = p_DU / (p_UD + p_DU); // 初期分布
	p_D = p_UD / (p_UD + p_DU); // 初期分布
	return (p_U*(p_UU-p_UD)-p_D*(p_DU-p_DD)-(p_U-p_D)*(p_U-p_D)) / (p_U+p_D-(p_U-p_D)*(p_U-p_D));
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
	double p_UU, p_UD, p_DU, p_DD, prob, EX, VX, start;
	int T, up, down, N;

	ofstream ofs((dirpath + "\\volatility_curve.csv").c_str());
	double VXarray[10] = {};
	int vxitr;
	ofs << "prob,50,100,150,200,250,300,350,400,450,500,correlation" << endl;
	for(double x=0.001; x < 1.0; x+=0.001) {
		start = 18000.0;
		vxitr = 0;
		for (T=50; T<=500; T+=50) {
			cout << x << ", " << T <<  endl;
			EX = 0.0;
			VX = 0.0;
			for(N=-T; N<=T; N+=2) {
				up = (T+N)/2;
				down = (T-N)/2;
				prob = calc_prob(up, down, 1.0-x, x, x, 1.0-x);
				EX += (log(start + 10.0*N) - log(start)) * prob;
				VX += (log(start + 10.0*N) - log(start)) * (log(start + 10.0*N) - log(start)) * prob;
			}
			VXarray[vxitr] = VX - EX*EX;
			++vxitr;
		}

		ofs << x << ","
			<< VXarray[0] << "," << VXarray[1] << ","
			<< VXarray[2] << "," << VXarray[3] << ","
			<< VXarray[4] << "," << VXarray[5] << ","
			<< VXarray[6] << "," << VXarray[7] << ","
			<< VXarray[8] << "," << VXarray[9] << ","
			<< calc_cor(1.0-x, x, x, 1.0-x)
			<< endl;
	}

	/*
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
				filepath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\source" + datayear + "\\probability3_" + unit + "pieces.csv"; //確率を取得する．
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
			double pp = 0.0;
			start = startval(rootpath + datayear + "\\sessionsep" + session + "\\" + date.substr(1,8) + "_.csv", true);
			if (T == 0) {
				VX = 0.0;
			} else {
				for (N=-T; N<=T; N+=2) {
					up = (T+N)/2;
					down = (T-N)/2;
					prob = calc_prob(up, down, p_UU, p_UD, p_DU, p_DD);
					pp += prob;
					EX += (log(start + 10.0*N) - log(start)) * prob;
					VX += (log(start + 10.0*N) - log(start)) * (log(start + 10.0*N) - log(start)) * prob;
				}
				VX = VX - EX * EX;
			}
			cout << pp << T << endl;
			ofs << rowdate
				<< ","
				<< VX
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
	*/

	return 0;
}
