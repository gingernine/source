/*
 * arrival_rate.cpp
 *
 *  Created on: 2016/12/29
 *      Author: kklab
 */

#include <iostream>
#include <vector>
#include "fileio.h"

using namespace std;

Fileio fio;
InfoManager im;

int convert_to_seconds(int time) {
	/* "hhmmss"(時分秒)表示される時間を始点0時0分0秒からの秒で表示する． */
	int h = (time - time % 10000) * 36 / 100;
	int m = (time % 10000 - time % 100) * 6 / 10;
	int s = time % 100;
	return h + m + s;
}

double* calc_arrival_rate(vector<string> time_vector, int interval) {
	/* interval 分間隔での到着率を計算する． */
	int size = time_vector.size();
	int stamp = convert_to_seconds(im.stoi(fio.split(time_vector[0], ',').at(0))); // time stamp for while loop
	int finish = convert_to_seconds(im.stoi(fio.split(time_vector[size-1], ',').at(0))); // closing time
	interval = interval * 60; // time interval in seconds
	int p1 = 1;
	int p2 = 0;
	int counter;
	double *retvec = new double[size];
	for (vector<string>::iterator itr=++time_vector.begin(); itr!=time_vector.end(); ++itr) {
		counter = 0;
		if (fio.split(*itr, ',').at(4)) continue;
		while (convert_to_seconds(im.stoi(fio.split(*itr, ',').at(0))) <= stamp + interval) {
			counter++;
			p1++;
			if (p1==size) break;
		}
		retvec[p2] = counter / interval;
		p2++;
		stamp += interval;
	}
	return retvec;
}

int main() {

	string maindir = "C:\\Users\\kklab\\Desktop\\yurispace\\board_fluctuation\\src\\nikkei_needs_output";
	string subdir = "\\statistics_of_the_limit_order_book\\arrival_time_series";
	string datayear = "\\2007";
	string branchs[4] = { "\\limit_buy", "\\limit_sell", "\\market_buy", "\\market_sell" };
	string sessions[2] = { "\\morning", "\\afternoon" };
	vector<string> filelist1, filelist2;
	string dirpath1, dirpath2, filepath, readline;
	stringstream ssfilepath;
	int size, p1, p2;
	bool boolean;

	for (int b=0; b<4; b++) {
		dirpath1 = maindir + subdir + datayear + branchs[b] + sessions[0];
		dirpath2 = maindir + subdir + datayear + branchs[b] + sessions[1];
		filelist1 = fio.listfiles(dirpath1, "csv");
		filelist2 = fio.listfiles(dirpath2, "csv");
		size = filelist1.size() + filelist2.size();
		//double rate_vec[size] = {};
		p1 = 0;
		p2 = 0;
		boolean = true;
		for(int i = 0; i < size; i++){
			filepath = "";
			if (i == 0) {
				filepath = dirpath1 + "\\" + filelist1[p1];
				p1++;
			} else if (boolean) {
				filepath = dirpath1 + "\\" + filelist1[p1];
				p1++;
				boolean = false;
			} else {
				filepath = dirpath2 + "\\" + filelist2[p2];
				p2++;
				boolean = true;
			}

			vector<string> data = fio.readcsv(filepath, false);
			for (vector<string>::iterator itr = data.begin(); itr!=data.end(); ++itr) {
				cout << convert_to_seconds(im.stoi(fio.split(*itr, ',').at(0))) << endl;
				cout << *itr << endl;
			}
		}
	}
	return 0;
}
