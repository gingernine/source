/*
 * gauss_integration_weights.cpp
 *
 *  Created on: 2017/01/12
 *      Author: kklab
 */

#include <cmath>
#include <iomanip>
#include "fileio.h"
#include "myMath.h"
using namespace std;

Fileio fio;
InfoManager im;

int main() {
	string dirpath = "C:\\Users\\kklab\\Desktop\\yurispace\\integration_cpp\\gauss_roots";
	string lagroots = "\\lagroots.csv";
	string legroots = "\\legroots.csv";
	string lagweights = "\\lagweights.csv";
	string legweights = "\\legweights.csv";
	vector<string> data = fio.readcsv(dirpath + lagweights, true);
	for (vector<string>::iterator itr = data.begin(); itr!=data.end(); ++itr) {
		if( 10 == im.stoi(fio.split(*itr, ',').at(0))) {
			cout << fixed << setprecision(40) << fio.split(*itr, ',').at(1) << endl;
		}
	}

	return 0;
}

