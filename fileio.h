/*
 * fileio.h
 *
 *  Created on: 2016/12/29
 *      Author: kklab
 */

#ifndef FILEIO_H_
#define FILEIO_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Windows.h>

using namespace std;

class InfoManager {
public:
	static int stoi(string str) {
		int ret;
		stringstream ss;
		ss << str;
		ss >> ret;
		return ret;
	}
};

class Fileio {
public:
	Fileio(void){} // constructor
	virtual ~Fileio(void){} // destructor

	static vector<string> listfiles(string dirpath, string exp) {
		// declaration
		vector<string> filelist;
		HANDLE hFind;
		WIN32_FIND_DATA fd;

		stringstream ss;
		ss << dirpath;
		string::iterator itr = dirpath.end();
		itr--;
		if(*itr != '\\') {
			ss << '\\';
		}
		if (exp != "") {
			ss << "*." << exp; // expansion, e.g.".csv", ".txt"
		} else {
			ss << "*.*";
		}

		// search files
		// FindFirstFile(filename, &fd)
		hFind = FindFirstFile(ss.str().c_str(), &fd);

		if(hFind == INVALID_HANDLE_VALUE) {
			cout << "Couldn't get list of files" << endl;
			exit(1);
		}

		do {
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
				&& !(fd.dwFileAttributes & FILE_ATTRIBUTE_HIDDEN) ) {
				char *file = fd.cFileName;
				string str = file;
				filelist.push_back(str);
			}
		} while (FindNextFile(hFind, &fd));

		FindClose(hFind);
		return filelist;
	}

	static vector<string> split(string& input, char delimiter) {
	    istringstream stream(input);
	    string field;
	    vector<string> result;
	    while (getline(stream, field, delimiter)) {
	        result.push_back(field);
	    }
	    return result;
	}

	static vector<string> readcsv(string filepath, bool header) {
		string readline;
		ifstream ifs(filepath.c_str());
		vector<string> strvec;
		int p = 0;
		while(getline(ifs, readline)){
			if (header && p==0) {
				p++;
				continue;
			}
			strvec.push_back(readline);
		}
		return strvec;
	}
};


#endif /* FILEIO_H_ */
