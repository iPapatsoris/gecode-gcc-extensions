#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <gecode/int.hh>
#include <math.h>

using namespace std;
using namespace Gecode;

double euclidianDistance2D(pair<double, double> pointA, 
													 pair<double, double> pointB) {
	double xDiff = pointB.first - pointA.first;
	double yDiff = pointB.second - pointA.second;
	return sqrt(xDiff*xDiff + yDiff*yDiff);
}

void readInput(string fileName, int& vars, IntSetArgs& domain, IntArgs& vals,
							 IntArgs& lowerBounds, IntArgs& upperBounds, IntArgs& costs) {

	string line;
  ifstream file(fileName);
	if (!file.is_open()) {
		throw "Could not open file";
	}

	getline(file, line);
	getline(file, line);
	getline(file, line);
	getline(file, line);
	int n;
	line.erase(std::remove_if(line.begin(), line.end(),
                            [](char c) { return !std::isdigit(c); }),
             line.end());
	stringstream stream(line);
	stream >> vars;	
	for (unsigned int i = 0; i < vars; i++) {
		vals << i;
	}
	for (unsigned int i = 0; i < vars; i++) {
		vector<int> d;
		for (unsigned int j = 0; j < vars; j++) {
			if (i != j) {
				d.push_back(j);
			}
		}
		IntSet tmp;
		IntSetInit<IntArgs>::init(tmp, IntArgs(d));
		domain << IntSet(tmp);
		upperBounds << 1;
		lowerBounds << 1;
	}

	/* 2d euclidian */
	getline(file, line);
	getline(file, line);

	vector<pair<double, double>> cord;
	while (getline(file, line)) {
		stringstream stream(line);
		stream >> n;
		pair<double, double> c(0, 0);
		stream >> c.first;
		stream >> c.second;
		cord.push_back(c);
	}

	for (unsigned int i = 0; i < vars; i++) {
		for (unsigned int j = 0; j < vars; j++) {
			costs << (int) (0.5 + euclidianDistance2D(cord[i], cord[j]));
		}
	}

	/* Full matrix */
	// getline(file, line);
	// getline(file, line);
	// getline(file, line);
	// getline(file, line);

	// for (unsigned int row = 0; row < vars; row++){ 
	// 	getline(file, line);
	// 	stringstream stream(line);
	// 	for (unsigned int i = 0; i < vars; i++) {
	// 		stream >> n;
	// 		costsD.push_back(n);
	// 		costs << n;
	// 	}
	// }

	file.close();
}