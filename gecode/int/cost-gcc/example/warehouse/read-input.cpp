#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <gecode/int.hh>

using namespace std;
using namespace Gecode;

enum {
	VARS, COST, DOMAIN, VALS, LOWER_BOUNDS, UPPER_BOUNDS, COSTS
} mode;

void readInput(string fileName, int& vars, IntSetArgs& domain, IntArgs& vals,
							 IntArgs& lowerBounds, IntArgs& upperBounds, IntArgs& costs, 
							 int& cost) {
// 	string line;
//   ifstream file(fileName);
// 	if (!file.is_open()) {
// 		throw "Could not open file";
// 	}

// 	getline(file, line);
// 	stringstream stream(line);
// 	vector<int> numbers;
// 	int n;
// 		while(stream >> n) {
// 			numbers.push_back(n);
// 		}

	
// 	vars = numbers.size();
// 	for (unsigned int i = 0; i < vars; i++) {
// 		vector<int> d;
// 		for (unsigned int j = 0; j < vars; j++) {
// 			if (i != j) {
// 				d.push_back(j);
// 			}
// 		}
// 		IntSet tmp;
// 		IntSetInit<IntArgs>::init(tmp, IntArgs(d));
// 		domain << tmp;
// 		vals << i;
// 		lowerBounds << 1;
// 		upperBounds << 1;
// 	}

// 	getline(file, line);
// 	stream = stringstream(line);
// 	unsigned int i = 0;
// 	vector<int> numbers2;
// 	while(stream >> n) {
// 		numbers2.push_back(n);
// 	}

// 	for (auto n1: numbers) {
// 		for (auto n2: numbers2) {
// 			costs << n1 * n2;
// 		}
// 	}

// 	costs << numbers;
// 	cost = 4330;

// 	file.close();


	string line;
  ifstream file(fileName);
	if (!file.is_open()) {
		throw "Could not open file";
	}

	getline(file, line);
	stringstream stream(line);
	vector<int> numbers;
	int n;
	stream >> n;
	cout << n << endl;

	getline(file, line);
	stream = stringstream(line);
	stream >> n;
	cout << n << endl;

	return;

	vars = n;
	for (unsigned int i = 0; i < vars; i++) {
		vector<int> d;
		for (unsigned int j = 0; j < vars; j++) {
			if (i != j) {
				d.push_back(j);
			}
		}
		IntSet tmp;
		IntSetInit<IntArgs>::init(tmp, IntArgs(d));
		domain << tmp;
		vals << i;
		lowerBounds << 1;
		upperBounds << 1;
	}

	getline(file, line);
	stream = stringstream(line);
	stream >> cost;

	getline(file, line);
	char c;
	stream = stringstream(line);
	while(stream >> n) {
		costs << n;
		stream >> c;
	}

	file.close();

}