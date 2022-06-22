#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <gecode/int.hh>

#include <unordered_map>
#include <unordered_set>
using namespace std;
using namespace Gecode;


void readInput(string fileName, int& vars, IntSetArgs& domain, IntArgs& vals,
							 IntArgs& lowerBounds, IntArgs& upperBounds, IntArgs& costs, int& fixed) {

	string line;
  ifstream file(fileName);
	if (!file.is_open()) {
		throw "Could not open file";
	}

	getline(file, line);
	getline(file, line);
	getline(file, line);
	vector<int> numbers;
	vector<unordered_map<int, int>> varToVals;
	stringstream stream(line);
	stream >> vars;
	stream >> fixed;
	for (int i = 1; i <= vars; i++) {
		vals << i;
		lowerBounds << 0;
		upperBounds << vars; 
		varToVals.push_back(unordered_map<int, int>());
	}

	getline(file, line);
	

	while (getline(file, line)) {
		stringstream stream(line);
		int warehouse;
		int store;
		int cost;
		stream >> warehouse;
		stream >> store;
		stream >> cost;
		varToVals[store-1].insert(pair<int, int>(warehouse, cost));
	}

	for (int var = 0; var < vars; var++) {
		vector<int> keys;
		for (auto k: varToVals[var]) {
			keys.push_back(k.first);
		}
		IntSet tmp;
		IntSetInit<IntArgs>::init(tmp, IntArgs(keys));
		domain << tmp;
		for (int val = 1; val <= vals.size(); val++) {
			auto res = varToVals[var].find(val);
			costs << (res != varToVals[var].end() ? res->second : 0); 
		}
	}

	file.close();
	return;

}