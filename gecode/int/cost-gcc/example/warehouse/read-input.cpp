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

enum {
	VARS, COST, DOMAIN, VALS, LOWER_BOUNDS, UPPER_BOUNDS, COSTS
} mode;

void readInput(string fileName, int& vars, IntSetArgs& domain, IntArgs& vals,
							 IntArgs& lowerBounds, IntArgs& upperBounds, IntArgs& costs, int& fixed, IntArgs& demands) {

	string line;
  ifstream file(fileName);
	if (!file.is_open()) {
		throw "Could not open file";
	}

	getline(file, line);
	getline(file, line);
	getline(file, line);
	getline(file, line);

	struct Assignment {
		int cost;
		int demand;

		Assignment(int cost, int demand) : cost(cost), demand(demand) {}
		Assignment() {}
	};

	vector<int> numbers;
	int n;
	vector<unordered_map<unsigned int, Assignment>> varToVals;
	stringstream stream(line);

	stream >> vars;
	stream >> fixed;
	stream >> n;
	for (unsigned int i = 1; i <= vars; i++) {
		vals << i;
		lowerBounds << 0;
		upperBounds << n; 
		varToVals.push_back(unordered_map<unsigned int, Assignment>());
	}

	getline(file, line);
	

	while (getline(file, line)) {
		stringstream stream(line);
		unsigned int warehouse;
		unsigned int store;
		int cost;
		double costD;
		int demand;
		stream >> warehouse;
		stream >> store;
		stream >> costD;
		cost = costD;
		stream >> demand;
		varToVals[store-1].insert(pair<unsigned int, Assignment>(warehouse, Assignment(cost, demand)));
	}


	for (unsigned int var = 0; var < vars; var++) {
		vector<int> keys;
		for (auto k: varToVals[var]) {
			keys.push_back(k.first);
		}
		IntSet tmp;
		IntSetInit<IntArgs>::init(tmp, IntArgs(keys));
		domain << tmp;
		for (unsigned int val = 1; val <= vals.size(); val++) {
			auto res = varToVals[var].find(val);
			costs << (res != varToVals[var].end() ? res->second.cost : 0); 
			demands << (res != varToVals[var].end() ? res->second.demand : 0); 
		}
	}

	for (unsigned int i = 0; i < vars; i++) {
		cout << "{";
		for (unsigned int j = 0; j < vals.size(); j++) {
			cout << demands[i*vals.size() + j] << (j == vals.size()-1 ? "" :",");
		}
		cout << "},\n";
	}
	exit(1);


	file.close();
	return;

}