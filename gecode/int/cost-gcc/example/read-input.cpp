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
	string line;
  ifstream file(fileName);
	if (!file.is_open()) {
		throw "Could not open file";
	}

	mode = VARS;
	int var = 0;
  while(getline(file, line)) {
		vector<int> numbers;
		stringstream stream(line);
		int n;
		char c;
		while(stream >> n) {
			// Ignore characters other than numbers
			stream >> c;
			numbers.push_back(n);
		}
		if (!numbers.empty()) {
			switch(mode) {
				case VARS:
					vars = numbers.front();
					mode = COST;
					break;
				case COST:
					cost = numbers.front();
					mode = DOMAIN;
					break;
				case DOMAIN: {
					IntSet tmp;
					IntSetInit<IntArgs>::init(tmp, IntArgs(numbers));
					domain << tmp;
					if (domain.size() == vars) {
						mode = VALS;
					}
					break;
				}
				case VALS:
					vals = numbers;
					mode = LOWER_BOUNDS;
					break;
				case LOWER_BOUNDS:
					lowerBounds = numbers;
					mode = UPPER_BOUNDS;
					break;
				case UPPER_BOUNDS:
					upperBounds = numbers;
					mode = COSTS;
					break;
				case COSTS:
					costs << numbers;
					break;
			}
		}
  }

	file.close();
}