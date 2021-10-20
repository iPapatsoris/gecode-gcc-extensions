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
	VARS, DOMAIN, VALS, LOWER_BOUNDS, UPPER_BOUNDS, LOWER_VAR_BOUNDS, 
	UPPER_VAR_BOUNDS
} mode;

void readInput(string fileName, int& vars, IntSetArgs& domain, IntArgs& vals,
							 IntArgs& lowerBounds, IntArgs& upperBounds, 
							 IntArgs& lowerVarBounds, IntArgs& upperVarBounds) {
	string line;
  ifstream file(fileName);
	if (!file.is_open()) {
		throw "Could not open file";
	}

	mode = VARS;
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
					mode = LOWER_VAR_BOUNDS;
					break;
				case LOWER_VAR_BOUNDS:
					lowerVarBounds = numbers;
					mode = UPPER_VAR_BOUNDS;
					break;
				case UPPER_VAR_BOUNDS:
					upperVarBounds = numbers;
					break;
			}
		}
  }

	file.close();
}