#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <gecode/int.hh>
#include "read-input.hpp"

using namespace std;
using namespace Gecode;

enum {
	VARS, DOMAIN, VALS, LOWER_BOUNDS, UPPER_BOUNDS, LOWER_VAR_BOUNDS, 
	UPPER_VAR_BOUNDS
} mode;

void readInput(string fileName, int& vars, vector<unordered_set<int> >& domain, IntArgs& vals,
							 IntArgs& lowerValBounds, IntArgs& upperValBounds, 
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
					unordered_set<int> d(numbers.begin(), numbers.end());
					
					domain.push_back(d);
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
					lowerValBounds = numbers;
					mode = UPPER_BOUNDS;
					break;
				case UPPER_BOUNDS:
					upperValBounds = numbers;
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