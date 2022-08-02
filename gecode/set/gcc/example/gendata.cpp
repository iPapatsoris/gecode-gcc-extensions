#include <iostream>
#include <fstream>
#include <getopt.h>
#include <exception>
#include <algorithm>
#include <random>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>

using namespace std;

enum {
	VARS, VALS, DENSITY, LSUBSET, LOWMIN, LOWMAX, UPMIN, UPMAX, LOWVARMIN,
	LOWVARMAX, UPVARMIN, UPVARMAX, SEED, OUTFILE
};

int main(int argc, char** argv) {
	int opt;
	int totalVars = 0, totalVals = 0, density = 100, valSubset = 100, lowMin = 0, 
		  lowMax = 100, upMin = 0, upMax = 100, lowVarMax = 100, lowVarMin = 0, 
			upVarMax = 100, upVarMin = 0;
	string fileName = "";
	unsigned int seed;
	bool seedFlag = false;
	char *endptr; // dummy
	const int base = 10;
	struct option options[] = {
		{"n", required_argument, NULL, VARS},
		{"m", required_argument, NULL, VALS},
		{"p", required_argument, NULL, DENSITY},
		{"lmin", required_argument, NULL, LOWMIN},
		{"lmax", required_argument, NULL, LOWMAX},
		{"umin", required_argument, NULL, UPMIN},
		{"umax", required_argument, NULL, UPMAX},
		{"lvarmin", required_argument, NULL, LOWVARMIN},
		{"lvarmax", required_argument, NULL, LOWVARMAX},
		{"uvarmin", required_argument, NULL, UPVARMIN},
		{"uvarmax", required_argument, NULL, UPVARMAX},
		{"f", required_argument, NULL, OUTFILE},
		{"s", required_argument, NULL, SEED},
		{0, 0, 0, 0}
	};

	while ((opt = getopt_long_only(argc, argv, "", options, NULL)) != -1) {
		switch (opt) {
			case VARS:
				totalVars = strtol(optarg, &endptr, base);
				break;
			case VALS:
				totalVals = strtol(optarg, &endptr, base);
				break;
			case DENSITY:
				density = strtol(optarg, &endptr, base);
				break;
			case LOWMIN:
				lowMin = strtol(optarg, &endptr, base);
				break;
			case LOWMAX:
				lowMax = strtol(optarg, &endptr, base);
				break;
			case UPMIN:
				upMin = strtol(optarg, &endptr, base);
				break;
			case UPMAX:
				upMax = strtol(optarg, &endptr, base);
				break;
			case LOWVARMIN:
				lowVarMin = strtol(optarg, &endptr, base);
				break;
			case LOWVARMAX:
				lowVarMax = strtol(optarg, &endptr, base);
				break;
			case UPVARMIN:
				upVarMin = strtol(optarg, &endptr, base);
				break;
			case UPVARMAX:
				upVarMax = strtol(optarg, &endptr, base);
				break;
			case SEED:
				seed = strtol(optarg, &endptr, base);
				seedFlag = true;
				break;
			case OUTFILE:
				fileName = string(optarg);
				break;
		}
	}
	
	auto assertPercentage = [](int p) {
		assert(p >= 0 && p <= 100);
	};

	assert(totalVars > 0 && totalVals > 0);
	assertPercentage(density);
	assertPercentage(lowMin);
	assertPercentage(lowMax);
	assertPercentage(upMin);
	assertPercentage(upMax);
	assertPercentage(lowVarMin);
	assertPercentage(lowVarMax);
	assertPercentage(upVarMin);
	assertPercentage(upVarMax);

	// Seed the random engine
	random_device rd;
	if (!seedFlag) {
		seed = rd() ^ (
    	(mt19937::result_type) chrono::duration_cast<chrono::seconds>(
      	chrono::system_clock::now().time_since_epoch()
			).count() +
      (mt19937::result_type) chrono::duration_cast<chrono::microseconds>(
        chrono::high_resolution_clock::now().time_since_epoch()
      ).count() 
		);	
	}
	mt19937 gen(seed);

	if (fileName.empty()) {
		throw invalid_argument("Option -f is missing: no file to write to");
	}
	ofstream file(fileName);

	// Write arguments to file
	file << "# ";
	for (int i = 0; i < argc; i++) {
		file << argv[i] << " ";
	}

	if (!seedFlag) {
		file << "-s " << seed;
	}
	uniform_int_distribution<int> rand(1, 100);
	uniform_int_distribution<int> randVal(0, totalVals - 1);

	// Gen domain
	vector<int> domain(totalVals);
	for (int i = 0; i < totalVals; i++) {
		domain[i] = i;
	}
	shuffle(domain.begin(), domain.end(), gen);

	// The actual values list. Can be a subset of domain
	vector<int> values;
	values.reserve(domain.size() / 100);

	file << "\n# Number of variables\n" << totalVars << "\n\n# Domains\n";

	auto updateValStuctures = [](vector<unordered_set<int>>& varToVals,
															 unordered_map<int, unordered_set<int>>& valToVars
														  ,vector<int>& values, int x, int v) {
		auto it = valToVars.find(v);
			if (it == valToVars.end()) {
				valToVars.insert({v, unordered_set<int>({x})});
				values.push_back(v);
			} else {
				it->second.insert(x);
			}
		varToVals[x].insert(v);
	};

	// Random Erdos Renyi Var->Val graph
	unordered_map<int, unordered_set<int>> valToVars;
	vector<unordered_set<int>> varToVals;
	for (int x = 0; x < totalVars; x++) {
		varToVals.push_back(unordered_set<int>());
		bool emptyDomain = true;
		for (auto v: domain) {
			if (rand(gen) <= density) {
				file << (emptyDomain ? "" : ",") << v;
				emptyDomain = false;
				updateValStuctures(varToVals, valToVars, values, x, v);
			}
		}
		if (emptyDomain) {
			// Make sure variable has at least 1 value in domain
			int val = values[randVal(gen)];
			file << val;
			updateValStuctures(varToVals, valToVars, values, x, val);
		}
		file << "\n";
	}

	if (valSubset != 100) {
		values.resize(valSubset * values.size() / 100);
	}

	file << "\n# Values\n";
	for (auto it = values.begin(); it != values.end(); it++) {
		file << *it << (it != values.end()-1 ? "," : "\n");
	}

	vector<unsigned int> upperBounds(values.size());
	vector<unsigned int> lowerBounds(values.size());

	for (int i = 0; i < values.size(); i++) {
		int edges = valToVars.find(values[i])->second.size();
		int min = upMin * edges / 100;
		if (!min) {
			min = 1;
		}
		int max = upMax * edges / 100;
		if (!max) {
			max = 1;
		}
		assert(min <= max);
		uniform_int_distribution<int> randUpperBound(min, max);
		upperBounds[i] = randUpperBound(gen);

		min = lowMin * edges / 100;
		max = lowMax * edges / 100;
		assert(min <= max);
		uniform_int_distribution<int> randLowerBound(min, max);
		lowerBounds[i] = randLowerBound(gen);
	}

	vector<unsigned int> upperVarBounds(totalVars);
	vector<unsigned int> lowerVarBounds(totalVars);

	for (int i = 0; i < totalVars; i++) {
		int edges = varToVals[i].size();
		int min = upVarMin * edges / 100;
		int max = upVarMax * edges / 100;
		assert(min <= max);
		uniform_int_distribution<int> randUpperBound(min, max);
		upperVarBounds[i] = randUpperBound(gen);

		min = lowVarMin * edges / 100;
		max = lowVarMax * edges / 100;
		assert(min <= max);
		uniform_int_distribution<int> randLowerBound(min, max);
		lowerVarBounds[i] = randLowerBound(gen);
	}

	file << "\n# Lower bounds\n";
	for (auto it = lowerBounds.begin(); it != lowerBounds.end(); it++) {
		file << *it << (it != lowerBounds.end()-1 ? "," : "\n");
	}

	file << "\n# Upper bounds\n";
	for (auto it = upperBounds.begin(); it != upperBounds.end(); it++) {
		file << *it << (it != upperBounds.end()-1 ? "," : "\n");
	}
	
	file << "\n# Lower var bounds\n";
	for (auto it = lowerVarBounds.begin(); it != lowerVarBounds.end(); it++) {
		file << *it << (it != lowerVarBounds.end()-1 ? "," : "\n");
	}

	file << "\n# Upper var bounds\n";
	for (auto it = upperVarBounds.begin(); it != upperVarBounds.end(); it++) {
		file << *it << (it != upperVarBounds.end()-1 ? "," : "\n");
	}

	return 0;
}