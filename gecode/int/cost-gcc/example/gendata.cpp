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
	VARS, VALS, DENSITY, LOWMIN, LOWMAX, UPMIN, UPMAX, CMIN, CMAX, 
	COST, SEED, OUTFILE
};

int main(int argc, char** argv) {
	int opt;
	int totalVars = 0, totalVals = 0, density = 100, lowMin = 0, 
		  lowMax = 100, upMin = 0, upMax = 100, costMin = 1, costMax = 10, cost = 5;
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
		{"cmin", required_argument, NULL, CMIN},
		{"cmax", required_argument, NULL, CMAX},
		{"c", required_argument, NULL, COST},
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
			case CMIN:
				costMin = strtol(optarg, &endptr, base);
				break;
			case CMAX:
				costMax = strtol(optarg, &endptr, base);
				break;
			case COST:
				cost = strtol(optarg, &endptr, base);
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
	assert (costMin <= costMax);
	assertPercentage(density);
	assertPercentage(lowMin);
	assertPercentage(lowMax);
	assertPercentage(upMin);
	assertPercentage(upMax);

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

	// The actual values list
	vector<int> values;
	values.reserve(domain.size() / 100);

	file << "\n# Number of variables\n" << totalVars << "\n\n# Cost upper bound\n"
			 << cost << "\n\n# Domains\n";

	auto updateValStuctures = [](unordered_map<int, unordered_set<int>>& valToVars
														  ,vector<int>& values, int x, int v) {
		auto it = valToVars.find(v);
			if (it == valToVars.end()) {
				valToVars.insert({v, unordered_set<int>({x})});
				values.push_back(v);
			} else {
				it->second.insert(x);
			}
	};

	// Random Erdos Renyi Var->Val graph
	unordered_map<int, unordered_set<int>> valToVars;
	for (int x = 0; x < totalVars; x++) {
		bool emptyDomain = true;
		for (auto v: domain) {
			if (rand(gen) <= density) {
				file << (emptyDomain ? "" : ",") << v;
				emptyDomain = false;
				updateValStuctures(valToVars, values, x, v);
			}
		}
		if (emptyDomain) {
			// Make sure variable has at least 1 value in domain
			int val = values[randVal(gen)];
			file << val;
			updateValStuctures(valToVars, values, x, val);
		}
		file << "\n";
	}

	file << "\n# Values\n";
	for (auto it = values.begin(); it != values.end(); it++) {
		file << *it << (it != values.end()-1 ? "," : "\n");
	}

	vector<unsigned int> upperBounds(values.size());
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
	}

	file << "\n# Lower bounds\n";
	for (int i = 0; i < values.size(); i++) {
		int min = lowMin * upperBounds[i] / 100;
		int max = lowMax * upperBounds[i] / 100;
		assert(min <= max);
		uniform_int_distribution<int> randLowerBound(min, max);
		file << randLowerBound(gen) << (i == values.size() - 1 ? "\n" : ",");
	}

	file << "\n# Upper bounds\n";
	for (auto it = upperBounds.begin(); it != upperBounds.end(); it++) {
		file << *it << (it != upperBounds.end()-1 ? "," : "\n");
	}

	uniform_int_distribution<int> randCost(costMin, costMax);
	file << "\n# Costs\n";
	for (int x = 0; x < totalVars; x++) {
		for (auto it = values.begin(); it != values.end(); it++) {
			// Not efficient, can fix by including a varToVals structure
			auto vars = valToVars.find(*it)->second;
			int cost = (vars.find(x) != vars.end() ? randCost(gen) : 0);
			file << cost << (it == values.end()-1 ? "\n" : ",");
		}
	}

	return 0;
}