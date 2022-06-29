#ifndef H_COST_GCC_POST
#define H_COST_GCC_POST

#include <gecode/int.hh>
#include <unordered_map>
#include <unordered_set>
#include "cost-gcc.hpp"

using namespace Gecode;
using namespace std;

void countCosts(Space& home, const IntVarArgs& vars, const IntArgs& vals,
								const IntArgs& lowerBounds, const IntArgs& upperBounds,
								const IntArgs& costs, IntVar costUpperBound, BestBranch* bestBranch = NULL,
								IntPropLevel ipl = IPL_DOM) {
								
	using namespace Int;

	int n = vars.size(); 
	if (!n) {
		throw TooFewArguments("Int::countCosts");
	}
	if (same(vars)) {
		throw ArgumentSame("Int::countCosts");
	}

	// Bounds arrays must be same size as vals array
	if (vals.size() != lowerBounds.size() || vals.size() != upperBounds.size()) {
		throw ArgumentSizeMismatch("Int::countCosts");
	}

	// Don't allow duplicates in vals
	// Map each value to its respective position in vals array
	unordered_map<int, int> valToIndex;
	for (int i = 0; i < vals.size(); i++) {
		if (valToIndex.insert({vals[i], i}).second == false) {
			throw ArgumentSizeMismatch("Int::countCosts");
		}
	}

	// Map values to their variables
	// Helps do fast lookups, for early pruning and FlowGraph creation
	MapToSet valToVars;
	// Map variables to their values
	// Is used to compare old domain with current Gecode domain, to find which
	// values got pruned between executions
	vector<unordered_set<int> > varToVals;
	for (int x = 0; x < vars.size(); x++) {
		varToVals.push_back(unordered_set<int>());
		for (IntVarValues i(vars[x]); i(); ++i) {
			if (valToIndex.find(i.val()) == valToIndex.end()) {
				throw ArgumentSizeMismatch("Int::costcGCC domain value doesn't "
																	 "exist in values array");
			}
			varToVals.back().insert(i.val());
			auto it = valToVars.find(i.val());
			if (it == valToVars.end()) {
				valToVars.insert({i.val(), 
														  unordered_set<int>({x})});
			} else {
				it->second.insert(x);
			}
		}
	}

	// Bounds must be nonnegative and lowerBound smaller or equal to upperBound
	// upperBound must also be at least 1
	for (auto i = 0; i < vals.size(); i++) {
		Int::Limits::nonnegative(lowerBounds[i], "Int::countCosts");
		Int::Limits::nonnegative(upperBounds[i], "Int::countCosts");
		if (upperBounds[i] < lowerBounds[i] || upperBounds[i] < 1) {
			throw OutOfLimits("Int::countCosts");
		}
	}

	if (costs.size() != n * vals.size()) {
		throw ArgumentSizeMismatch("Int::countCosts");
	}

	if (costUpperBound.max() < 0) {
		home.fail();
		return;
	}

	int minCost = 0;
	for (auto cost: costs) {
		if (cost < minCost) {
			minCost = cost;
		}
	}
	
	// If there are negative costs, transform them and costUpperBound
	IntVar positiveCostUpperBound(home, 0, Int::Limits::max);
	IntArgs positiveCosts;
	if (minCost < 0) {
		minCost *= -1;
		rel(home, positiveCostUpperBound == costUpperBound + vars.size() * minCost);
		for (auto cost: costs) {
			positiveCosts << cost + minCost;
		}
	} else {
		positiveCostUpperBound = costUpperBound;
		positiveCosts = costs;
	}

	ViewArray<Int::IntView> views(home, vars);
	GECODE_POST;
	GECODE_ES_FAIL(CostGcc::post(home, views, varToVals, valToVars, vals, 
															 lowerBounds, upperBounds, positiveCosts, 
															 positiveCostUpperBound, bestBranch, ipl
															));
}

#endif