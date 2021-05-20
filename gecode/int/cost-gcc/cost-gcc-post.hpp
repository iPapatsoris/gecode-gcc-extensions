#include <gecode/int.hh>
#include <unordered_map>
#include <unordered_set>
#include "cost-gcc.hpp"

using namespace Gecode;
using namespace std;

void countCosts(Space& home, const IntVarArgs& vars, const IntArgs& vals,
								const IntArgs& lowerBounds, const IntArgs& upperBounds,
								const IntArgs& costs, int costUpperBound, 
								IntPropLevel ipl) {

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

	// No duplicates in vals
	unordered_set<int> distinctVals;
	for (auto& val: vals) {
		if (distinctVals.insert(val).second == false) {
			throw ArgumentSizeMismatch("Int::countCosts");
		}
	}

	// Map values to their variables
	MapToSet valToVars;
	for (int x = 0; x < vars.size(); x++) {
		for (IntVarValues i(vars[x]); i(); ++i) {
			auto it = valToVars.find(i.val());
			if (it == valToVars.end()) {
				valToVars.insert({i.val(), unordered_set<int>({x})});
			} else {
				it->second.insert(x);
			}
		}
	}

	// Values in vals argument must belong to at least one variable domain
	for (auto& v: distinctVals) {
		if (valToVars.find(v) == valToVars.end()) {
			throw OutOfLimits("Int::countCosts");
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
	
	ViewArray<Int::IntView> views(home, vars);
	GECODE_POST;
	GECODE_ES_FAIL(CostGcc::post(home, views, valToVars, vals, distinctVals, 
															 lowerBounds, upperBounds, costs, costUpperBound
															));
}