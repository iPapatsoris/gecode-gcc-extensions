#ifndef H_SYM_GCC_POST
#define H_SYM_GCC_POST

#include <gecode/int.hh>
#include <unordered_map>
#include <unordered_set>
#include "util.hpp"
#include "sym-gcc.hpp"
#include "example/LI.hpp"

using namespace Gecode;
using namespace std;



void symmetricGCC(Space& home, const BoolVarArgs& x, const vector<unordered_set<int> > domains, const IntArgs& vals,
								const IntArgs& lowerValBounds, const IntArgs& upperValBounds,
								const IntArgs& lowerVarBounds, const IntArgs& upperVarBounds, 
								const VarUtil& varUtil,
								LI* li, IntPropLevel ipl) {
								
	using namespace Int;

	int n = x.size(); 
	if (!n) {
		throw TooFewArguments("Int::symmetricGCC");
	}
	if (same(x)) {
		throw ArgumentSame("Int::symmetricGCC");
	}

	// Val bounds arrays must be same size as vals array
	if (vals.size() != lowerValBounds.size() || vals.size() != upperValBounds.size()) {
		throw ArgumentSizeMismatch("Int::symmetricGCC");
	}

	// Var bounds arrays must be same size as vars array
	if (varUtil.inputVarToXIndex.size() != lowerVarBounds.size() || varUtil.inputVarToXIndex.size() != upperVarBounds.size()) {
		throw ArgumentSizeMismatch("Int::symmetricGCC");
	}

	// Don't allow duplicates in vals
	// Map each value to its respective position in vals array
	unordered_map<int, unsigned int> valToIndex;
	for (int i = 0; i < vals.size(); i++) {
		if (valToIndex.insert({vals[i], i}).second == false) {
			throw ArgumentSizeMismatch("Int::symmetricGCC");
		}
	}

	unordered_set<int> valsSet;
	for (auto v: vals) {
		valsSet.insert(v);
	}

	// Map variables to their values
	// Is used to compare old domain with current Gecode domain, to find which
	// values got pruned between executions
	MapToSet<int, unsigned int> valToVars;
	for (int i = 0; i < domains.size(); i++) {
		for (auto v: domains[i]) {
			if (valsSet.find(v) == valsSet.end()) {
				throw ArgumentSizeMismatch("Int::symmetricGCC domain value doesn't exist in values array");
			}
			auto it = valToVars.map.find(v);
			if (it == valToVars.map.end()) {
				valToVars.map.insert({v, 
														  unordered_set<unsigned int>({(unsigned int) i})});
			} else {
				it->second.insert(i);
			}
		}
	}

	// Bounds must be nonnegative and lowerBound smaller or equal to upperBound
	// upperBound must also be at least 1
	for (auto i = 0; i < vals.size(); i++) {
		Int::Limits::nonnegative(lowerValBounds[i], "Int::symmetricGCC");
		Int::Limits::nonnegative(upperValBounds[i], "Int::symmetricGCC");
		if (upperValBounds[i] < lowerValBounds[i] || upperValBounds[i] < 1) {
			throw OutOfLimits("Int::symmetricGCC");
		}
	}

	// Variable bounds must be nonnegative and lowerBound smaller or equal to 
	// upperBounb
	for (auto i = 0; i < varUtil.inputVarToXIndex.size(); i++) {
		Int::Limits::nonnegative(lowerVarBounds[i], "Int::symmetricGCC");
		Int::Limits::nonnegative(upperVarBounds[i], "Int::symmetricGCC");
		if (upperVarBounds[i] < lowerVarBounds[i]) {
			throw OutOfLimits("Int::symmetricGCC");
		}
	}
	
	ViewArray<Int::BoolView> views(home, x);


	GECODE_POST;
	GECODE_ES_FAIL(SymGcc::post(home, views, valToVars, domains, vals, lowerValBounds, 
															upperValBounds, lowerVarBounds, upperVarBounds, 
															varUtil,
															li, ipl
															));
}

#endif