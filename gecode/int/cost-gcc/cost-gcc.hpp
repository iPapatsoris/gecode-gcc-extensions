#include <gecode/int.hh>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>

using namespace Gecode;
using namespace std;

typedef unordered_map<int, unordered_set<int>> MapToSet;

class CostGcc : public NaryPropagator<Int::IntView, Int::PC_INT_DOM> {
protected:

public:
	CostGcc(Space& home, ViewArray<Int::IntView> x)
			: NaryPropagator(home, x) {}

	static ExecStatus post(Space& home, ViewArray<Int::IntView>& vars,
												 MapToSet& valToVars,
												 const IntArgs& inputVals, 
												 const unordered_set<int>& inputValsSet,
												 const IntArgs& lowerBounds, const IntArgs& upperBounds,
												 const IntArgs& costs, int costUpperBound) {

		// Prune values that are not mentioned in the inputVals array/set
		unordered_set<int> prunedVals;
		for (auto& v: valToVars) {
			auto value = v.first;
			if (inputValsSet.find(value) == inputValsSet.end()) {
				for (auto& x: v.second) {
					GECODE_ME_CHECK(vars[x].nq(home, value));
					cout << "Prunning " << value << " from " << x << 
									" (not mentioned in values array)\n";
				}
				prunedVals.insert(value);
			}
		}

		// Also remove them from valToVars
		for (auto& val: prunedVals) {
			valToVars.erase(val);
		}

		// Check lower bounds for early propagation:
		// - A value that has lower bound greater than the variables that can it can
		//   be assigned to, means that the bounds restriction will always fail
		// - A value that has lower bound equal to the variables than it can be 
		//   assigned to, means that these variables can be assigned with it already
		for (int i = 0; i < inputVals.size(); i++) {
			auto value = inputVals[i];
			auto it = valToVars.find(value);
			assert(it != valToVars.end());
			int totalVarsWithThisVal = it->second.size();
			if (lowerBounds[i] > totalVarsWithThisVal) {
				cout << "Lower bound of value" << value << " is greater than remaining " 
						 << "edges: " << totalVarsWithThisVal << endl;
				return ES_FAILED;
			} else if (lowerBounds[i] == totalVarsWithThisVal) {
				for (auto x: it->second) {
					cout << "Assigning var " << x << " to " << value << endl;
					assignValToVar(x, vars[x], value, valToVars);
					GECODE_ME_CHECK(vars[x].eq(home, value));
				}
			}
		}

		#ifndef NDEBUG
			cout << vars << endl;
			for (auto &val : valToVars) {
				cout << "Val " << val.first << " vars: ";
				for (auto var : val.second) {
					cout << var << " ";
				}
				cout << "\n";
			}
			assertCorrectDomains(vars, valToVars);
		#endif

		(void)new (home) CostGcc(home, vars);
		return ES_OK;
	}

	CostGcc(Space& home, CostGcc& p)
			: NaryPropagator<Int::IntView, Int::PC_INT_DOM>(home, p) {}

	virtual Propagator *copy(Space& home) {
		return new (home) CostGcc(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::cubic(PropCost::LO, x.size());
	}

	virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
		return ES_NOFIX;
	}

private:
	// When a value is assigned to a variable: for every other value in its 
	// domain, we need to find its corresponding variables set in valToVars,
	// and remove the assigned variable from it
	void static assignValToVar(int xIndex, Int::IntView x, int val, 
											       MapToSet& valToVars) {
		for (IntVarValues v(x); v(); ++v) {
			if (v.val() != val) {
				auto otherValVars = valToVars.find(v.val());
				assert(otherValVars != valToVars.end());
				otherValVars->second.erase(xIndex);
				if (otherValVars->second.empty()) {
					valToVars.erase(otherValVars);
				}
			}
		}
	}

	// Make sure Gecode variable domains and valToVars are in sync
	void static assertCorrectDomains(const ViewArray<Int::IntView>& vars, 
															 		 const MapToSet& valToVars) {
		for (int x = 0; x < vars.size(); x++) {
			for (IntVarValues v(vars[x]); v(); ++v) {
				auto it = valToVars.find(v.val());
				assert(it != valToVars.end());
				assert(it->second.find(x) != it->second.end());
			} 
		}

		for (auto& v: valToVars) {
			for (auto x: v.second) {
				assert(x < vars.size());
				assert(vars[x].in(v.first));
			}
		}
	}
};
