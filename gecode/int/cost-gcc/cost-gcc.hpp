#include <gecode/int.hh>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include "flow-graph-algorithms.hpp"

using namespace Gecode;
using namespace std;

class CostGcc : public NaryPropagator<Int::IntView, Int::PC_INT_DOM> {
protected:

public:
	CostGcc(Space& home, ViewArray<Int::IntView> x)
			: NaryPropagator(home, x) {}

	static ExecStatus post(Space& home, ViewArray<Int::IntView>& vars,
												MapToSet& valToVars,
												const IntArgs& inputVals, 
												const unordered_map<int, unsigned int>& inputValToIndex,
												const IntArgs& lowerBounds, const IntArgs& upperBounds,
												const IntArgs& costs, int costUpperBound) {

		if (pruneOmittedVales(home, vars, valToVars, inputValToIndex) == ES_FAILED) {
			return ES_FAILED;
		}
		if (assignUsingLowerBounds(home, vars, valToVars, inputVals, 
															 inputValToIndex, lowerBounds) == ES_FAILED) {
			return ES_FAILED;
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

		FlowGraph graph(vars, valToVars, inputVals, lowerBounds, upperBounds, costs,
										true);
		#ifndef NDEBUG
			graph.print();
		#endif

		FlowGraphAlgorithms graphAlgorithms(graph);
		if (!graphAlgorithms.findMinCostFlow()) {
			return ES_FAILED;
		}
		if (graph.getFlowCost() > costUpperBound) {
			cout << "Cost constraint failed!" << endl;
			return ES_FAILED;
		}

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

	// Prune values that belong to a domain, but are not mentioned in the 
	// inputVals array/set
	ExecStatus static pruneOmittedVales(
											Space& home, ViewArray<Int::IntView>& vars, 
											MapToSet& valToVars,
											const unordered_map<int, unsigned int>& inputValToIndex) {
		
		unordered_set<int> prunedVals;
		for (auto& v: valToVars) {
			auto value = v.first;
			if (inputValToIndex.find(value) == inputValToIndex.end()) {
				for (auto& x: v.second) {
					GECODE_ME_CHECK(vars[x].nq(home, value));
					cout << "Prunning " << value << " from " << x << 
									" (not mentioned in values array)\n";
				}
				prunedVals.insert(value);
			}
		}

		// Also mirror the changes to valToVars
		for (auto& val: prunedVals) {
			valToVars.erase(val);
		}

		return ES_OK;
	}

	// Check lower bounds for early propagation:
	// - A value that has lower bound greater than the variables that can it can
	//   be assigned to, means that the bounds restriction will always fail
	// - A value that has lower bound equal to the variables than it can be
	//   assigned to, means that these variables can be assigned with it already
	ExecStatus static assignUsingLowerBounds(
												Space& home, ViewArray<Int::IntView>& vars, 
												MapToSet& valToVars, const IntArgs& inputVals,
												const unordered_map<int, unsigned int>& inputValToIndex,
												const IntArgs& lowerBounds) {

		// Assigning a value to a variable means that we prune the rest of its 
		// values from its domain. Subsequently, those values can be assigned to
		// less variables than before. Thus, we can repeat our checks, until
		// no more assignments are possible.
		// 
		// We hold a pair of prunedVals structure, and pointers to them, to be 
		// able to swap reading/writing from/to them, during the fixpoint loop.
		unordered_set<int> prunedVals1, prunedVals2;
		unordered_set<int>* prunedValsRead, *prunedValsWrite;
		prunedValsWrite = &prunedVals1;
		prunedValsRead = &prunedVals2;

		// Initial check for assignments
		for (int i = 0; i < inputVals.size(); i++) {
			auto value = inputVals[i];
			auto status = checkLowerBounds(home, vars, value, lowerBounds[i],
																		 valToVars, *prunedValsWrite,
																		 *prunedValsRead);
			if (status == ES_FAILED) {
				return status;
			}
		}
		// Keep checking until we reach a fixpoint
		while (!prunedValsWrite->empty()) {
			swap(prunedValsRead, prunedValsWrite);
			prunedValsWrite->clear();
			for (auto value: *prunedValsRead) {
				auto itValToIndex = inputValToIndex.find(value);
				assert(itValToIndex != inputValToIndex.end());
				unsigned int i = itValToIndex->second;
				auto status = checkLowerBounds(home, vars, value, lowerBounds[i], 
																			 valToVars, *prunedValsWrite, 
																			 *prunedValsRead);
				if (status == ES_FAILED) {
					return status;
				}
			}
		}

		return ES_OK;
	}

	// Check the lower bound for specific value, according to the rules mentioned
	// in assignUsingLowerBounds() function. If a value has been pruned early
	// from all possible variables, we fail if its lower bound is greater than 0.
	ExecStatus static checkLowerBounds(Space& home, ViewArray<Int::IntView>& vars, 
																	 	 int value, unsigned int lowerBound, 
																		 MapToSet& valToVars,
																	 	 unordered_set<int>& prunedVals, 
																		 const unordered_set<int>& alreadyProcessing
																		 ) {
		auto it = valToVars.find(value);
		unsigned int totalVarsWithThisVal = (it != valToVars.end() ? 
																				 it->second.size() : 0);
		if (lowerBound > totalVarsWithThisVal) {
			cout << "Lower bound of value" << value << " is greater than remaining " 
					<< "edges: " << totalVarsWithThisVal << endl;
			return ES_FAILED;
		} else if (it != valToVars.end() && lowerBound == totalVarsWithThisVal) {
			for (auto x: it->second) {
				cout << "Assigning var " << x << " to " << value << endl;
				assignValToVar(x, vars[x], value, valToVars, prunedVals, 
												alreadyProcessing);
				GECODE_ME_CHECK(vars[x].eq(home, value));
			}
		}
	return ES_OK;
	}

	// When a value is assigned to a variable: for every other value in its 
	// domain, we need to find its corresponding variables set in valToVars,
	// and remove the assigned variable from it. 
	// Also put those values in prunedValues, unless they exist in 
	// alreadyProcessing set. If they do, we are already iterating through them
	void static assignValToVar(int xIndex, Int::IntView x, int val, 
											       MapToSet& valToVars, 
														 unordered_set<int>& prunedVals, 
														 const unordered_set<int>& alreadyProcessing) {
		for (IntVarValues v(x); v(); ++v) {
			if (v.val() != val) {
				if (alreadyProcessing.find(v.val()) == alreadyProcessing.end()) {
					prunedVals.insert(v.val());
				}
				auto otherValVars = valToVars.find(v.val());
				assert(otherValVars != valToVars.end());
				otherValVars->second.erase(xIndex);
				if (otherValVars->second.empty()) {
					valToVars.erase(otherValVars);
				}
			}
		}
	}

	// Assert that Gecode variable domains and valToVars are in sync
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
