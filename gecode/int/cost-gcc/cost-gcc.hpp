#include <gecode/int.hh>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include "flow-graph-algorithms.hpp"
//#include "flow-graph-seq.hpp"

using namespace Gecode;
using namespace std;

typedef NaryPropagator<Int::IntView, Int::PC_INT_NONE> CostGccBase;

class CostGcc : public CostGccBase {

protected:
		class ViewAdvisor : public Advisor {
			public:
				Int::IntView x;
				unsigned int xIndex;
				ViewAdvisor(Space& home, Propagator& p, Council<ViewAdvisor>& c, 
										Int::IntView x, unsigned int xIndex) 
					: Advisor(home, p, c), x(x), xIndex(xIndex) {
					x.subscribe(home, *this);
				}
				ViewAdvisor(Space& home, ViewAdvisor& a)
					: Advisor(home, a) {
					x.update(home, a.x);
					xIndex = a.xIndex;
				}
				void dispose(Space& home, Council<ViewAdvisor>& c) {
					x.cancel(home, *this);
					Advisor::dispose(home, c);
				}
		};
		Council<ViewAdvisor> c;
		FlowGraph* graph;
		vector<EdgeNodes> updatedEdges;
		// TODO: do not store, instead use different post functions?
		IntPropLevel ipl;

public:
	CostGcc(Space& home, ViewArray<Int::IntView> x, FlowGraph* graph, 
					const vector<EdgeNodes>& updatedEdges, IntPropLevel ipl)
			: NaryPropagator(home, x), c(home), graph(graph), 
				updatedEdges(updatedEdges), ipl(ipl) {
		for (int i = 0; i < x.size(); i++) {
			(void)new (home) ViewAdvisor(home, *this, c, x[i], i);
		}
		home.notice(*this, AP_DISPOSE);
	}

	static ExecStatus post(Space& home, ViewArray<Int::IntView>& vars,
												MapToSet<unsigned int, int>& varToVals,
												MapToSet<int, unsigned int>& valToVars,
												const IntArgs& inputVals, 
												const unordered_map<int, unsigned int>& inputValToIndex,
												const IntArgs& lowerBounds, const IntArgs& upperBounds,
												const IntArgs& costs, int costUpperBound,
												IntPropLevel ipl) {

		if (pruneOmittedVales(home, vars, varToVals, valToVars, 
													inputValToIndex) == ES_FAILED) {
			return ES_FAILED;
		}
		if (assignUsingLowerBounds(home, vars, varToVals, valToVars, inputVals, 
															 inputValToIndex, lowerBounds) == ES_FAILED) {
			return ES_FAILED;
		}

		#ifndef NDEBUG
			assertCorrectDomains(vars, varToVals, valToVars);
		#endif
		FlowGraph* graph = new FlowGraph(vars, varToVals, valToVars, inputVals, 
																		 lowerBounds, upperBounds, costs, 
																		 costUpperBound, true);

		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(*graph);

		//graphAlgorithms.benchmark();
		//return ES_FAILED;


		if (!graphAlgorithms.findMinCostFlow()) {
			return ES_FAILED;
		}

		vector<pair<unsigned int, unsigned int>> updatedEdges;
		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, vars, updatedEdges) != ES_OK) {
				return ES_FAILED;
		}

		(void)new (home) CostGcc(home, vars, graph, updatedEdges, ipl);
		return ES_OK;
	}

	CostGcc(Space& home, CostGcc& p) : CostGccBase(home, p) {
		c.update(home, p.c);
    x.update(home, p.x);
		graph = new FlowGraph(*(p.graph));
		updatedEdges = p.updatedEdges;
		ipl = p.ipl;
  }

	virtual Propagator *copy(Space& home) {
		return new (home) CostGcc(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::cubic(PropCost::LO, x.size());
	}

	virtual size_t dispose(Space& home) {
		home.ignore(*this, AP_DISPOSE);
		delete graph;
		updatedEdges.~vector();
    c.dispose(home);
    (void) CostGccBase::dispose(home);
    return sizeof(*this);
  }

	virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(*graph);
		if (!graphAlgorithms.updateMinCostFlow(updatedEdges)) {
			return ES_FAILED;
		}
		updatedEdges.clear();

		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, x, updatedEdges) != ES_OK) {
				return ES_FAILED;
		}
		return ES_FIX;
	}

	virtual ExecStatus advise(Space&, Advisor& a, const Delta&) {
		int xIndex = static_cast<ViewAdvisor&>(a).xIndex;
		graph->updatePrunedValues(x[xIndex], xIndex, updatedEdges);
		return ES_NOFIX;
	}

private:
	// Prune values that belong to a domain, but are not mentioned in the 
	// inputVals array/set
	ExecStatus static pruneOmittedVales(
											Space& home, ViewArray<Int::IntView>& vars, 
											MapToSet<unsigned int, int>& varToVals,
											MapToSet<int, unsigned int>& valToVars,
											const unordered_map<int, unsigned int>& inputValToIndex) {
		
		unordered_set<int> prunedVals;
		for (auto& v: valToVars.map) {
			auto value = v.first;
			if (inputValToIndex.find(value) == inputValToIndex.end()) {
				for (auto& x: v.second) {
					GECODE_ME_CHECK(vars[x].nq(home, value));
					varToVals.map.find(x)->second.erase(value);
					cout << "Prunning " << value << " from " << x << 
									" (not mentioned in values array)\n";
				}
				prunedVals.insert(value);
			}
		}

		// Also mirror the changes to valToVars
		for (auto& val: prunedVals) {
			valToVars.map.erase(val);
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
												MapToSet<unsigned int, int>& varToVals,
												MapToSet<int, unsigned int>& valToVars, 
												const IntArgs& inputVals,
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
																		 varToVals, valToVars, *prunedValsWrite,
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
																			 varToVals, valToVars, *prunedValsWrite, 
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
																		 MapToSet<unsigned int, int>& varToVals,
																	   MapToSet<int, unsigned int>& valToVars,
																	 	 unordered_set<int>& prunedVals, 
																		 const unordered_set<int>& alreadyProcessing
																		 ) {
		auto it = valToVars.map.find(value);
		unsigned int totalVarsWithThisVal = (it != valToVars.map.end() ? 
																				 it->second.size() : 0);
		if (lowerBound > totalVarsWithThisVal) {
			cout << "Lower bound of value" << value << " is greater than remaining " 
					<< "edges: " << totalVarsWithThisVal << endl;
			return ES_FAILED;
		} else if (it != valToVars.map.end() && lowerBound == totalVarsWithThisVal) {
			for (auto x: it->second) {
				cout << "Assigning var " << x << " to " << value << endl;
				assignValToVar(x, vars[x], value, varToVals, valToVars, prunedVals, 
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
	// alreadyProcessing set. If they do, we are already iterating through them.
	// In addition, update varToVals to include the assignment
	void static assignValToVar(unsigned int xIndex, Int::IntView x, int val, 
											       MapToSet<unsigned int, int>& varToVals,
														 MapToSet<int, unsigned int>& valToVars, 
														 unordered_set<int>& prunedVals, 
														 const unordered_set<int>& alreadyProcessing) {
		for (IntVarValues v(x); v(); ++v) {
			if (v.val() != val) {
				if (alreadyProcessing.find(v.val()) == alreadyProcessing.end()) {
					prunedVals.insert(v.val());
				}
				auto otherValVars = valToVars.map.find(v.val());
				assert(otherValVars != valToVars.map.end());
				otherValVars->second.erase(xIndex);
				if (otherValVars->second.empty()) {
					valToVars.map.erase(otherValVars);
				}
			}
		}
		auto varToValsEntry = varToVals.map.find(xIndex);
		assert(varToValsEntry != varToVals.map.end());
		varToValsEntry->second.clear();
		varToValsEntry->second.insert(val);
	}

	#ifndef NDEBUG
	// Assert that Gecode variable domains and valToVars/varToVals are in sync
	void static assertCorrectDomains(const ViewArray<Int::IntView>& vars, 
															 		 const MapToSet<unsigned int, int>& varToVals,
																	 const MapToSet<int, unsigned int>& valToVars
																	) {
		for (int x = 0; x < vars.size(); x++) {
			auto varToValsEntry = varToVals.map.find(x);
			assert(varToValsEntry != varToVals.map.end());
			for (IntVarValues v(vars[x]); v(); ++v) {
				assert(varToValsEntry->second.find(v.val()) 
							!= varToValsEntry->second.end());
							
				auto it = valToVars.map.find(v.val());
				assert(it != valToVars.map.end());
				assert(it->second.find(x) != it->second.end());
			} 
		}

		for (auto& x: varToVals.map) {
			assert((int) x.first < vars.size());
			assert(x.second.size() == vars[x.first].size());
			for (auto v: x.second) {
				assert(vars[x.first].in(v));
			}
		}

		for (auto& v: valToVars.map) {
			for (auto x: v.second) {
				assert((int) x < vars.size());
				assert(vars[x].in(v.first));
			}
		}
	}
	#endif
};
