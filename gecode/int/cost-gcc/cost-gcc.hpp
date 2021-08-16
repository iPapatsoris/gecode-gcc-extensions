#ifndef H_COST_GCC
#define H_COST_GCC

#include <gecode/int.hh>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include "flow-graph-algorithms.hpp"

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
		LI li;
		// TODO: do not store, instead use different post functions?
		IntPropLevel ipl;

public:
	CostGcc(Space& home, ViewArray<Int::IntView> x, FlowGraph* graph, 
					const vector<EdgeNodes>& updatedEdges, LI& li, IntPropLevel ipl)
			: NaryPropagator(home, x), c(home), graph(graph), 
				updatedEdges(updatedEdges), ipl(ipl) {
		for (int i = 0; i < x.size(); i++) {
			(void)new (home) ViewAdvisor(home, *this, c, x[i], i);
		}
		this->li = li;
		home.notice(*this, AP_DISPOSE);
	}

	static ExecStatus post(Space& home, ViewArray<Int::IntView>& vars,
												MapToSet<unsigned int, int>& varToVals,
												MapToSet<int, unsigned int>& valToVars,
												const IntArgs& inputVals, 
												const unordered_map<int, unsigned int>& inputValToIndex,
												const IntArgs& lowerBounds, const IntArgs& upperBounds,
												const IntArgs& costs, int costUpperBound, LI& li,
												IntPropLevel ipl) {

		if (pruneOmittedVales(home, vars, varToVals, valToVars, 
													inputValToIndex) == ES_FAILED) {
			return ES_FAILED;
		}

		#ifndef NDEBUG
			assertCorrectDomains(vars, varToVals, valToVars);
		#endif
		FlowGraph* graph = new FlowGraph(vars, varToVals, valToVars, inputVals, 
																		 lowerBounds, upperBounds, costs, 
																		 costUpperBound);

		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(*graph);

		if (!graphAlgorithms.findMinCostFlow(li)) {
			return ES_FAILED;
		}

		vector<pair<unsigned int, unsigned int>> updatedEdges;
		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, vars, updatedEdges) != ES_OK) {
				return ES_FAILED;
		}

		(void)new (home) CostGcc(home, vars, graph, updatedEdges, li, ipl);
		return ES_OK;
	}

	CostGcc(Space& home, CostGcc& p) : CostGccBase(home, p) {
		c.update(home, p.c);
    x.update(home, p.x);
		li.update(home, p.li);
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
		if (!graphAlgorithms.updateMinCostFlow(updatedEdges, li)) {
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
		return graph->getOldFlowIsFeasible() ? ES_FIX : ES_NOFIX;
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

#endif