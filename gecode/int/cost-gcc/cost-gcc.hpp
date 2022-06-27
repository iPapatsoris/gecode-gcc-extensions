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
				int xIndex;
				ViewAdvisor(Space& home, Propagator& p, Council<ViewAdvisor>& c, 
										Int::IntView x, int xIndex) 
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
		Int::IntView costUpperBound;
		FlowGraph* graph;
		vector<EdgeInfo> updatedEdges;
		BestBranch bestBranch;
		bool usingBestBranch;
		IntPropLevel ipl;

public:
	CostGcc(Space& home, ViewArray<Int::IntView> x, FlowGraph* graph, 
					const vector<EdgeInfo>& updatedEdges, BestBranch* bestBranch, IntPropLevel ipl, 
					Int::IntView costUpperBound)
			: NaryPropagator(home, x), c(home), costUpperBound(costUpperBound), 
				graph(graph), updatedEdges(updatedEdges), 
				usingBestBranch(bestBranch != NULL), ipl(ipl) {
		for (int i = 0; i < x.size(); i++) {
			(void)new (home) ViewAdvisor(home, *this, c, x[i], i);
		}
		if (usingBestBranch) {
			this->bestBranch = *bestBranch;
		}
		home.notice(*this, AP_DISPOSE);
	}

	static ExecStatus post(Space& home, ViewArray<Int::IntView>& vars,
												vector<unordered_set<int> >& varToVals,
												MapToSet& valToVars,
												const IntArgs& inputVals, 
												const IntArgs& lowerBounds, const IntArgs& upperBounds,
												const IntArgs& costs, Int::IntView costUpperBound,
												BestBranch* bestBranch, IntPropLevel ipl) {

		#ifndef NDEBUG
			assertCorrectDomains(vars, varToVals, valToVars);
		#endif
		FlowGraph* graph = new FlowGraph(vars, varToVals, valToVars, inputVals, 
																		 lowerBounds, upperBounds, costs);

		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(*graph);
		if (!graphAlgorithms.findMinCostFlow(bestBranch, costUpperBound)) {
			return ES_FAILED;
		}
		graph->addTResidualEdges();
		vector<EdgeInfo> updatedEdges;
		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, vars, 
																										 costUpperBound) != ES_OK) {
				return ES_FAILED;
		}

		(void)new (home) CostGcc(home, vars, graph, updatedEdges, bestBranch, ipl, 
														 costUpperBound);
		return ES_OK;
	}

	CostGcc(Space& home, CostGcc& p) : CostGccBase(home, p) {
		c.update(home, p.c);
    x.update(home, p.x);
		costUpperBound.update(home, p.costUpperBound);
		usingBestBranch = p.usingBestBranch;
		if (usingBestBranch) {
			bestBranch.update(home, p.bestBranch);
		}
		graph = new FlowGraph(*(p.graph));
		updatedEdges = p.updatedEdges;
		ipl = p.ipl;
  }

	virtual Propagator *copy(Space& home) {
		return new (home) CostGcc(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::cubic(PropCost::HI, x.size());
	}

	virtual size_t dispose(Space& home) {
		home.ignore(*this, AP_DISPOSE);
		delete graph;
		updatedEdges.~vector();
    c.dispose(home);
		// todo: delete cost?
    (void) CostGccBase::dispose(home);
    return sizeof(*this);
  }
	virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(*graph);
		if (!graphAlgorithms.updateMinCostFlow(updatedEdges, 
															 						 usingBestBranch ? &bestBranch : NULL, 
																					 costUpperBound)) {
			return ES_FAILED;
		}
		updatedEdges.clear();

		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, x, 
																										 costUpperBound) != ES_OK) {
				return ES_FAILED;
		}
		return ES_FIX;
	}

	virtual ExecStatus advise(Space&, Advisor& a, const Delta&) {
		int xIndex = static_cast<ViewAdvisor&>(a).xIndex;
		bool isFeasible = graph->updatePrunedValues(x[xIndex], xIndex, 
																							  updatedEdges);
		return isFeasible ? ES_FIX : ES_NOFIX;
	}

private:
	#ifndef NDEBUG
	// Assert that Gecode variable domains and valToVars/varToVals are in sync
	void static assertCorrectDomains(const ViewArray<Int::IntView>& vars, 
															 		 vector<unordered_set<int> >& varToVals,
																	 const MapToSet& valToVars
																	) {
		assert(varToVals.size() == vars.size());
		for (int x = 0; x < vars.size(); x++) {
			auto varToValsEntry = varToVals[x];
			for (IntVarValues v(vars[x]); v(); ++v) {
				assert(varToValsEntry.find(v.val()) 
							!= varToValsEntry.end());
							
				auto it = valToVars.find(v.val());
				assert(it != valToVars.end());
				assert(it->second.find(x) != it->second.end());
			} 
		}

		for (int x = 0; x < varToVals.size(); x++) {
			assert(varToVals[x].size() == vars[x].size());
			for (auto v: varToVals[x]) {
				assert(vars[x].in(v));
			}
		}

		for (auto& v: valToVars) {
			for (auto x: v.second) {
				assert((int) x < vars.size());
				assert(vars[x].in(v.first));
			}
		}
	}
	#endif
};

#endif