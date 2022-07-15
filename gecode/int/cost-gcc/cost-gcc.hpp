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

/**
 * The base algorithm is according to "J-C. Régin, Cost-Based Arc Consistency 
 * for Global Cardinality Constraints, Constraints 7, 2002, pp. 387-405"
 * The following changes and optimizations have been applied not mentioned in 
 * the paper:
 * - Shortest Path Faster Algorithm instead of Dijkstra.
 * - Efficient structure to backtrack the graph by copying
 *   only an integer instead of the whole graph on each branch,
 *   taken from "P. Nightingale The extended global cardinality constraint: An 
 *   empirical survey, Artificial Intelligence 175(2), pp. 586-614, 2011"
 * - Custom branching heuristic utilizing the current min cost flow
 *   to guide the search straight to the next solution, minimizing 
 *   failures.
 */
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
		FlowGraph graph;
		BestBranch bestBranch;
		bool usingBestBranch;
		IntPropLevel ipl;

public:
	CostGcc(Space& home, ViewArray<Int::IntView> x, FlowGraph graph, 
					const vector<EdgeInfo>& updatedEdges, BestBranch* bestBranch, IntPropLevel ipl, 
					Int::IntView costUpperBound)
			: NaryPropagator(home, x), c(home), costUpperBound(costUpperBound), 
				graph(graph), usingBestBranch(bestBranch != NULL), ipl(ipl) {
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
		FlowGraph graph = FlowGraph(vars, varToVals, valToVars, inputVals, 
											lowerBounds, upperBounds, costs);

		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(graph);
		if (!graphAlgorithms.findMinCostFlow(bestBranch, costUpperBound)) {
			return ES_FAILED;
		}
		graph.addTResidualEdges();
		vector<EdgeInfo> updatedEdges;
		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, vars, 
																										 costUpperBound) != ES_OK) {
				return ES_FAILED;
		}

		(void)new (home) CostGcc(home, vars, graph, updatedEdges, bestBranch, ipl, 
														 costUpperBound);
		return ES_OK;
	}

	CostGcc(Space& home, CostGcc& p) : CostGccBase(home, p), graph(p.graph) {
		c.update(home, p.c);
    x.update(home, p.x);
		costUpperBound.update(home, p.costUpperBound);
		usingBestBranch = p.usingBestBranch;
		if (usingBestBranch) {
			bestBranch.update(home, p.bestBranch);
		}
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
		graph.~FlowGraph();
    c.dispose(home);
		// todo: delete cost?
    (void) CostGccBase::dispose(home);
    return sizeof(*this);
  }

	virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(graph);
		if (!graphAlgorithms.updateMinCostFlow(usingBestBranch ? &bestBranch : NULL,
																			     costUpperBound)) {
			return ES_FAILED;
		}

		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, x, 
																										 costUpperBound) != ES_OK) {
				return ES_FAILED;
		}
		return ES_FIX;
	}

	virtual ExecStatus advise(Space&, Advisor& a, const Delta&) {
		int xIndex = static_cast<ViewAdvisor&>(a).xIndex;
		bool isFeasible = graph.updatePrunedValues(x[xIndex], xIndex);
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