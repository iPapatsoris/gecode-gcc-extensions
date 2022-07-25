#ifndef H_SYM_GCC
#define H_SYM_GCC

#include <gecode/set.hh>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>
#include "util.hpp"
#include "flow-graph-algorithms.hpp"

using namespace Gecode;
using namespace std;

typedef NaryPropagator<Set::SetView, Set::PC_SET_NONE> SymGccBase;
/**
 * The base algorithm is according to "W. Kocjan, P. Kreuger, Filtering Methods 
 * for Symmetric Cardinality Constraint, Integration of AI and OR
 * Techniques in Constraint Programming for Combinatorial Optimization 
 * Problems, First International Conference, CPAIOR 2004"
 * 
 * The following optimization has been applied not mentioned in 
 * the paper:
 * - Efficient structure to backtrack the graph by copying
 *   only an integer instead of the whole graph on each branch,
 *   taken from "P. Nightingale The extended global cardinality constraint: An 
 *   empirical survey, Artificial Intelligence 175(2), pp. 586-614, 2011"
 */
class SymGcc : public SymGccBase {

protected:
		class ViewAdvisor : public Advisor {
			public:
				Set::SetView x;
				unsigned int xIndex;
				ViewAdvisor(Space& home, Propagator& p, Council<ViewAdvisor>& c, 
										Set::SetView x, unsigned int xIndex) 
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
		FlowGraph graph;
		IntPropLevel ipl;

public:
	SymGcc(Space& home, ViewArray<Set::SetView> x, FlowGraph graph,
				 IntPropLevel ipl)
			: NaryPropagator(home, x), c(home), graph(graph), 
				ipl(ipl) {
		for (int i = 0; i < x.size(); i++) {
			(void)new (home) ViewAdvisor(home, *this, c, x[i], i);
		}
		home.notice(*this, AP_DISPOSE);
	}

	static ExecStatus post(Space& home, ViewArray<Set::SetView>& vars,
												const vector<unordered_set<int> >& varToVals,
												MapToSet& valToVars,
												const IntArgs& inputVals, 
												const IntArgs& lowerValBounds, 
												const IntArgs& upperValBounds,
												const IntArgs& lowerVarBounds, 
												const IntArgs& upperVarBounds, 
												IntPropLevel ipl) {

		#ifndef NDEBUG
			//assertCorrectDomains(vars, valToVars);
		#endif
		FlowGraph graph = FlowGraph(vars, varToVals, valToVars, inputVals, 
																		 lowerValBounds, upperValBounds, 
																		 lowerVarBounds, upperVarBounds);

		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(graph);

		if (!graphAlgorithms.findFlow(vars)) {
			return ES_FAILED;
		}

		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, vars) 
				!= ES_OK) {
				return ES_FAILED;
		}

		(void)new (home) SymGcc(home, vars, graph, ipl);
		return ES_OK;
	}

	SymGcc(Space& home, SymGcc& p) : SymGccBase(home, p), graph(p.graph) {
		c.update(home, p.c);
    x.update(home, p.x);
		ipl = p.ipl;
  }

	virtual Propagator *copy(Space& home) {
		return new (home) SymGcc(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::cubic(PropCost::HI, x.size());
	}

	virtual size_t dispose(Space& home) {
		home.ignore(*this, AP_DISPOSE);
		graph.~FlowGraph();
    c.dispose(home);
    (void) SymGccBase::dispose(home);
    return sizeof(*this);
  }

	virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(graph);
		if (!graphAlgorithms.updateFlow(x)) {
			return ES_FAILED;
		}

		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, x) 
				!= ES_OK) {
				return ES_FAILED;
		}
		return ES_FIX;
	}

	virtual ExecStatus advise(Space&, Advisor& a, const Delta&) {
		int xIndex = static_cast<ViewAdvisor&>(a).xIndex;
		bool isFeasible = graph.updatePrunedValues(x[xIndex], xIndex);
		return isFeasible ? ES_FIX : ES_NOFIX;
	}
};

#endif