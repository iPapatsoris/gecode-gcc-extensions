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

typedef NaryPropagator<Set::SetView, Set::PC_SET_ANY> SymGccBase;

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
		FlowGraph* graph;
		vector<EdgeNodes> updatedEdges;
		LI li;
		bool usingLocalHandle;
		// TODO: do not store, instead use different post functions?
		IntPropLevel ipl;

public:
	SymGcc(Space& home, ViewArray<Set::SetView> x, FlowGraph* graph, 
					const vector<EdgeNodes>& updatedEdges, LI* li, 
					IntPropLevel ipl)
			: NaryPropagator(home, x), c(home), graph(graph), 
				updatedEdges(updatedEdges), usingLocalHandle(li != NULL), 
				ipl(ipl) {
		for (int i = 0; i < x.size(); i++) {
			(void)new (home) ViewAdvisor(home, *this, c, x[i], i);
		}
		if (usingLocalHandle) {
			this->li = *li;
		}
		home.notice(*this, AP_DISPOSE);
	}

	static ExecStatus post(Space& home, ViewArray<Set::SetView>& vars,
												MapToSet<int, unsigned int>& valToVars,
												const IntArgs& inputVals, 
												const IntArgs& lowerValBounds, 
												const IntArgs& upperValBounds,
												const IntArgs& lowerVarBounds, 
												const IntArgs& upperVarBounds, LI* li,
												IntPropLevel ipl) {

		#ifndef NDEBUG
			//assertCorrectDomains(vars, valToVars);
		#endif
		FlowGraph* graph = new FlowGraph(vars, valToVars, inputVals, 
																		 lowerValBounds, upperValBounds, 
																		 lowerVarBounds, upperVarBounds);

		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(*graph);

		if (!graphAlgorithms.findMinCostFlow(li)) {
			return ES_FAILED;
		}

		vector<pair<unsigned int, unsigned int>> updatedEdges;
		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, vars, updatedEdges) != ES_OK) {
				return ES_FAILED;
		}

		(void)new (home) SymGcc(home, vars, graph, updatedEdges, li, ipl);
		return ES_OK;
	}

	SymGcc(Space& home, SymGcc& p) : SymGccBase(home, p) {
		c.update(home, p.c);
    x.update(home, p.x);
		usingLocalHandle = p.usingLocalHandle;
		if (usingLocalHandle) {
			li.update(home, p.li);
			//cout << "SymGcc copy: " << p.li[0] << endl;
		}
		graph = new FlowGraph(*(p.graph));
		updatedEdges = p.updatedEdges;
		ipl = p.ipl;
  }

	virtual Propagator *copy(Space& home) {
		return new (home) SymGcc(home, *this);
	}

	virtual PropCost cost(const Space&, const ModEventDelta&) const {
		return PropCost::cubic(PropCost::LO, x.size());
	}

	virtual size_t dispose(Space& home) {
		home.ignore(*this, AP_DISPOSE);
		delete graph;
		updatedEdges.~vector();
    c.dispose(home);
    (void) SymGccBase::dispose(home);
    return sizeof(*this);
  }

	virtual ExecStatus propagate(Space& home, const ModEventDelta&) {
		/*
		GECODE_ME_CHECK(x[4].cardMin(home, 2));
		GECODE_ME_CHECK(x[4].include(home, 3));
	  GECODE_ME_CHECK(x[4].exclude(home, 6));


			cout << "unkown values\n";
			for (SetVarUnknownValues i(x[4]); i(); ++i)
				std::cout << i.val() << " ";
			cout << "\nlowest upper bound set\n";
			for (SetVarLubValues i(x[4]); i(); ++i)
				std::cout << i.val() << " ";
			cout << "\ngreater lower bound set\n";
			for (SetVarGlbValues i(x[4]); i(); ++i)
				std::cout << i.val() << " ";
			cout << "\ncard " << x[4].cardMin() << " " << x[4].cardMax() << endl;

		//GECODE_ME_CHECK(x[4].cardMin(home, 2));
	 // GECODE_ME_CHECK(x[4].exclude(home, 4));
	 	  


		cout << " \n\nchange " << endl;
		cout << "\nunkown values\n";
			for (SetVarUnknownValues i(x[4]); i(); ++i)
				std::cout << i.val() << " ";
			cout << "\nlowest upper bound set\n";
			for (SetVarLubValues i(x[4]); i(); ++i)
				std::cout << i.val() << " ";
			cout << "\ngreater lower bound set\n";
			for (SetVarGlbValues i(x[4]); i(); ++i)
				std::cout << i.val() << " ";
							cout << "\ncard " << x[4].cardMin() << " " << x[4].cardMax() << endl;

			cout << "\n\n";

			exit(1);
		*/
		//cout << "propagate" << endl;
		FlowGraphAlgorithms graphAlgorithms = FlowGraphAlgorithms(*graph);
		if (!graphAlgorithms.updateMinCostFlow(updatedEdges, 
																					 usingLocalHandle ? &li : NULL
			 )) {
			//cout << "yo updateMinCostFlow fail lmao" << endl;
			return ES_FAILED;
		}
		updatedEdges.clear();

		if (ipl == IPL_DOM && graphAlgorithms.performArcConsistency(home, x, updatedEdges) != ES_OK) {
				return ES_FAILED;
		}

		//graph->print();
		return ES_FIX;
	}

	virtual ExecStatus advise(Space&, Advisor& a, const Delta&) {
		int xIndex = static_cast<ViewAdvisor&>(a).xIndex;
		//cout << "\nadvisor on " << xIndex << endl;
		graph->updatePrunedValues(x[xIndex], xIndex, updatedEdges, 
															usingLocalHandle ? &li : NULL);
		/*for (auto e: updatedEdges) {
			cout << e.first << "->" << e.second << endl;
		}*/
		return graph->getOldFlowIsFeasible() ? ES_FIX : ES_NOFIX;
	}

private:
	/*#ifndef NDEBUG
	// Assert that Gecode variable domains and valToVars are in sync
	void static assertCorrectDomains(const ViewArray<Set::SetView>& vars, 
																	 const MapToSet<int, unsigned int>& valToVars
																	) {
		for (int x = 0; x < vars.size(); x++) {
				auto it = valToVars.map.find(v.val());
				assert(it != valToVars.map.end());
				assert(it->second.find(x) != it->second.end());
			} 
		}

		for (auto& v: valToVars.map) {
			for (auto x: v.second) {
				assert((int) x < vars.size());
				assert(vars[x].in(v.first));
			}
		}
	}
	#endif*/
};

#endif