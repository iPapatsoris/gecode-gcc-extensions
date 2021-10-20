#ifndef H_FLOW_GRAPH
#define H_FLOW_GRAPH

#include "graph-base-components.hpp"
#include "LI.hpp"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

using namespace Gecode;
using namespace std;

template <class T1, class T2>
struct MapToSet {
	unordered_map<T1, unordered_set<T2>> map;
	MapToSet() {}
	MapToSet(const MapToSet &c) {
		map = c.map;
	}
};

// Edge containing both source node and destination info
// Normally source is not included in Edge class, because node ID N corresponds
// to the N-th position in the node list arrays
// It is typically used in updatedEdges vectors, which hold which edges got 
// updated due to some pruning or assignment, and need to be checked on when
// we next update the residual graph. It is important to have NormalEdge object
// and not NormalEdge*, because when search finds a solution or fails, it will
// clone the graph and destroy the original, thus invalidating the pointer.
typedef pair<unsigned int, unsigned int> EdgeNodes;

/**
 * Graph used to solve min cost flow problem, to prove consistency of
 * costgcc, and also achieve arc consistency later on.
 * There is a node for each variable, for each value, and an S and T node.
 * There are value->variable edges according to the respective domains,
 * variable->T edges, S->value edges and a T->S edge.
 * The lower/upper bounds on the edges are chosen in a way to respect 
 * the original contraints, while also making sure that each variable gets
 * eventually assigned to a value.
 * 
 * The graph is built once, and later on any changes in domain such as value 
 * prunes or variable assignments, are expressed by setting value->variable 
 * upper bounds to 0 (prune) or lower bounds to 1 (assignment)
 * 
 * The graph's residual graph exists through Node's residualEdgeList member.
 * It is calculated once on construction, and later on when we have a new flow,
 * it is only modified in the edges that differ from the previous flow.
 */
class FlowGraph {
	friend class FlowGraphAlgorithms;
	private:

		// The node order in nodeList is variables,values,S,T.
		vector<Node> nodeList;
		unsigned int totalVarNodes; 
		// Fast lookup of what node a value corresponds to
		unordered_map<int, unsigned int> *valToNode; 			  // TODO: memory leak, use shared object
		// Fast lookup of what value a node corresponds to
		unordered_map<unsigned int, int> *nodeToVal;				 // TODO: same
		// Holds each variable's domain. Is needed to find out which values got
		// pruned during an iteration, by comparing with Gecode's variable domains.
		// We know which variables got changed by using advisors, so domain
		// comparison is only done among those.
		// Should be kept up to date with assignments and pruning.
		MapToSet<unsigned int, int> varToVals;

		// Total flow through the graph, starts at 0. Is calculated at once using
		// appropriate function, not gradually
		int flowCost;

		// Cost upper bound as defined by the constraint input
		int costUpperBound;

		// When values are pruned or variables are assigned, we update the 
		// bounds of the corresponding edges. At that point, set this variable
		// to note whether our old flow is still feasible, or if we need to find a 
		// new one. It is much more efficient to check this when we update 
		// the bounds of specific edges, than to scan the whole graph later to see 
		// if the old flow still stands
		bool oldFlowIsFeasible;

		// Position of S node
		unsigned int sNode() const { return nodeList.size() - 2; }
		// Position of T node
		unsigned int tNode() const { return nodeList.size() - 1; }

		// Mark the first time we find a solution. The first solution is of minimal
		// cost, use this variable to print it for debugging or testing
		bool firstTimeValidCost;

		// Search for an edge flow violating lower bounds
		// Return false if none exists
		bool getLowerBoundViolatingEdge(pair<unsigned int, unsigned int>& violation) 
			const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				for (auto& edge: nodeList[i].edgeList) {
					if (edge.flow < edge.lowerBound) {
						violation = {i, edge.destNode};
						return true;
					}
				}
			}
			return false;
		}

		// Search for source->dest edge, return pointer to it or NULL if it 
		// doesn't exist
		NormalEdge* getEdge(unsigned int source, unsigned int dest) {
			for (auto& edge: nodeList[source].edgeList) {
				if (edge.destNode == dest) {
					return &edge;
				}
			}
			return NULL;
		}
		
		// Search for source->dest residual edge, return pointer to it or NULL if it
		// doesn't exists. If index is not NULL, also return its position in that
		// node's residual edges list 
		ResidualEdge* getResidualEdge(unsigned int source, unsigned int dest, 
																  unsigned int *index = NULL) {
			for (unsigned int i=0; i<nodeList[source].residualEdgeList.size(); i++) {
				ResidualEdge& edge = nodeList[source].residualEdgeList[i];
				if (edge.destNode == dest) {
					if (index != NULL) {
						*index = i;
					}
					return &edge;
				}
			}
			return NULL;
		}

		// If existingEdge is NULL, create edge from source to newEdge.dest
		// If existingEdge is not NULL, change it to newEdge
		void setOrCreateResidualEdge(ResidualEdge* existingEdge, 
															   unsigned int source, 
																 const ResidualEdge& newEdge) {
			if (existingEdge != NULL) {
				*existingEdge = newEdge;
			} else {
				nodeList[source].residualEdgeList.push_back(newEdge);
			}
		}

		// Iterate through each edge that has flow, to find its total cost
		int calculateFlowCost(LI &lii);

		bool checkFlowCost() {
			if (firstTimeValidCost && flowCost <= costUpperBound) {
				firstTimeValidCost = false;
				cout << flowCost << endl;
			}
			return flowCost <= costUpperBound;
		}

		void calculateReducedCosts(const vector<int>& distances) {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				for (auto& edge : nodeList[i].residualEdgeList) {
					edge.reducedCost = distances[i] + edge.cost - distances[edge.destNode];
					//cout << i << "->" << edge.destNode << " " << edge.reducedCost << " = " << distances[i] << " + " << edge.cost << " - " << distances[edge.destNode] << endl;
				}
			}
		}

		#ifndef NDEBUG
		// Assert varToVals is synchronized with Gecode variable X domain
		void assertVarToValsInSync(Int::IntView x, int xIndex) const {
			auto vals = varToVals.map.find(xIndex)->second;
			assert(vals.size() == x.size());
			for (IntVarValues v(x); v(); ++v) {
				assert(vals.find(v.val()) != vals.end());
			}
			for (auto val: vals) {
				assert(x.in(val));
			}
		}
		#endif

	public:

		FlowGraph(
			const ViewArray<Int::IntView>& vars, 
 			const MapToSet<unsigned int, int>& varToVals,
			const MapToSet<int, unsigned int>& valToVars,
			const IntArgs& inputVals, const IntArgs& lowerBounds, 
			const IntArgs& upperBounds, const IntArgs& costs, int costUpperBound);

		// Update graph state to match variable X domain pruning/assignment.
		// Update is made by tightening the bounds of edge V->X as follows:
		// - If X got assigned to value V, set the lower bound to 1.
		// - For every value V that has been pruned off X, set the upper bound 
		//   to 0. 
	  // If we prune a value that is used by current flow, or assign a value 
		// that is not used by it, set oldFlowIsFeasible to false.
		// Populate updatedEdges, so we know where we should update the old residual
		// graph later on
		void updatePrunedValues(Int::IntView x, unsigned int xIndex, 
													  vector<EdgeNodes>& updatedEdges); 

		void print() const;

		void printResidual() const; 


		bool getOldFlowIsFeasible() const {
			return oldFlowIsFeasible;
		}

		void addTResidualEdges() {
			for (unsigned int var = 0; var < totalVarNodes; var++) {
				nodeList[tNode()].residualEdgeList.push_back(ResidualEdge(var, 1, 0));
			}
		}

		void removeTResidualEdges() {
			nodeList[tNode()].residualEdgeList.clear();
		}
};

#endif