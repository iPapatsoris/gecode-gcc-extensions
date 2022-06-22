#ifndef H_FLOW_GRAPH
#define H_FLOW_GRAPH

#include "node.hpp"
#include "example/LI.hpp"
#include "bt-vector.hpp"
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
struct EdgeUpdate {
	int src;
	int dest;

	EdgeUpdate(int src, int dest) : src(src), dest(dest)
						{}
};

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
		//vector<int> debug;
		int totalVarNodes; 
		// Fast lookup of what node a value corresponds to
		unordered_map<int, int> *valToNode; 			  // TODO: memory leak, use shared object
		// Fast lookup of what value a node corresponds to
		unordered_map<int, int> *nodeToVal;				 // TODO: same
		// Holds each variable's domain. Is needed to find out which values got
		// pruned during an iteration, by comparing with Gecode's variable domains.
		// We know which variables got changed by using advisors, so domain
		// comparison is only done among those.
		// Should be kept up to date with assignments and pruning.
		vector<BtVector<int>> varToVals;

	//	vector<int> *dist;

		// Total flow through the graph, starts at 0. Is calculated at once using
		// appropriate function, not gradually
		int *flowCost;

		// Cost upper bound as defined by the constraint input
		// int costUpperBound;

		// When values are pruned or variables are assigned, we update the 
		// bounds of the corresponding edges. At that point, set this variable
		// to note whether our old flow is still feasible, or if we need to find a 
		// new one. It is much more efficient to check this when we update 
		// the bounds of specific edges, than to scan the whole graph later to see 
		// if the old flow still stands
		bool *oldFlowIsFeasible;

		// Position of S node
		int sNode() const { return nodeList.size() - 2; }
		// Position of T node
		int tNode() const { return nodeList.size() - 1; }

		// Mark the first time we find a solution. The first solution is of minimal
		// cost, use this variable to print it for debugging or testing
		bool firstTimeValidCost;

		// Search for an edge flow violating lower bounds
		// Return false if none exists
		// TODO: can be optimized to look for less?
		bool getLowerBoundViolatingEdge(pair<int, int>& violation) 
			const {
			for (int i = 0; i < nodeList.size(); i++) {
				/*if (i == totalVarNodes) {
					i = sNode();				// if at init some values are already pruned,
															// we might have tightened var->val bounds,
															// so we need to check through every node
				}*/
				auto& edges = nodeList[i].edgeList;
				for (int e = 0; e < edges.listSize; e++) {
					auto& edge = (*edges.list)[e];
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
		NormalEdge* getEdge(int source, int dest) {
			return nodeList[source].edgeList.getVal(dest);
		}

		void deleteEdge(int source, int dest) {
			nodeList[source].edgeList.deleteVal(dest);
		}
		
		void deleteResidualEdge(int source, int dest) {
			bool found = false;
			auto residual = nodeList[source].residualEdgeList;
			for (auto it = residual->begin(); it != residual->end(); it++) {
				if (it->destNode == dest) {
					residual->erase(it);
					found = true;
					break;
				}
			}

			residual = nodeList[dest].residualEdgeList;
			for (auto it = residual->begin(); it != residual->end(); it++) {
				if (it->destNode == source) {
					residual->erase(it);
					found = true;
					break;
				}
			}
			if (!found)
				cout << "RESIDUAL EDGE NOT FOUND" << endl;
		}

		// Search for source->dest residual edge, return pointer to it or NULL if it
		// doesn't exists. If index is not NULL, also return its position in that
		// node's residual edges list 
		ResidualEdge* getResidualEdge(int source, int dest, 
																  int *index = NULL) {
			for (int i=0; i < nodeList[source].residualEdgeList->size(); i++) {
				ResidualEdge& edge = (*nodeList[source].residualEdgeList)[i];
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
															   int source, 
																 const ResidualEdge& newEdge) {
			if (existingEdge != NULL) {
				*existingEdge = newEdge;
			} else {
				nodeList[source].residualEdgeList->push_back(newEdge);
			}
		}

		// Iterate through each edge that has flow, to find its total cost
		int calculateFlowCost(LI &lii);

		bool checkFlowCost(Int::IntView costUpperBound) {
			if (firstTimeValidCost && *flowCost <= costUpperBound.max()) {
				firstTimeValidCost = false;
			}
			// cout << *flowCost << " " << costUpperBound.max() << endl;
			return *flowCost <= costUpperBound.max();
		}

		void calculateReducedCosts(const vector<int>& distances) {
			for (int i = 0; i < nodeList.size(); i++) {
				for (auto& edge : (*nodeList[i].residualEdgeList)) {
					edge.reducedCost = distances[i] + edge.cost - distances[edge.destNode];
					// cout << i << "->" << edge.destNode << " " << edge.reducedCost << " = " << distances[i] << " + " << edge.cost << " - " << distances[edge.destNode] << endl;
				}
			}
		}

		#ifndef NDEBUG
		// Assert varToVals is synchronized with Gecode variable X domain
		/*void assertVarToValsInSync(Int::IntView x, int xIndex) const {
			auto vals = (*varToVals)[xIndex];
			assert(vals.size() == x.size());
			for (IntVarValues v(x); v(); ++v) {
				assert(vals.find(v.val()) != vals.end());
			}
			for (auto val: vals) {
				assert(x.in(val));
			}
		}*/
		#endif

	public:

		FlowGraph(
			const ViewArray<Int::IntView>& vars, 
 			const vector<unordered_set<int> >& varToVals,
			const MapToSet<int, int>& valToVars,
			const IntArgs& inputVals, const IntArgs& lowerBounds, 
			const IntArgs& upperBounds, const IntArgs& costs);

		// Update graph state to match variable X domain pruning/assignment.
		// Update is made by tightening the bounds of edge V->X as follows:
		// - If X got assigned to value V, set the lower bound to 1.
		// - For every value V that has been pruned off X, set the upper bound 
		//   to 0. 
	  // If we prune a value that is used by current flow, or assign a value 
		// that is not used by it, set oldFlowIsFeasible to false.
		// Populate updatedEdges, so we know where we should update the old residual
		// graph later on
		bool updatePrunedValues(Int::IntView x, int xIndex, 
													  vector<EdgeUpdate>& updatedEdges); 

		void print() const;

		void printResidual() const; 


		bool getOldFlowIsFeasible() const {
			return *oldFlowIsFeasible;
		}

		void addTResidualEdges() {
			for (int var = 0; var < totalVarNodes; var++) {
				nodeList[tNode()].residualEdgeList->push_back(ResidualEdge(var, 1, 0));
			}
		}

		void removeTResidualEdges() {
			nodeList[tNode()].residualEdgeList->clear();
		}

		void debug() const {
			if (nodeList[0].residualEdgeList == nodeList[15].residualEdgeList) {
				cout << "SAME!" << endl;
				exit(1);
			}
		}
};

#endif