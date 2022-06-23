#ifndef H_FLOW_GRAPH
#define H_FLOW_GRAPH

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include "node.hpp"
#include "example/LI.hpp"
#include "bt-vector.hpp"
#include "util.hpp"

using namespace Gecode;
using namespace std;

/**
 * Graph used to solve min cost flow problem, to prove consistency of
 * costgcc, and to achieve arc consistency.
 * There is a node for each variable, for each value, and an S and T node.
 * There are value->variable edges according to the respective domains,
 * variable->T edges, S->value edges and a T->S edge.
 * The graph is built once and is backtracked in an efficient way using 
 * BtVector. As values get pruned, the corresponding edges are removed and
 * the flow is incrementally repaired to a feasible one. The flow itself is 
 * not backtracked. Thus when the a graph is backtracked, we need to make sure
 * that we are still on a min cost flow, and if not repair it, since the 
 * addition of previously deleted edges might remove the optimality of our flow.
 * A pointer is used for members that don't need to be backtracked, to have them
 * reside in the same memory location and not copy on each iteration.
 */
class FlowGraph {
	friend class FlowGraphAlgorithms;
	private:

		// The node order in nodeList is variables,values,S,T.
		vector<Node> nodeList;
		
		// The total number of nodes that correspond to a variable
		int totalVarNodes; 
		
		// Fast lookup of what node a value corresponds to
		unordered_map<int, int> *valToNode;
		
		// Fast lookup of what value a node corresponds to
		unordered_map<int, int> *nodeToVal;
		
		// Holds each variable's domain. Is needed to find out which values got
		// pruned between iterations, by comparing with Gecode's variable domains.
		// We know which variables got changed by using advisors, so domain
		// comparison is only done among those. We need it for fast lookup,
		// because the graph we hold has Val->Var edges and not the other direction.
		vector<BtVector<int>> varToVals;

		// Total flow through the graph, starts at 0. Since the flow is not 
		// backtracked, its cost should also not be, thus the use of pointer.
		int *flowCost;

		// Position of S node
		int sNode() const { return nodeList.size() - 2; }
		// Position of T node
		int tNode() const { return nodeList.size() - 1; }

		// Mark the first time we find a solution. The first solution is of minimal
		// cost, use this variable to print it for debugging or testing
		bool firstTimeValidCost;

		// Search for an edge flow violating lower bounds
		// Return false if none exists
		bool getLowerBoundViolatingEdge(EdgeInfo& violation) const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				auto& edges = nodeList[i].edgeList;
				for (int e = 0; e < edges.listSize; e++) {
					auto& edge = (*edges.list)[e];
					if (edge.flow < edge.lowerBound) {
						violation.src = i;
						violation.dest = edge.destNode;
						return true;
					}
				}
			}
			return false;
		}

		// Search for src->dest edge, return pointer to it or NULL if it 
		// doesn't exist
		NormalEdge* getEdge(int src, int dest) {
			return nodeList[src].edgeList.getVal(dest);
		}

		void deleteEdge(int src, int dest) {
			nodeList[src].edgeList.deleteVal(dest);
		}
		
		void deleteResidualEdge(int src, int dest) {
			bool found = false;
			auto residual = nodeList[src].residualEdgeList;
			for (auto it = residual->begin(); it != residual->end(); it++) {
				if (it->destNode == dest) {
					residual->erase(it);
					found = true;
					break;
				}
			}
			assert(found);
			if (!found) {
				cout << "Internal error: did not find residual edge" << endl;
				exit(1);
			}
		}

		// Search for src->dest residual edge, return pointer to it or NULL if it
		// doesn't exists. If index is not NULL, also return its position in that
		// node's residual edges list 
		ResidualEdge* getResidualEdge(int src, int dest, 
																  int *index = NULL) {
			for (unsigned int i=0; i < nodeList[src].residualEdgeList->size(); i++) {
				ResidualEdge& edge = (*nodeList[src].residualEdgeList)[i];
				if (edge.destNode == dest) {
					if (index != NULL) {
						*index = i;
					}
					return &edge;
				}
			}
			return NULL;
		}

		// If existingEdge is NULL, create edge from src to newEdge.dest
		// If existingEdge is not NULL, change it to newEdge
		void setOrCreateResidualEdge(ResidualEdge* existingEdge, int src, 
																 const ResidualEdge& newEdge) {
			if (existingEdge != NULL) {
				*existingEdge = newEdge;
			} else {
				nodeList[src].residualEdgeList->push_back(newEdge);
			}
		}

		// Check flow cost validity against upper bound
		bool checkFlowCost(Int::IntView costUpperBound) {
			if (firstTimeValidCost && *flowCost <= costUpperBound.max()) {
				firstTimeValidCost = false;
				cout << *flowCost << "\n";
			}
			return *flowCost <= costUpperBound.max();
		}

		void calculateReducedCosts(const vector<int>& distances) {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				for (auto& edge : (*nodeList[i].residualEdgeList)) {
					edge.reducedCost = distances[i] + edge.cost - 
														 distances[edge.destNode];
				}
			}
		}

	public:

		FlowGraph(
			const ViewArray<Int::IntView>& vars, 
 			const vector<unordered_set<int> >& varToVals,
			const MapToSet& valToVars,
			const IntArgs& inputVals, const IntArgs& lowerBounds, 
			const IntArgs& upperBounds, const IntArgs& costs);

		// Update graph state to match variable X domain prunings.
		// If a variable-value pair is pruned that has no flow, delete it on the 
		// spot. If it has flow, insert it on updatedEdges, for the flow repair
		// algorithm to fix it on the next propagation, and mark flow as infeasible.
		// If a variable is assigned but has no flow, mark flow as infeasible. 
		// Returns whether the current flow is still feasible or not.
		bool updatePrunedValues(Int::IntView x, int xIndex, 
													  vector<EdgeInfo>& updatedEdges); 

		void print() const;

		void printResidual() const; 

		// Add T->Var residual edges. After the initial min cost flow is
		// established, normally the T->Var residual edges are removed
		// (because Var->T edge has flow which is equal to both its lower and upper 
		// bounds). We add them back to be able to calculate the distances 
		// from T to all other nodes during arc consistency, for calculating 
		// the reduced costs.
		void addTResidualEdges() {
			for (int var = 0; var < totalVarNodes; var++) {
				nodeList[tNode()].residualEdgeList->push_back(ResidualEdge(var, 1, 0));
			}
		}
};

#endif