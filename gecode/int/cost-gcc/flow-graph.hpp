#ifndef H_FLOW_GRAPH
#define H_FLOW_GRAPH

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include "node.hpp"
#include "example/BestBranch.hpp"
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
 */
class FlowGraph {
	friend class FlowGraphAlgorithms;
	private:
		/**
		 * Holds content that does not need to be backtracked, thus should not be
		 * copied when the space is cloned.
		 */
		struct BacktrackStableContent {
			// Value network graph including edge bounds and flow values.
			// The node order in nodeList is variables,values,S,T.
			vector<Node> nodeList;

			// Fast lookup of what node a value corresponds to
			unordered_map<int, int> valToNode;
			
			// Fast lookup of what value a node corresponds to
			unordered_map<int, int> nodeToVal;
			
			// Total flow through the graph. Since the flow is not 
			// backtracked, its cost should also not be.
			int flowCost;

			// Holds each variable's domain. Is needed to find out which values got
			// pruned between iterations, by comparing with Gecode's variable domains.
			// We know which variables got changed by using advisors, so domain
			// comparison is only done among those. We need it for fast lookup,
			// because the graph we hold has Val->Var edges and not the inverse.
			vector<BtVector<int>> varToVals;

			// The total number of nodes that correspond to a variable
			int totalVarNodes; 
		
			// Mark the first time we find a solution. The first solution is of minimal
			// cost, use this variable to print it for debugging or testing
			bool firstTimeValidCost;

			BacktrackStableContent() : flowCost(0) {}
		};

		// Heap allocation to avoid deep copies, with smart pointer so that at 
		// the end of the program the memory will be deallocated. Cannot be done 
		// manually with raw pointer without a memory leak. We pack all the content
		// in one structure, to maintain only one shared_ptr, as they can be
		// expensive.
		shared_ptr<BacktrackStableContent> backtrackStable;

		// These two fields normally belong to BtVector, but we place them 
		// externally because they need to backtracked, while the rest of BtVector
		// does not. We need one "size" int field for each BtVector instance that 
		// exists in the program. On deletion of a BtVector element, we
		// swap the element with the last one, and decrement the respective 
		// size field accordingly (the actual vector size remains unchanged).
		// This allows us to backtrack to previous graph states just by recovering
		// the old value of the size, without needing to copy the graph each time. 
		vector<int> edgeListSize;
		vector<int> varToValsSize;

		// Each time the advisor is executed, record the edges with upper bound
		// violations that need to be fixed at the next execution of the propagator.
		vector<EdgeInfo> updatedEdges;
		
		// Position of S node
		int sNode() const { return backtrackStable->nodeList.size() - 2; }
		// Position of T node
		int tNode() const { return backtrackStable->nodeList.size() - 1; }

		// Search for an edge flow violating lower bounds
		// Return false if none exists
		bool getLowerBoundViolatingEdge(EdgeInfo& violation) const {
			for (unsigned int i = 0; i < backtrackStable->nodeList.size(); i++) {
				auto& edges = backtrackStable->nodeList[i].edgeList;
				for (int e = 0; e < edgeListSize[i]; e++) {
					auto& edge = (edges.list)[e];
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
			return backtrackStable->nodeList[src].edgeList.getVal(dest, 
																													  edgeListSize[src]);
		}

		void deleteEdge(int src, int dest) {
			backtrackStable->nodeList[src].edgeList.deleteVal(dest, 
																												&edgeListSize[src]);
		}
		
		void deleteResidualEdge(int src, int dest) {
			bool found = false;
			auto& residual = backtrackStable->nodeList[src].residualEdgeList;
			for (auto it = residual.begin(); it != residual.end(); it++) {
				if (it->destNode == dest) {
					residual.erase(it);
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
			for (unsigned int i = 0; 
					 i < backtrackStable->nodeList[src].residualEdgeList.size(); 
					 i++) {
				ResidualEdge& edge = backtrackStable->nodeList[src].residualEdgeList[i];
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
				backtrackStable->nodeList[src].residualEdgeList.push_back(newEdge);
			}
		}

		// Check flow cost validity against upper bound
		bool checkFlowCost(Int::IntView costUpperBound) {
			if (backtrackStable->firstTimeValidCost && backtrackStable->flowCost <= 
															  costUpperBound.max()) {
				backtrackStable->firstTimeValidCost = false;
				cout << backtrackStable->flowCost << "\n";
			}
			return backtrackStable->flowCost <= costUpperBound.max();
		}

		void calculateReducedCosts(const vector<int>& distances) {
			for (unsigned int i = 0; i < backtrackStable->nodeList.size(); i++) {
				for (auto& edge : backtrackStable->nodeList[i].residualEdgeList) {
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
		bool updatePrunedValues(Int::IntView x, int xIndex); 

		void print() const;

		void printResidual() const; 

		// Add T->Var residual edges. After the initial min cost flow is
		// established, normally the T->Var residual edges are removed
		// (because Var->T edge has flow which is equal to both its lower and upper 
		// bounds). We add them back to be able to calculate the distances 
		// from T to all other nodes during arc consistency, for calculating 
		// the reduced costs.
		void addTResidualEdges() {
			for (int var = 0; var < backtrackStable->totalVarNodes; var++) {
				backtrackStable->nodeList[tNode()].residualEdgeList.push_back(
																					 ResidualEdge(var, 1, 0));
			}
		}
};

#endif