#ifndef H_FLOW_GRAPH
#define H_FLOW_GRAPH

#include "node.hpp"
#include "BestBranch.hpp"
#include "util.hpp"
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
/**
 * Graph used to solve a feasible flow problem, to prove consistency of
 * symgcc, and to achieve arc consistency.
 * There is a node for each variable, for each value, and an S and T node.
 * There are value->variable edges according to the respective domains,
 * variable->T edges, S->value edges and a T->S edge.
 * The graph is built once and is backtracked in an efficient way using 
 * BtVector. As values get pruned, the corresponding edges are removed and
 * the flow is incrementally repaired to a feasible one. The flow itself is 
 * not backtracked.
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
			
			// Holds each variable's domain. Is needed to find out which values got
			// pruned between iterations, by comparing with Gecode's variable domains.
			// We know which variables got changed by using advisors, so domain
			// comparison is only done among those. We need it for fast lookup,
			// because the graph we hold has Val->Var edges and not the inverse.
			vector<BtVector<int>> varToVals;
			int totalVarNodes; 

			BacktrackStableContent() {}
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

		// Each time the advisor is executed, record the edges with bound
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

		// Search for source->dest edge, return pointer to it or NULL if it 
		// doesn't exist
		NormalEdge* getEdge(int source, int dest) {
			return backtrackStable->nodeList[source].edgeList.getVal(dest, 
																										      edgeListSize[source]);
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

	public:

		FlowGraph(
			const ViewArray<Set::SetView>& vars, 
			const vector<unordered_set<int> >& varToVals,
			const MapToSet& valToVars,
			const IntArgs& inputVals, const IntArgs& lowerValBounds, 
			const IntArgs& upperValBounds, const IntArgs& lowerVarBounds, 
			const IntArgs& upperVarBounds);

		// If a variable-value pair is pruned that has no flow, delete it on the 
		// spot. If it has flow, insert it on updatedEdges, for the flow repair
		// algorithm to fix it on the next propagation, and mark flow as infeasible.
		// If a variable is assigned but has no flow, mark flow as infeasible. 
		// Returns whether the current flow is still feasible or not.
		bool updatePrunedValues(Set::SetView x, int xIndex); 

		void print() const;
		void printResidual() const; 
		void printBounds(int x) const;
};

#endif