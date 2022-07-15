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

		int totalVarNodes; 
		bool debug; 

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

		// TODO: update this comment

		// Update graph state to match variable X domain pruning/assignment.
		// Update is made by tightening the bounds of edge V->X as follows:
		// - If X got assigned to value V, set the lower bound to 1.
		// - For every value V that has been pruned off X, set the upper bound 
		//   to 0. 
	  // If we prune a value that is used by current flow, or assign a value 
		// that is not used by it, set oldFlowIsFeasible to false.
		// Populate updatedEdges, so we know where we should update the old residual
		// graph later on
		bool updatePrunedValues(Set::SetView x, int xIndex, 
													  vector<EdgeInfo>& updatedEdges); 

		void print() const;
		void printResidual() const; 
		void printBounds(int x) const;
};

#endif