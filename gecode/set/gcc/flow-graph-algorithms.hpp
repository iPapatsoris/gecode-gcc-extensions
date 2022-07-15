#ifndef FLOW_GRAPH_ALGORITHMS
#define FLOW_GRAPH_ALGORITHMS

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <climits>
#include "flow-graph.hpp"

#define INF INT_MAX
#define NONE INF - 1

using namespace std;

class FlowGraphAlgorithms {
	private:
		FlowGraph& graph;

		// Add / update / delete residual edges related to the original graph edge 
		// source->dest. Updates / deletions are needed because in each iteration,
		// instead of building the residual graph from scratch, we modify the 
		// previous one only in the edges that change
		void updateResidualGraph(const ViewArray<Set::SetView>& x, int src, 
														 int dest, 
														 NormalEdge edge) {
			auto& nodeList = graph.backtrackStable->nodeList;	
			int residualEdgeIndex;
			int residualBackwardsEdgeIndex;
			ResidualEdge* residualEdgeSearch = graph.getResidualEdge(
																							 	src, 
																								dest, 
																								&residualEdgeIndex);
			ResidualEdge* residualBackwardsEdgeSearch = graph.getResidualEdge(
																									dest,
																									src,
																								 	&residualBackwardsEdgeIndex);
			int lowerBound = edge.lowerBound;
			int upperBound = edge.upperBound;
			if (dest == graph.tNode()) {
				lowerBound = x[src].cardMin();
				upperBound = x[src].cardMax();
			}																						 
			if (edge.flow < upperBound) {
				// Add / update forward residual edge
				graph.setOrCreateResidualEdge(residualEdgeSearch, src, 
																		  ResidualEdge(dest, 
																								   upperBound - edge.flow 
																									 ));
			} else if (residualEdgeSearch != NULL) {
				// Delete forward residual edge that should no longer exist
				auto it = nodeList[src].residualEdgeList.begin() + 
																					residualEdgeIndex;
					graph.backtrackStable->nodeList[src].residualEdgeList.erase(it);
			}

			if (edge.flow > lowerBound) {
				// Add / update backward residual edge
				graph.setOrCreateResidualEdge(residualBackwardsEdgeSearch, dest, 
																		  ResidualEdge(src, 
																									 edge.flow - lowerBound 
																									));
			} else if (residualBackwardsEdgeSearch != NULL) {
				// Delete backward residual edge that should no longer exist
				auto it = nodeList[dest].residualEdgeList.begin() + 
																			 residualBackwardsEdgeIndex;
				nodeList[dest].residualEdgeList.erase(it);
			}
		}

		// Clear residual graph and build it again. Do not clear T->Var edges.
		// Since we iterate through the graph, this step is a good opportunity
		// to also update BestBranch with the latest flow assignment.
		void buildResidualGraph(const ViewArray<Set::SetView>& x, 
													  BestBranch *bestBranch) {
			auto& nodeList = graph.backtrackStable->nodeList;
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				nodeList[i].residualEdgeList.clear();
				if (bestBranch != NULL && i < (unsigned int) graph.totalVarNodes) {
					(*bestBranch)[i].clear();
				}
			}

			for (unsigned int i = 0; i < nodeList.size(); i++) {
				auto& node = nodeList[i];
				for (int j = 0; j < graph.edgeListSize[i]; j++) {
					auto& edge = (node.edgeList.list)[j];
					int lowerBound = edge.lowerBound;
					int upperBound = edge.upperBound;
					if (edge.destNode == graph.tNode()) {
						lowerBound = x[i].cardMin();
						upperBound = x[i].cardMax();
					}	else if (edge.destNode < graph.totalVarNodes 
										 && x[edge.destNode].contains(
										 graph.backtrackStable->nodeToVal[i])) {
						lowerBound = 1; 
					}																					 
					if (edge.flow < upperBound) {
						node.residualEdgeList.push_back(ResidualEdge(
																						  edge.destNode, 
																						  upperBound - edge.flow));
					}
					if (edge.flow > lowerBound) {
						nodeList[edge.destNode].residualEdgeList.push_back(
																						 ResidualEdge(
																							 i, 
																							 edge.flow - lowerBound));
					}
					if (bestBranch != NULL && edge.flow && edge.destNode < 
																								 graph.totalVarNodes) {
						// If edge is of type Val->Var and has flow, update BestBranch
						(*bestBranch)[edge.destNode].insert(graph.backtrackStable->nodeToVal[i]);
					}
				}
			}
		}

		void sendFlow(const ViewArray<Set::SetView>& x, const EdgeInfo& violation, vector<int>& shortestPath, 
								  BestBranch* bestBranch) {
			// Send flow through the path edges and update residual graph
			int prev = violation.src;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				NormalEdge *edge = graph.getEdge(prev, *it);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow++;
					if (edge->destNode < graph.totalVarNodes && bestBranch != NULL) {
					//	cout << "Adding " << (*graph.nodeToVal)[prev] << " to bestBranch" << endl; 
						(*bestBranch)[edge->destNode].insert(graph.backtrackStable->nodeToVal[prev]);
					}
					updateResidualGraph(x, prev, *it, *edge);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(*it, prev);
					edge->flow--;
					if (edge->destNode < graph.totalVarNodes && bestBranch != NULL) {
						(*bestBranch)[edge->destNode].erase(graph.backtrackStable->nodeToVal[prev]);
					}
					updateResidualGraph(x, *it, prev, *edge);
				}
				prev = *it;
			}
		}

		bool minCostFlowIteration(const ViewArray<Set::SetView>& x, 
															const EdgeInfo& violation, BestBranch* bestBranch
		) {		
			//graph.printResidual();
			// if (graph.debug) cout << "Violation " << violation.src << " " << violation.dest 
			// 		 << endl;
			vector<int> path;
			if (!findPath(violation.dest, violation.src, path)) {
				// Constraint is not consistent
				// if (graph.debug) cout << "false lol" << endl; 
				return false;
			}
			// if (graph.debug) {
			// for (vector<int>::reverse_iterator i = path.rbegin(); 
      //      i != path.rend(); 
			// 		 ++i ) { 			
			// 	cout << *i << (i != path.rend()-1 ? "->" : "");
			// }
			// cout << endl;
			// exit(1);
			// }

			sendFlow(x, violation, path, bestBranch);

		/*	if (bestBranch != NULL) {
				for (int i = 0; i < graph.totalVarNodes; i++) {
					for (auto& v: (*bestBranch)[i])
						cout << v << endl;
				}
				for (int i = graph.totalVarNodes; i < graph.sNode(); i++) {
					for (auto& e: graph.nodeList[i].edgeList) {
						if (e.flow) {
							cout << "flow for var " << e.destNode << " with val " << (*graph.nodeToVal)[i] << endl;
						}
					}
				}
			}*/

			return true;
		}

		// Shortest path from source node to dest node
		// Return the shortest path and cost through parameters,
		// or false as return value in case of no path
		bool findPath(int source, int dest, 
									vector<int>& path) const {
			vector<int> prev;
			DFS(source, dest, prev);
	
			// No path exists
			if (prev[dest] == NONE) {
				return false;
			}
			// Trace back shortest path
			int node = dest;
			while (node != source) {
				path.push_back(node);
				node = prev[node];
			}

			path.push_back(node);
			return true;
		}

		void DFS(int source, int dest, vector<int>& prev) 
		const { 
			auto& nodeList = graph.backtrackStable->nodeList;
			prev.assign(nodeList.size(), NONE);
			stack<int> frontier;
			frontier.push(source);
			while (!frontier.empty()) {
				int node = frontier.top();
				frontier.pop();
				for (auto& e: nodeList[node].residualEdgeList) {
					if ((prev[e.destNode] != NONE) || 
							(node == source && e.destNode == dest)) {
						continue;
					}
					prev[e.destNode] = node;					
					if (e.destNode == dest) {
						return;
					}
					frontier.push(e.destNode);
				}
			} 
		}

		void BFS(int source, int dest, vector<int>& prev) const { 
			auto& nodeList = graph.backtrackStable->nodeList;
			prev.assign(nodeList.size(), NONE);
			queue<int> frontier;
			frontier.push(source);
			while (!frontier.empty()) {
				int node = frontier.front();
				frontier.pop();
				for (auto& e: nodeList[node].residualEdgeList) {
					if ((prev[e.destNode] != NONE) || 
							(node == source && e.destNode == dest)) {
						continue;
					}
					prev[e.destNode] = node;					
					if (e.destNode == dest) {
						return;
					}
					frontier.push(e.destNode);
				}
			} 
		}

	public:
		FlowGraphAlgorithms(FlowGraph& graph) : graph(graph) {}

		bool findMinCostFlow(const ViewArray<Set::SetView>& x, BestBranch* bestBranch) {
			EdgeInfo violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
				if (!minCostFlowIteration(x, violation, bestBranch)) {
					return false;
				}
			}
			// graph.print();
			return true;
		}

		// Given updatedEdges contains the edges whose bounds have been tightened
		// since last execution, do the following:
		// - Update the residual graph to match the changes
		// - If the old flow is not still feasible, find a new one, using the 
		//   incremental algorithm from the publication
		bool updateMinCostFlow(const ViewArray<Set::SetView>& x, vector<EdgeInfo>& updatedEdges, BestBranch* bestBranch) {
			vector<NormalEdge*> edgeReference;
			buildResidualGraph(x, bestBranch);
			assert(updatedEdges.size());
			if (graph.debug) graph.printResidual();
			// Repair upper bound violations
			for (auto& e: updatedEdges) {
				auto res = graph.getEdge(e.src, e.dest);
				auto& var = (e.dest == graph.tNode() ? x[e.src] : x[e.dest]);
				if (res == NULL || 
					 (e.lowerBoundViolation && res->flow >= e.violationBound) ||
					 (!e.lowerBoundViolation && res->flow <= e.violationBound)) {					// Check if edge has already been removed because it existed twice 
					if (graph.debug) cout << "skipping" << endl; 
					// inside updatedEdges.
					// It is possible to have duplicates inside updatedEdges, because
					// when the propagator is scheduled, it is not necessary that it 
					// will execute right away. It is possible that another propagator
					// may take precedence, removing some values, and then costgcc
					// will check again for domain changes. If the change is on a variable
					// that has already been inserted in updatedEdges but has not been 
					// processed yet, then the old variable-value pair will be inserted
					// again, along with the new one. This could be fixed by using a set
					// instead of vector and check for existence before inserting, but 
					// it would be more expensive. So we just allow duplicates sometimes
					// and skip them.
					continue;
				}
				do {
					int src = e.src;
					int dest = e.dest;
					if (!e.lowerBoundViolation) {
						swap(src, dest);
					}
					if (!minCostFlowIteration(x, {src, dest}, bestBranch)) {
						return false;
					}
					if (!e.lowerBoundViolation && e.dest != graph.tNode()) {
						graph.deleteEdge(e.src, e.dest);
						// We have already removed flow from this Val->Var edge and updated 
						// residual graph, so we know for sure that the only residual edge 
						// direction will be the same as the original edge.
						graph.deleteResidualEdge(e.src, e.dest);
						int val = graph.backtrackStable->nodeToVal[e.src];
						graph.backtrackStable->varToVals[e.dest].deleteVal(
																										val, 
																										&graph.varToValsSize[e.dest]
																										);
					}
				} while (e.dest == graph.tNode() && !res->isFeasible(var));
			}
			return true;
		}
	
		void findOneSCC(int src, vector<int>& ids, vector<int>& low, 
										stack<int>& localVisited, vector<bool>& onLocalVisited, 
										int* id, int* sccCount) const {
			stack<pair<int, int>> frontier;
			auto& nodeList = graph.backtrackStable->nodeList;
			frontier.push({src, 0});
			localVisited.push(src);
			onLocalVisited[src] = true;
			ids[src] = low[src] = *id;
			(*id)++;
			while (!frontier.empty()) {
				int node = frontier.top().first;
				int curEdgeIndex = frontier.top().second;
				auto& edges = nodeList[node].residualEdgeList;
				int destNode; 
				if (curEdgeIndex < (int) edges.size()) {
					destNode = edges[curEdgeIndex].destNode;
					if (ids[destNode] == NONE) {
						// start of call here
						localVisited.push(destNode);
						onLocalVisited[destNode] = true;
						ids[destNode] = low[destNode] = *id;
						(*id)++;
						frontier.push({destNode, 0});
					} else {
						if (onLocalVisited[destNode]) {
							low[node] = min(low[node], low[destNode]);
						}
						frontier.pop();
						frontier.push({node, ++curEdgeIndex});
					}
					continue;
				}

				// last part here
				if (ids[node] == low[node]) {
					while (!localVisited.empty()) {
						int tmp = localVisited.top();
						localVisited.pop();
						onLocalVisited[tmp] = false;
						low[tmp] = ids[node];
						if (tmp == node) {
							break;
						}
					}
					(*sccCount)++;
				}

				frontier.pop();
				if (!frontier.empty()) {
					node = frontier.top().first;
					curEdgeIndex = frontier.top().second;
					destNode = nodeList[node].residualEdgeList[curEdgeIndex].destNode;
					
				  // process here as at, to
					if (onLocalVisited[destNode]) {
						low[node] = min(low[node], low[destNode]);
					}
					
					frontier.pop();
					frontier.push({node, ++curEdgeIndex});
				}
			}
		}


		void findSCC(vector<int>& scc) const {
			auto& nodeList = graph.backtrackStable->nodeList;
			vector<int> ids;
			vector<bool> onLocalVisited;
			stack<int> localVisited;
			ids.assign(nodeList.size(), NONE);
			scc.assign(nodeList.size(), NONE);
			onLocalVisited.assign(nodeList.size(), false);

			int id = 0;
			int sccCount = 0;
			
			for (unsigned int src = 0; src < nodeList.size(); src++) {
				if (ids[src] == NONE) {
					findOneSCC(src, ids, scc, localVisited, onLocalVisited, &id, &sccCount);
				}
			}
			
			/*for (unsigned int i = 0; i < graph.nodeList.size(); i++) {
				cout << "Node " << i << " in SCC " << low[i] << "\n";
			}*/
		}
	
		// In addition to pruning, hold the affected Val->Var edges in updatedEdges
		// so we can update the residual graph later accordingly. We do not update
		// it here, because in case the home space fails due to another
		// constraint, or if we find a solution from this pruning, we would have 
		// updated it for no reason, as we wouldn't need to re-check the validity of 
		// costgcc, the search would backtrack to previous instances.
		// The reason why 

		ExecStatus performArcConsistency(Space& home, ViewArray<Set::SetView>& x) {
			// Edge nodes, along with the actual value the src node
			// corresponds to 
			struct EdgeWithVal {
				int src;
				int dest;
				int val;
				EdgeWithVal(const int src, const int dest, const int val)
					: src(src), dest(dest), val(val) {}
			};

			// Hold the edges we decide to prune during arc consistency
			// We do the actual pruning at the end of this function's iterations
			vector<EdgeWithVal> edgesToPrune;
			vector<int> scc;
			findSCC(scc);
			auto& nodeList = graph.backtrackStable->nodeList;

			for (int n = graph.totalVarNodes; n < graph.sNode(); n++) {
				for (int i = 0; i < graph.edgeListSize[n]; i++) {
					auto& e = nodeList[n].edgeList.list[i];
					if (!e.flow && scc[n] != scc[e.destNode]) {
						int val = (graph.backtrackStable->nodeToVal)[n];
						edgesToPrune.push_back(EdgeWithVal(n, e.destNode, val));
					}
				}
			}

			// Do the actual pruning and update data structures
			for (auto& edge: edgesToPrune) {
				GECODE_ME_CHECK(x[edge.dest].exclude(home, edge.val));
			}
			return ES_OK;
		}
};

#endif