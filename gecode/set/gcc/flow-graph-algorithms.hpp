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

#define INF_INT INT_MAX
#define INF_UINT UINT_MAX
#define NONE_UINT INF_UINT - 1

using namespace std;

class FlowGraphAlgorithms {
	private:
		FlowGraph& graph;

		// Add / update / delete residual edges related to the original graph edge 
		// source->dest. Updates / deletions are needed because in each iteration,
		// instead of building the residual graph from scratch, we modify the 
		// previous one only in the edges that change
		void updateResidualGraph(unsigned int source, unsigned int dest, 
														 NormalEdge edge) {
			unsigned int residualEdgeIndex;
			unsigned int residualBackwardsEdgeIndex;
			ResidualEdge* residualEdgeSearch = graph.getResidualEdge(source, dest, 
																														&residualEdgeIndex);
			ResidualEdge* residualBackwardsEdgeSearch = graph.getResidualEdge(dest, source,
																								 &residualBackwardsEdgeIndex);
			if (edge.flow < edge.upperBound) {
				// Add / update forward residual edge
				graph.setOrCreateResidualEdge(residualEdgeSearch, source, 
																		  ResidualEdge(dest, 
																								   edge.upperBound - edge.flow 
																									 ));
			} else if (residualEdgeSearch != NULL) {
				// Delete forward residual edge that should no longer exist
				auto it = graph.nodeList[source].residualEdgeList.begin() + 
																				 residualEdgeIndex;
				graph.nodeList[source].residualEdgeList.erase(it);
			}

			if (edge.flow > edge.lowerBound) {
				// Add / update backward residual edge
				graph.setOrCreateResidualEdge(residualBackwardsEdgeSearch, dest, 
																		  ResidualEdge(source, 
																									 edge.flow - edge.lowerBound 
																									));
			} else if (residualBackwardsEdgeSearch != NULL) {
				// Delete backward residual edge that should no longer exist
				auto it = graph.nodeList[dest].residualEdgeList.begin() + 
																			 residualBackwardsEdgeIndex;
				graph.nodeList[dest].residualEdgeList.erase(it);
			}
		}

		void sendFlow(pair<unsigned int, unsigned int>& violation, 
								  vector<unsigned int>& shortestPath, 
									unsigned int minUpperBound, LI* li
		) {
			// Send flow through the path edges and update residual graph
			unsigned int prev = violation.first;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				NormalEdge *edge = graph.getEdge(prev, *it);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow += minUpperBound;
					if (edge->destNode < graph.totalVarNodes && li != NULL) {
						cout << "Adding " << (*graph.nodeToVal)[prev] << " to li" << endl; 
						(*li)[edge->destNode].insert((*graph.nodeToVal)[prev]);
					}
					updateResidualGraph(prev, *it, *edge);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(*it, prev);
					edge->flow -= minUpperBound;
					updateResidualGraph(*it, prev, *edge);
				}
				prev = *it;
			}
		}

		unsigned int findMinUpperBound(pair<unsigned int, unsigned int>& violation, 
																	 vector<unsigned int>& shortestPath) {
			// Find min upper bound along shortest path
			unsigned int prev = violation.first;
			unsigned int minUpperBound = INF_UINT;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				// Bellman returns the path in reverse, so traverse it in reverse
				ResidualEdge *edge = graph.getResidualEdge(prev, *it);
				minUpperBound = min(minUpperBound, edge->upperBound);
				prev = *it;
			}
			return minUpperBound;
		}

		bool minCostFlowIteration(pair<unsigned int, unsigned int> violation, LI* li
		) {		
			//graph.printResidual();
			//cout << "Violation " << violation.first << " " << violation.second 
			//		 << endl;
			vector<unsigned int> path;
			if (!findPath(violation.second, violation.first, path)) {
				// Constraint is not consistent
				return false;
			}
			/*for (vector<unsigned int>::reverse_iterator i = path.rbegin(); 
           i != path.rend(); 
					 ++i ) { 			
				cout << *i << (i != path.rend()-1 ? "->" : "");
			}
			cout << endl;*/

			unsigned int minUpperBound = findMinUpperBound(violation, path);		
			sendFlow(violation, path, minUpperBound, li);

			if (li != NULL) {
				for (int i = 0; i < graph.totalVarNodes; i++) {
					for (auto& v: (*li)[i])
						cout << v << endl;
				}
				for (int i = graph.totalVarNodes; i < graph.sNode(); i++) {
					for (auto& e: graph.nodeList[i].edgeList) {
						if (e.flow) {
							cout << "flow for var " << e.destNode << " with val " << (*graph.nodeToVal)[i] << endl;
						}
					}
				}
			}

			return true;
		}

		// Shortest path from source node to dest node
		// Return the shortest path and cost through parameters,
		// or false as return value in case of no path
		bool findPath(unsigned int source, unsigned int dest, 
									vector<unsigned int>& path) const {
			vector<unsigned int> prev;
			DFS(source, dest, prev);
	
			// No path exists
			if (prev[dest] == NONE_UINT) {
				return false;
			}
			// Trace back shortest path
			unsigned int node = dest;
			while (node != source) {
				path.push_back(node);
				node = prev[node];
			}

			path.push_back(node);
			return true;
		}

		void DFS(unsigned int source, unsigned int dest, vector<unsigned int>& prev) 
		const { 
			prev.assign(graph.nodeList.size(), NONE_UINT);
			stack<unsigned int> frontier;
			frontier.push(source);
			while (!frontier.empty()) {
				unsigned int node = frontier.top();
				frontier.pop();
				for (auto& e: graph.nodeList[node].residualEdgeList) {
					if ((prev[e.destNode] != NONE_UINT) || 
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

		bool findMinCostFlow(LI* li
		) {
			pair<unsigned int, unsigned int> violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
				if (!minCostFlowIteration(violation, li)) {
					return false;
				}
			}
			graph.oldFlowIsFeasible = true;
			return true;
		}

		// Given updatedEdges contains the edges whose bounds have been tightened
		// since last execution, do the following:
		// - Update the residual graph to match the changes
		// - If the old flow is not still feasible, find a new one, using the 
		//   incremental algorithm from the publication
		bool updateMinCostFlow(vector<EdgeNodes>& updatedEdges, LI* li) {
			vector<NormalEdge*> edgeReference;
			for (auto& edge: updatedEdges) {
				edgeReference.push_back(graph.getEdge(edge.first, edge.second));
				updateResidualGraph(edge.first, edge.second, *edgeReference.back());
			}
			if (graph.oldFlowIsFeasible) {
				return true;
			}
			for (unsigned int i = 0 ; i < updatedEdges.size(); i++) {
				// Among the edges that changed, look for one violating bounds
				// Assume we are violating lower bound on init
				auto e = edgeReference[i];
				unsigned int src = updatedEdges[i].first;
				unsigned int dest = updatedEdges[i].second;
				if (e->flow >= e->lowerBound && e->flow <= e->upperBound) {
					// All bounds satisfied, try another edge
					continue;
				}
				if (e->flow > e->upperBound) {
					// Violating upper bound, swap direction of initial violating edge
					std::swap(src, dest);
				}
				if (!minCostFlowIteration({src, dest}, li)) {
					return false;
				}
			}
			graph.oldFlowIsFeasible = true;
			return true;
		}
	
	void findOneSCC(unsigned int src, vector<unsigned int>& ids, vector<unsigned int>& low, 
										stack<unsigned int>& localVisited, vector<bool>& onLocalVisited, 
										unsigned int* id, unsigned int* sccCount) const {
			stack<pair<unsigned int, unsigned int>> frontier;
			frontier.push({src, 0});
			localVisited.push(src);
			onLocalVisited[src] = true;
			ids[src] = low[src] = *id;
			(*id)++;
			while (!frontier.empty()) {
				unsigned int node = frontier.top().first;
				unsigned int curEdgeIndex = frontier.top().second;
				auto& edges = graph.nodeList[node].residualEdgeList;
				unsigned int destNode; 
				if (curEdgeIndex < edges.size()) {
					destNode = edges[curEdgeIndex].destNode;
					if (ids[destNode] == NONE_UINT) {
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
						unsigned int tmp = localVisited.top();
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
					destNode = graph.nodeList[node].residualEdgeList[curEdgeIndex].destNode;
					
				  // process here as at, to
					if (onLocalVisited[destNode]) {
						low[node] = min(low[node], low[destNode]);
					}
					
					frontier.pop();
					frontier.push({node, ++curEdgeIndex});
				}
			}
		}


		void findSCC(vector<unsigned int>& scc) const {
			vector<unsigned int> ids;
			vector<bool> onLocalVisited;
			stack<unsigned int> localVisited;
			ids.assign(graph.nodeList.size(), NONE_UINT);
			scc.assign(graph.nodeList.size(), NONE_UINT);
			onLocalVisited.assign(graph.nodeList.size(), false);

			unsigned int id = 0;
			unsigned int sccCount = 0;
			
			for (unsigned int src = 0; src < graph.nodeList.size(); src++) {
				if (ids[src] == NONE_UINT) {
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

		ExecStatus performArcConsistency(Space& home, ViewArray<Set::SetView>& x, 
															       vector<EdgeNodes>& updatedEdges) {
			// Edge nodes, along with the actual value the src node
			// corresponds to 
			struct EdgeWithVal {
				unsigned int src;
				unsigned int dest;
				int val;
				EdgeWithVal(const unsigned int src, const unsigned int dest, const int val)
					: src(src), dest(dest), val(val) {}
			};

			// Hold the edges we decide to prune during arc consistency
			// We do the actual pruning at the end of this function's iterations
			vector<EdgeWithVal> edgesToPrune;
			vector<unsigned int> scc;
			findSCC(scc);

			for (unsigned int n = graph.totalVarNodes; n < graph.sNode(); n++) {
				for (auto& e: graph.nodeList[n].edgeList) {
					if (!e.flow && scc[n] != scc[e.destNode]) {
						int val = (*graph.nodeToVal)[n];
						edgesToPrune.push_back(EdgeWithVal(n, e.destNode, val));
					}
				}
			}


			// Do the actual pruning and update data structures
			for (auto& edge: edgesToPrune) {
				NormalEdge* actualEdge = graph.getEdge(edge.src, edge.dest);
				assert(actualEdge != NULL);
				// Push to updatedEdges so we can modify the residual graph accordingly
				// on the next min cost flow computation
				updatedEdges.push_back(EdgeNodes(edge.src, edge.dest));
				// Prune
				cout << "Prunning val " << edge.val << " from " << edge.dest << endl;
				GECODE_ME_CHECK(x[edge.dest].exclude(home, edge.val));
				// Also remove from varToVals
				auto& vals = graph.varToLub[edge.dest];
				vals.erase(edge.val);
				// Update upper bound
				actualEdge->upperBound = 0;
				assert(!actualEdge->flow);
				if (x[edge.dest].assigned()) {
					// If a variable got assigned by pruning, set corresponding edge
					// lower bound to 1
					for (SetVarGlbValues v(x[edge.dest]); v(); ++v) {
						auto valNode = graph.valToNode->find(v.val())->second;
						graph.getEdge(valNode, edge.dest)->lowerBound = 1;
					}
				}
			}
			return ES_OK;
		}

};

#endif