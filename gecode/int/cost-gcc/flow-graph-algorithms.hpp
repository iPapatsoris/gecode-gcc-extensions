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
// unsigned int countt = 0;

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
																								   edge.upperBound - edge.flow, 
																									 edge.cost));
			} else if (residualEdgeSearch != NULL) {
				// Delete forward residual edge that should no longer exist
				auto it = graph.nodeList[source].residualEdgeList.begin() + 
																				 residualEdgeIndex;
				graph.nodeList[source].residualEdgeList.erase(it);
			//	graph.orderGraph.removeEdge(source, dest);
			}

			if (edge.flow > edge.lowerBound) {
				// Add / update backward residual edge
				graph.setOrCreateResidualEdge(residualBackwardsEdgeSearch, dest, 
																		  ResidualEdge(source, 
																									 edge.flow - edge.lowerBound, 
																									 -edge.cost));
				//graph.orderGraph.addEdge(dest, source);
			} else if (residualBackwardsEdgeSearch != NULL) {
				// Delete backward residual edge that should no longer exist
				auto it = graph.nodeList[dest].residualEdgeList.begin() + 
																			 residualBackwardsEdgeIndex;
				graph.nodeList[dest].residualEdgeList.erase(it);
			}
		}

		// Bellman-Ford algorithm for shortest paths with negative costs.
		// If dest is not NULL, ignore any direct source->dest edge.
		// This is needed when searching for shortest path to a specific 
		// destination, by the min cost flow algorithm.
     void bellmanFordShortestPaths(unsigned int source, 
																	vector<unsigned int>& prev, vector<int>& dist, 
																	unsigned int* dest = NULL) const {
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist.assign(graph.nodeList.size(), INF_INT);
			dist[source] = 0;

			vector<unsigned int> updatedNodes1 = {source};
			vector<unsigned int> updatedNodes2;
			vector<unsigned int>* updatedNodesOld = &updatedNodes1;
			vector<unsigned int> *updatedNodesNew = &updatedNodes2;
			while (!updatedNodesOld->empty()) {
				bool foundUpdate = false;
				for (auto node: *updatedNodesOld) {
					for (auto &edge : graph.nodeList[node].residualEdgeList) {
						if ((dest == NULL || 
							!(node == source && edge.destNode == *dest)) && 
							dist[node] != INF_INT && dist[node] + edge.cost < 
																			 dist[edge.destNode]) {
							// Ignore direct source->dest edge when looking for shortest path 
							// to specific requested destination
							dist[edge.destNode] = dist[node] + edge.cost;
							prev[edge.destNode] = node;
							foundUpdate = true;
							updatedNodesNew->push_back(edge.destNode);
						}
					}
				}
				if (!foundUpdate) {
					// No updates between two iterations; early termination
					break;
				}
				swap(updatedNodesNew, updatedNodesOld);
				updatedNodesNew->clear();
			}
		}

		/*
		void updatePotentials(vector<unsigned int>& visitedNodes, vector<unsigned int>& dist, int pathCost) {
			for (unsigned int n = 0; n < visitedNodes.size(); n++) {
				if (visitedNodes[n]) {
					//assert(dist[n] != INF_UINT);
					graph.nodeList[n].potential += (pathCost - dist[n]);
				}
			}
		}

		void updateCosts() {
		for (unsigned int i = 0; i < graph.nodeList.size(); i++) { // can opt to not iterate all
				for (auto& e: graph.nodeList[i].residualEdgeList) {
												// SIGOURA cost edw kai OXI reduced. alliws sto paper bgianoun arnitika
					e.reducedCost = (e.isBackwards ? 0 : e.cost - graph.nodeList[i].potential + graph.nodeList[e.destNode].potential);
				}
			}
		}*/

			void sendFlow(pair<unsigned int, unsigned int>& violation, vector<int>& shortestPath, unsigned int minUpperBound, LI& li) {
				// Send flow through the path edges and update residual graph
			unsigned int prev = violation.first;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				NormalEdge *edge = graph.getEdge(prev, *it);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow += minUpperBound;
					graph.flowCost += edge->cost;
					if (edge->destNode < graph.totalVarNodes) {
						li[edge->destNode] = (*graph.nodeToVal)[prev];
					}
					updateResidualGraph(prev, *it, *edge);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(*it, prev);
					edge->flow -= minUpperBound;
					if (!edge->flow) {
						graph.flowCost -= edge->cost;
					}
					updateResidualGraph(*it, prev, *edge);
				}
				prev = *it;
			}
		}

		unsigned int findMinUpperBound(pair<unsigned int, unsigned int>& violation, vector<int>& shortestPath) {
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

		bool minCostFlowIteration(pair<unsigned int, unsigned int> violation, LI& li) {		
			vector<int> shortestPath, dist;
			int pathCost; 
			if (!findShortestPathNegativeCosts(violation.second, violation.first, 
																				 shortestPath, dist, pathCost)) {
				// Constraint is not consistent
				return false;
			}

		//	updatePotentials(visitedNodes, dist, pathCost);
		//	updateCosts();

			unsigned int minUpperBound = findMinUpperBound(violation, shortestPath);		
			sendFlow(violation, shortestPath, minUpperBound, li);
			return true;
		}

		// Shortest path from source node to dest node
		// Return the shortest path and cost through parameters,
		// or false as return value in case of no path
		bool findShortestPathNegativeCosts(unsigned int source, unsigned int dest, 
																			 vector<int>& path, vector<int>& dist,
																			int& cost) 
																			 const {
			vector<unsigned int> prev;
			
			bellmanFordShortestPaths(source, prev, dist, &dest);
	
			// No path exists
			if (dist[dest] == INF_INT) {
				return false;
			}
			// Trace back shortest path
			unsigned int node = dest;
			while (node != source) {
				path.push_back(node);
				node = prev[node];
			}

			path.push_back(node);
			cost = dist[dest];
			return true;
		}

		// Disjktra's algorithm for shortest paths.
		// Early termination when we have found paths to targetNodes, we do not 
		// compute the rest of the nodes.
		// Nodes are removed from targetNodes as we finalize their paths. 
		// Normally, targetNodes would be empty when function returns.
		// Due to another early termination condition involving costLowerBound, 
		// it is possible to terminate without visiting all targetNodes.
		// In this case, the nodes remaining in targetNodes are inconsistent 
		// and can be pruned.
		void findShortestPathsReducedCosts(unsigned int source, 
																			 unordered_set<unsigned int>& targetNodes, 
																			 vector<unsigned int>& dist, 
																			 unsigned int costLowerBound) const {
			vector<bool> visited;
			visited.assign(graph.nodeList.size(), false);
			dist.assign(graph.nodeList.size(), INF_UINT);
			dist[source] = 0;
			if (targetNodes.empty()) {
				return;
			}

			struct HeapItem {
				unsigned int node;
				unsigned int dist;

				HeapItem(unsigned int node, unsigned int dist) : node(node), dist(dist) 
				{}
				HeapItem() {}
				// Overload so we can use STL min-heap
				int operator() (const HeapItem& a, const HeapItem& b) { 
					return a.dist > b.dist; 
				}
			};		

			priority_queue<HeapItem, vector<HeapItem>, HeapItem> heap;
			heap.push(HeapItem(source, 0));

			while (!heap.empty()) {
				struct HeapItem curItem = heap.top();
				heap.pop();
				unsigned int node = curItem.node;
				visited[node] = true;

				// Early termination. If condition holds for current node, it will also 
				// hold for every unexplored one. This condition is enough to prune the 
				// associated values of the remaining unexplored targetNodes early,
				// without having to find the remaining shortest paths
				if (dist[node] > costLowerBound) {
					return;
				}

				// Stop when we find paths to all targetNodes of interest
				targetNodes.erase(node);
				if (targetNodes.empty()) {
					return;
				}

				if (dist[node] < curItem.dist) {
					// We have already found a better path than the one this heap item 
					// suggests
					continue;
				}

				for (auto& edge: graph.nodeList[node].residualEdgeList) {
					if (visited[edge.destNode]) {
						continue;
					}
					unsigned int newDist = dist[node] + edge.reducedCost;;
					if (newDist < dist[edge.destNode]) {
						// Found path of lower cost
						dist[edge.destNode] = newDist;
						heap.push(HeapItem(edge.destNode, newDist));
					}
				}
			}
		}


		// Optimization to skip finding shortest paths for GAC, according to 
		// Practical Improvements section of the research paper
		bool earlyPrune(unsigned int a, unsigned int b, unsigned int y, 
										unsigned int m) const {

			auto mFactor = [](unsigned int b, unsigned int y, FlowGraph& graph) {
				unsigned int min = INF_UINT;
				for (unsigned int z = 0 ; z < graph.totalVarNodes ; z++) {
					ResidualEdge *edgeZB;
					if (z == y || ((edgeZB = graph.getResidualEdge(z, b)) == NULL)) {
						continue;
					}
					for (auto& edgeZC: graph.nodeList[z].residualEdgeList) {
						min = std::min(min, edgeZB->reducedCost + edgeZC.reducedCost);
					}
				}
				return min;
			};

			NormalEdge* edgeSA = graph.getEdge(graph.sNode(), a);
			NormalEdge* edgeSB = graph.getEdge(graph.sNode(), b);
			ResidualEdge* edgeAY = graph.getResidualEdge(a, y);
			ResidualEdge* edgeYB = graph.getResidualEdge(y, b);

			if (edgeSA->flow < edgeSA->upperBound 
					&& edgeSB->flow > edgeSB->lowerBound 
					&& edgeAY->reducedCost > m - edgeYB->reducedCost) {
		//		cout << "\tCONDITION 1 EARLY PRUNNING " << a << " from " << y << endl; 
				return true;
			} 
			unsigned int mB = mFactor(b, y, graph);
			if (mB != INF_UINT && edgeSA->flow < edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mB > m) {
		//		cout << "\tCONDITION 2 EARLY PRUNNING " << a << " from " << y << endl;
				cout << edgeSA->flow << " " << edgeSA->upperBound << " " << edgeSB->flow << " " << 
				edgeSB->lowerBound << " " << edgeAY->reducedCost << " " << mB << " " << m << endl; 
				return true;
			}
			unsigned int mA = mFactor(a, y, graph);
			if (mA != INF_UINT && edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow > edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA > m - edgeYB->reducedCost) {
		//		cout << "\tCONDITION 3 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			if (mA != INF_UINT && mB != INF_UINT && edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA + mB > m) {
		//		cout << "\tCONDITION 4 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			return false;
		}

	public:
		FlowGraphAlgorithms(FlowGraph& graph) : graph(graph) {}

		bool findMinCostFlow(LI& li) {
			pair<unsigned int, unsigned int> violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
				if (!minCostFlowIteration(violation, li)) {
					return false;
				}
			}
			graph.oldFlowIsFeasible = true;
			//graph.calculateFlowCost(li);
			return graph.checkFlowCost();
		}

		// Given updatedEdges contains the edges whose bounds have been tightened
		// since last execution, do the following:
		// - Update the residual graph to match the changes
		// - If the old flow is not still feasible, find a new one, using the 
		//   incremental algorithm from the publication
		bool updateMinCostFlow(vector<EdgeNodes>& updatedEdges, LI& li) {
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
			return graph.checkFlowCost();
		}
	
		// In addition to pruning, hold the affected Val->Var edges in updatedEdges
		// so we can update the residual graph later accordingly. We do not update
		// it here, because in case the home space fails due to another
		// constraint, or if we find a solution from this pruning, we would have 
		// updated it for no reason, as we wouldn't need to re-check the validity of 
		// costgcc, the search would backtrack to previous instances.
		// The reason why 

		ExecStatus performArcConsistency(Space& home, ViewArray<Int::IntView>& vars, 
															       vector<EdgeNodes>& updatedEdges) {
			graph.addTResidualEdges();
			vector<int> distances;
			vector<unsigned int> prev;
			bellmanFordShortestPaths(graph.tNode(), prev, distances, NULL);
		  graph.calculateReducedCosts(distances);
 
			// Edge nodes, along with the actual value the src node
			// corresponds to 
			struct EdgeWithVal {
				unsigned int src;
				unsigned int dest;
				int val;
				EdgeWithVal(const unsigned int src, const unsigned int dest, const int val)
					: src(src), dest(dest), val(val) {}
			};

			struct MaxRegret {
				unsigned int var = NONE_UINT;
				unsigned int val = NONE_UINT;
				unsigned int regret = 0;
			};

			struct MinDist {
				unsigned int bestVal = NONE_UINT;
				unsigned int val = NONE_UINT;
				unsigned int dist = INF_UINT;
			};

			// Hold the edges we decide to prune during arc consistency
			// We do the actual pruning at the end of this function's iterations
			vector<EdgeWithVal> edgesToPrune;

			// Gather the targetNodes we want to find shortests paths to from B,
			// and check early prune conditions to skip finding some
			for (auto& edge: graph.nodeList[graph.sNode()].edgeList) {
				if (edge.flow > 0) {
					vector<unsigned int> yList;
					unsigned int b = edge.destNode;
					unordered_set<unsigned int> targetNodes;
					vector< pair<unsigned int, unsigned int>> ayList;
					for (auto& edgeBY: graph.nodeList[b].edgeList) {
						if (edgeBY.flow == 1) {
							unsigned int y = edgeBY.destNode;
						//	minDist[y].bestVal = b;
							for (IntVarValues v(vars[y]); v(); ++v) {
								unsigned int a = (*graph.valToNode)[v.val()];
								if (a != b) {
									if (earlyPrune(a, b, y, graph.costUpperBound - 
															   graph.flowCost)) {
										edgesToPrune.push_back(EdgeWithVal(a, y, v.val()));
										continue;
									}
									targetNodes.insert(a);
									ayList.push_back({a, y});
								}
							}
						}
					}

					vector<unsigned int> reducedDistances;
					if (targetNodes.empty()) {
						continue;
					}
					findShortestPathsReducedCosts(b, targetNodes, reducedDistances, 
																				graph.costUpperBound - graph.flowCost);
					
					for (const auto& ay: ayList) {
						const auto a = ay.first;
						const auto y = ay.second;
						if (targetNodes.find(a) != targetNodes.end()) {
							// Shortest paths function normally removes targetNodes as it
							// computes them. For any targetNodes that remain, we know for 
							// sure that they can be pruned without checking shortest paths,
							// because we have explored a node before them with
							// (dist > costUpperBound - flowCost).
							// Since all reduced costs are non-negative, we already know
							// that (dist > costUpperBound - flowCost - costAY - costYB)
					//		cout << "\tEARLY PRUNING " << a << " FROM " << y << endl;
							edgesToPrune.push_back(EdgeWithVal(a, y, 
																						  graph.nodeToVal->find(a)->second));
							continue;
						}
						ResidualEdge *residualEdge = graph.getResidualEdge(a, y);
						unsigned int costAY = residualEdge->reducedCost;
						unsigned int costYB = graph.getResidualEdge(y, b)->reducedCost;
						if ((int)reducedDistances[a] > (graph.costUpperBound - graph.flowCost 
																		  - (int)costAY - (int)costYB)) {
							edgesToPrune.push_back(EdgeWithVal(a, y,
																							graph.nodeToVal->find(a)->second));
						} 
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
				GECODE_ME_CHECK(vars[edge.dest].nq(home, edge.val));
				// Also remove from varToVals
				auto& vals = graph.varToVals.map.find(edge.dest)->second;
				vals.erase(edge.val);
				// Update upper bound
				actualEdge->upperBound = 0;
				assert(!actualEdge->flow);
				if (vars[edge.dest].assigned()) {
					// If a variable got assigned by pruning, set corresponding edge
					// lower bound to 1
					int assignedVal = vars[edge.dest].val();
					assert(*vals.begin() == assignedVal);
					auto valNode = graph.valToNode->find(assignedVal)->second;
					graph.getEdge(valNode, edge.dest)->lowerBound = 1;
				}
			}

    	graph.removeTResidualEdges();
			return ES_OK;
		}

};

#endif