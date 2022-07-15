#ifndef FLOW_GRAPH_ALGORITHMS
#define FLOW_GRAPH_ALGORITHMS

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <climits>
#include <chrono>
#include "flow-graph.hpp"

#define INF_INT INT_MAX
#define NONE_INT INF_INT - 1

using namespace std;
class FlowGraphAlgorithms {
	private:
		FlowGraph& graph;

		// Add / update / delete residual edges related to the original graph edge 
		// source->dest.
		void updateResidualGraph(int source, int dest, NormalEdge edge) {
			auto& nodeList = graph.backtrackStable->nodeList;
			int residualEdgeIndex;
			int residualBackwardsEdgeIndex;
			ResidualEdge* residualEdgeSearch = graph.getResidualEdge(source, dest, 
																														&residualEdgeIndex);
			ResidualEdge* residualBackwardsEdgeSearch = graph.getResidualEdge(
																									 dest, 
																									 source,
																								   &residualBackwardsEdgeIndex);
			if (edge.flow < edge.upperBound) {
				// Add / update forward residual edge
				graph.setOrCreateResidualEdge(residualEdgeSearch, source, 
																		  ResidualEdge(dest, 
																								   edge.upperBound - edge.flow, 
																									 edge.cost));
			} else if (residualEdgeSearch != NULL) {
				// Delete forward residual edge that should no longer exist
				auto it = nodeList[source].residualEdgeList.begin() + residualEdgeIndex;
				nodeList[source].residualEdgeList.erase(it);
			}

			if (edge.flow > edge.lowerBound) {
				// Add / update backward residual edge
				graph.setOrCreateResidualEdge(residualBackwardsEdgeSearch, dest, 
																		  ResidualEdge(source, 
																									 edge.flow - edge.lowerBound, 
																									 -edge.cost));
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
		void buildResidualGraph(BestBranch *bestBranch) {
			auto& nodeList = graph.backtrackStable->nodeList;
			for (int i = 0; i < graph.tNode(); i++) {
				nodeList[i].residualEdgeList.clear();
			}

			for (int i = graph.backtrackStable->totalVarNodes; i < graph.tNode(); i++) {
				auto& node = nodeList[i];
				for (int j = 0; j < graph.edgeListSize[i]; j++) {
					auto& edge = (node.edgeList.list)[j];
					if (edge.flow < edge.upperBound) {
						node.residualEdgeList.push_back(ResidualEdge(
																						  edge.destNode, 
																						  edge.upperBound - edge.flow, 
																							edge.cost));
					}
					if (edge.flow > edge.lowerBound) {
						nodeList[edge.destNode].residualEdgeList.push_back(
																						 ResidualEdge(
																							 i, 
																							 edge.flow - edge.lowerBound, 
																							 -edge.cost));
					}
					if (bestBranch != NULL && edge.flow && edge.destNode < 
																				 graph.backtrackStable->totalVarNodes) {
						// If edge is of type Val->Var and has flow, update BestBranch
						(*bestBranch)[edge.destNode] = graph.backtrackStable->nodeToVal[i];
					}
				}
			}
		}

		// Bellman-Ford algorithm for shortest paths with negative costs.
		// Based on the Shortest Path Faster Algorithm optimization.
		// If isCycle is NULL, it does NOT check for cycles. If it is not, it puts
		// the first cycle that it finds in "cycle" vector parameter.
		// If the graph can contain a cycle, always call with isCycle not NULL, 
		// because otherwise it will run infinitely.
		// If dest is not NULL, ignore any direct source->dest edge.
		// This is needed when searching for shortest path to a specific 
		// destination, by the min cost flow algorithm. When searching without 
		// a specific destination, pass as NULL.
    // Under a cycle, normally the distance of the nodes within and reachable
		// by it would become infinite. We keep track of the length of the 
		// shortest path to each node. Under no cycles, the maximum possible
		// would be equal to all the graph nodes. So if we exceed that limit for a
		// node, it means that there is a cycle. Places it in the appropriate 
		// parameter and also sets isCycle to true.
		void bellmanFordShortestPaths(int source, 
																	vector<int>& prev, vector<int>& dist,
																	int* dest = NULL,
																	vector<int> *cycle = NULL, 
																	bool *isCycle = NULL
																	) const {
			auto& nodeList = graph.backtrackStable->nodeList;
			prev.assign(nodeList.size(), NONE_INT);
			dist.assign(nodeList.size(), INF_INT);
			dist[source] = 0;
			vector<int> len;
			if (isCycle != NULL) {
				len.assign(nodeList.size(), 0);
			}
			vector<int> updatedNodes1 = {source};
			vector<int> updatedNodes2;
			vector<int>* updatedNodesOld = &updatedNodes1;
			vector<int> *updatedNodesNew = &updatedNodes2;
			while (!updatedNodesOld->empty()) {
				bool foundUpdate = false;
				for (auto node: *updatedNodesOld) {
					for (auto &edge : nodeList[node].residualEdgeList) {
						if ((dest == NULL || 
							!(node == source && edge.destNode == *dest)) && 
							dist[node] != INF_INT && dist[node] + edge.cost < 
																			 dist[edge.destNode]) {
							// Ignore direct source->dest edge when looking for shortest path 
							// to specific requested destination
							dist[edge.destNode] = dist[node] + edge.cost;
							prev[edge.destNode] = node;
							if (isCycle != NULL && 
									++len[edge.destNode] == (int) nodeList.size()) {
								*isCycle = true;
								traceCycle(prev, edge.destNode, cycle);
								return;
							}
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

		// Given that a cycle containing node has been detected, find cycle nodes 
		void traceCycle(const vector<int>& prev, int node, vector<int> *cycle) 
			const {
			unordered_set<int> stackContent;
			stack<int> stack;

			while (stackContent.find(node) == stackContent.end()) {
				stack.push(node);
				stackContent.insert(node);
				node = prev[node];
			}
			while (stack.top() != node) {
				cycle->push_back(stack.top());
				stack.pop();
			}
			cycle->push_back(node);
		}

		// Send flow through the edge and update residual graph.
		// Also update BestBranch with latest flow changes
		void sendEdgeFlow(int *prev, int cur, BestBranch *bestBranch) {
			NormalEdge *edge = graph.getEdge(*prev, cur);
			if (edge != NULL) {
				// Path residual edge is a forward edge in the original graph
				edge->flow++;
				graph.backtrackStable->flowCost += edge->cost;
				if (edge->destNode < graph.backtrackStable->totalVarNodes 
					  && bestBranch != NULL) {
					(*bestBranch)[edge->destNode] = graph.backtrackStable->
																				  nodeToVal[*prev];
				}
				updateResidualGraph(*prev, cur, *edge);
			}	else {
				// Path residual edge is a backward edge in the original graph
				edge = graph.getEdge(cur, *prev);
				edge->flow--;
				if (!edge->flow) {
					graph.backtrackStable->flowCost -= edge->cost;
				}
				updateResidualGraph(cur, *prev, *edge);
			}
			*prev = cur;
		}

		// Send flow along the path. In the case of a normal path, we traverse
		// in reverse because we got it reversed from the shortest path algorithm.
		// In case of cycle, traverse normally.
		void sendFlow(int src, const vector<int>& path, BestBranch* bestBranch, 
									bool isCycle) {
			int prev = src;
			if (!isCycle) {;
				for(auto it = path.rbegin(); it != path.rend(); it++) {
					sendEdgeFlow(&prev, *it, bestBranch);
				}
			} else {
				for (unsigned int i = 0; i < path.size(); i++) {
					sendEdgeFlow(&prev, path[i], bestBranch);
				}
			}
		}

		// Iterate through the shortest path and calculate the new total flow cost,
		// without actually changing the actual total flow cost. Is useful for
		// knowing if this flow change will result in a feasible one, before 
		// applying it.
		int calculateFlowCostChange(int src, vector<int>& shortestPath) {
			int prev = src;
			int flowCost = graph.backtrackStable->flowCost;
			for (auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				// Bellman returns the path in reverse, so traverse it in reverse
				NormalEdge *e = graph.getEdge(prev, *it);
				if (e != NULL) {
					flowCost += e->cost;
				} else {
					e = graph.getEdge(*it, prev);
					flowCost -= e->cost;
				}
				prev = *it;
			}
			return flowCost;
		}

		// Repair upper bound violation by sending flow along a shortest path from
		// violation.dest to violation.src.
		// If called with non-NULL isCycle, cycles will be taken into account in
		// the search. If a cycle is found, repair the cycle instead of
		// the violation, and isCycle will be set to true.
		bool minCostFlowIteration(const EdgeInfo& violation, bool *isCycle,
														  BestBranch* bestBranch, 
															Int::IntView costUpperBound) {		
			vector<int> shortestPath;
			vector<int> dist;
			if (!findShortestPathNegativeCosts(violation.dest, violation.src, 
																				 &shortestPath, dist, isCycle)) {
				// Constraint is not consistent
				return false;
			}

			if (isCycle != NULL && *isCycle) {
				// Repair cycle
				sendFlow(shortestPath[shortestPath.size() - 1], shortestPath,
								 bestBranch, true);
				return true;		
			}

			// Make sure we won't violate cost upper bound by sending flow
			int flowCost = calculateFlowCostChange(violation.src, shortestPath);		
			if (flowCost > costUpperBound.max()) {
				return false;
			}
			sendFlow(violation.src, shortestPath, bestBranch, false);
			return true;
		}

		// Shortest path from source node to dest node
		// Return the shortest path through parameters,
		// or false as return value in case of no path.
		// If isCycle is not NULL, take cycles into account in the search.
		// If a cycle is found, isCycle will be set to true and the cycle 
		// will be placed inside path. If a cycle is not found, or if isCycle
		// is NULL, the regular shortest path will be placed instead.
		bool findShortestPathNegativeCosts(int source, int dest, 
																			 vector<int> *path, vector<int>& dist,
																			 bool *isCycle) const {
			vector<int> prev;
			bellmanFordShortestPaths(source, prev, dist, &dest, path, isCycle);
			if (isCycle != NULL && isCycle) {
				return true;
			} 

			// No path exists
			if (dist[dest] == INF_INT) {
				return false;
			}
			// Trace back shortest path
			int node = dest;
			while (node != source) {
				path->push_back(node);
				node = prev[node];
			}

			path->push_back(node);
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
		void findShortestPathsReducedCosts(int source, 
																			 unordered_set<int>& targetNodes, 
																			 vector<int>& dist, 
																			 int costLowerBound) const {
			auto& nodeList = graph.backtrackStable->nodeList;
			vector<bool> visited;
			visited.assign(nodeList.size(), false);
			dist.assign(nodeList.size(), INF_INT);
			dist[source] = 0;
			if (targetNodes.empty()) {
				return;
			}

			struct HeapItem {
				int node;
				int dist;

				HeapItem(int node, int dist) : node(node), dist(dist) 
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
				int node = curItem.node;
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

				for (auto& edge: nodeList[node].residualEdgeList) {
					if (visited[edge.destNode]) {
						continue;
					}
					int newDist = dist[node] + edge.reducedCost;
					if (newDist < dist[edge.destNode]) {
						// Found path of lower cost
						dist[edge.destNode] = newDist;
						heap.push(HeapItem(edge.destNode, newDist));
					}
				}
			}
		}

		// Optimization to skip finding shortest paths for GAC, according to 
		// Practical Improvements section of Regin's publication
		bool earlyPrune(int a, int b, int y, 
										int m) const {

			auto mFactor = [](int b, int y, FlowGraph& graph) {
				auto& nodeList = graph.backtrackStable->nodeList;
				int min = INF_INT;
				for (int z = 0 ; z < graph.backtrackStable->totalVarNodes ; z++) {
					ResidualEdge *edgeZB;
					if (z == y || ((edgeZB = graph.getResidualEdge(z, b)) == NULL)) {
						continue;
					}
					for (auto& edgeZC: nodeList[z].residualEdgeList) {
						if (edgeZC.destNode == b) {
							continue;
						}
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
				return true;
			}
			int mB = mFactor(b, y, graph);
			if (mB != INF_INT && edgeSA->flow < edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mB > m) {
				return true;
			}
			int mA = mFactor(a, y, graph);
			if (mA != INF_INT && edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow > edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA > m - edgeYB->reducedCost) {
				return true;
			}
			if (mA != INF_INT && mB != INF_INT && edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA + mB > m) {
				return true;
			}
			return false;
		}

	public:
		FlowGraphAlgorithms(FlowGraph& graph) : graph(graph) {}

		// Send flow along all lower bound violating edges, to find a feasible
		// flow.
		bool findMinCostFlow(BestBranch* bestBranch, Int::IntView costUpperBound) {
			EdgeInfo violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
				if (!minCostFlowIteration(violation, NULL, bestBranch, 
																	costUpperBound)) {
					return false;
				}
			}
			return graph.checkFlowCost(costUpperBound);
		}

		// The edges in updatedEdges are edges that have been issued for pruning, 
		// but they also contain flow, so there is an upper bound violation. 
		// Update the min cost flow to fix the violations first, and delete the 
		// edges afterwards.
		// But first, it is essential to ensure that the current flow is still of 
		// min cost.
		// When backtracking happens, some previously deleted edges are restored,
		// which means that our current flow may not necessarily be minimal still.
		// This is because the flow is not backtracked. If the flow is not minimal,
		// then there will exist at least one cycle of negative cost. So it is
		// necessary to search for such cycles first and send flow along them, 
		// to establish a min cost flow.
		// Then, repair violation and delete the edge in question.
		bool updateMinCostFlow(BestBranch* bestBranch, 
												   Int::IntView costUpperBound) {
			buildResidualGraph(bestBranch);
			bool isCycle;
			do {
				// Send flow along all cycles of negative cost, to re-establish a min 
				// cost flow
				isCycle = false;
				vector<int> prev, path;
				vector<int> dist;
				bellmanFordShortestPaths(graph.tNode(), prev, dist, NULL, &path, 
																 &isCycle);
				if (isCycle) {
					sendFlow(path[path.size() - 1], path, bestBranch, true);
				}
			} while (isCycle);

			// Repair upper bound violations
			for (auto& e: graph.updatedEdges) {
				auto res = graph.getEdge(e.src, e.dest);
				if (res == NULL) {
					// Check if edge has already been removed because it existed twice 
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
				if (res->flow) {
					// This check is needed because it is possible that this upper bound 
					// violation has already been fixed by the cycle repair algorithm 
					// above.
					if (!minCostFlowIteration({e.dest, e.src}, NULL, bestBranch, 
																		 costUpperBound)) {
						return false;
					}
				}
				
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
			graph.updatedEdges.clear();
			return graph.checkFlowCost(costUpperBound);
		}
	
		// Algorithm based on J-C. Régin, Cost-Based Arc Consistency for Global 
		// Cardinality Constraints, Constraints 7, 2002, pp. 387-405
		ExecStatus performArcConsistency(Space& home, ViewArray<Int::IntView>& vars, 
															       Int::IntView costUpperBound) {
			vector<int> distances, prev;
			bellmanFordShortestPaths(graph.tNode(), prev, distances);
		  graph.calculateReducedCosts(distances);
 
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

			// Gather the targetNodes we want to find shortests paths to from B,
			// and check early prune conditions to skip finding some
			auto& sNode = graph.backtrackStable->nodeList[graph.sNode()];
			for (int i = 0; i < graph.edgeListSize[graph.sNode()]; i++) {
				auto& edge = (sNode.edgeList.list)[i];
				if (edge.flow > 0) {
					vector<int> yList;
					int b = edge.destNode;
					unordered_set<int> targetNodes;
					vector< pair<int, int>> ayList;
					auto& bNode = graph.backtrackStable->nodeList[b];
					for (int j = 0; j < graph.edgeListSize[b]; j++) {
						auto& edgeBY = (bNode.edgeList.list)[j];
						if (edgeBY.flow == 1) {
							int y = edgeBY.destNode;
							for (IntVarValues v(vars[y]); v(); ++v) {
								int a = (graph.backtrackStable->valToNode)[v.val()];
								if (a != b) {
									if (earlyPrune(a, b, y, costUpperBound.max() - 
															   graph.backtrackStable->flowCost)) {
										edgesToPrune.push_back(EdgeWithVal(a, y, v.val()));
										continue;
									}
									targetNodes.insert(a);
									ayList.push_back({a, y});
								}
							}
						}
					}

					vector<int> reducedDistances;
					if (targetNodes.empty()) {
						continue;
					}
					findShortestPathsReducedCosts(b, targetNodes, reducedDistances, 
																			  costUpperBound.max() - 
																				graph.backtrackStable->flowCost
																			 );
					
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
							edgesToPrune.push_back(EdgeWithVal(
														 a, 
														 y, 
														 graph.backtrackStable->nodeToVal.find(a)->second));
							continue;
						}
						ResidualEdge *residualEdge = graph.getResidualEdge(a, y);
						int costAY = residualEdge->reducedCost;
						int costYB = graph.getResidualEdge(y, b)->reducedCost;
						if ((int)reducedDistances[a] > (costUpperBound.max() - 
																					  graph.backtrackStable->flowCost 
																		        - (int)costAY - (int)costYB)) {
							edgesToPrune.push_back(EdgeWithVal(
														 a, 
														 y,
														 graph.backtrackStable->nodeToVal.find(a)->second));
						} 
					}
				}
			}

			// Do the actual pruning
			for (auto& edge: edgesToPrune) {
				GECODE_ME_CHECK(vars[edge.dest].nq(home, edge.val));
			}
			return ES_OK;
		}
};

#endif