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
		bool isMinCost;

		// Add / update / delete residual edges related to the original graph edge 
		// source->dest. Updates / deletions are needed because in each iteration,
		// instead of building the residual graph from scratch, we modify the 
		// previous one only in the edges that change
		void updateResidualGraph(int source, int dest, 
														 NormalEdge edge) {
			int residualEdgeIndex;
			int residualBackwardsEdgeIndex;
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
				auto it = graph.nodeList[source].residualEdgeList->begin() + 
																				 residualEdgeIndex;
				graph.nodeList[source].residualEdgeList->erase(it);
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
				auto it = graph.nodeList[dest].residualEdgeList->begin() + 
																			 residualBackwardsEdgeIndex;
				graph.nodeList[dest].residualEdgeList->erase(it);
			}
		}

void buildResidualGraph(LI *li) {
			//cout << "Building res" << endl;
			for (int i = 0; i < graph.tNode(); i++) {
				graph.nodeList[i].residualEdgeList->clear();
			}

			for (int i = graph.totalVarNodes; i < graph.tNode(); i++) {
				auto& node = graph.nodeList[i];
				for (int j = 0; j < node.edgeList.listSize; j++) {
					auto& edge = (*node.edgeList.list)[j];
					if (edge.flow < edge.upperBound) {
						node.residualEdgeList->push_back(ResidualEdge(edge.destNode, edge.upperBound - edge.flow, edge.cost));
					}
					if (edge.flow > edge.lowerBound) {
						graph.nodeList[edge.destNode].residualEdgeList->push_back(ResidualEdge(i, edge.flow - edge.lowerBound, -edge.cost));
					}
					if (li != NULL && edge.flow && edge.destNode < graph.totalVarNodes) {
						(*li)[edge.destNode] = (*graph.nodeToVal)[i];
					}
				}
			}
		//	graph.printResidual();
		}

		

		// Bellman-Ford algorithm for shortest paths with negative costs.
		// If dest is not NULL, ignore any direct source->dest edge.
		// This is needed when searching for shortest path to a specific 
		// destination, by the min cost flow algorithm.
     void bellmanFordShortestPaths(int source, 
																	vector<int>& prev, vector<int>& dist, 
																	int* dest = NULL) const {
			prev.assign(graph.nodeList.size(), NONE_INT);
			dist.assign(graph.nodeList.size(), INF_INT);
			dist[source] = 0;
			// if (dest == NULL) {
			// 	debugCounter++;
			// 	if (debugCounter == 4569161 || debugCounter == 4569171 || debugCounter == 4569439 || debugCounter == 4569449) {
			// 		graph.print();
			// 		graph.printResidual();
			// 		return;
			// 	}
			// } //else 
				//cout << "look for " << source << "->" << *dest << endl;
/*			if (source == 21 && *dest == 5) {
				graph.print();
				graph.printResidual();
				//exit(1);
				debug = true;
			}		*/
			vector<int> updatedNodes1 = {source};
			vector<int> updatedNodes2;
			vector<int>* updatedNodesOld = &updatedNodes1;
			vector<int> *updatedNodesNew = &updatedNodes2;
			while (!updatedNodesOld->empty()) {
				bool foundUpdate = false;
				for (auto node: *updatedNodesOld) {
					for (auto &edge : (*graph.nodeList[node].residualEdgeList)) {
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
					//		cout << "dist[" << edge.destNode << "] = " << dist[node] << "+" << edge.cost << endl;
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

		void traceCycle(const vector<int>& prev, int node, vector<int>& cycle) const {
				unordered_set<int> stackContent;
				stack<int> stack;

				while (stackContent.find(node) == stackContent.end()) {
					stack.push(node);
					stackContent.insert(node);
					node = prev[node];
				}
				cycle.push_back(node);
				while (stack.top() != node) {
					cycle.push_back(stack.top());
					stack.pop();
				}
				cycle.push_back(node);
		}

		void bellmanFordShortestPathsCycles(int source, 
																	vector<int>& prev, vector<int>& dist,
																	vector<int>& cycle, bool *isCycle,
																	int* dest = NULL) const {
			// if (dest != NULL)
			// cout << "look for " << source << "->" << *dest << endl;
			// else 
			// cout << "look for " << source << "->" << "all" << endl;
				
			prev.assign(graph.nodeList.size(), NONE_INT);
			dist.assign(graph.nodeList.size(), INF_INT);
			dist[source] = 0;
			vector<int> len;
			len.assign(graph.nodeList.size(), 0);
		
			vector<int> updatedNodes1 = {source};
			vector<int> updatedNodes2;
			vector<int>* updatedNodesOld = &updatedNodes1;
			vector<int> *updatedNodesNew = &updatedNodes2;
			while (!updatedNodesOld->empty()) {
				bool foundUpdate = false;
				for (auto node: *updatedNodesOld) {
					for (auto &edge : (*graph.nodeList[node].residualEdgeList)) {
						if ((dest == NULL || 
							!(node == source && edge.destNode == *dest)) && 
							dist[node] != INF_INT && dist[node] + edge.cost < 
																			 dist[edge.destNode]) {
							// Ignore direct source->dest edge when looking for shortest path 
							// to specific requested destination
							dist[edge.destNode] = dist[node] + edge.cost;
							prev[edge.destNode] = node;
							if (++len[edge.destNode] == (int) graph.nodeList.size() - 1) {
	//							cout << "cycle lol" << endl;
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
			//cout << "done " << endl;
		}

		void sendFlow(const EdgeInfo& violation, const vector<int>& shortestPath, int minUpperBound, LI* li) {
			// Send flow through the path edges and update residual graph
			int prev = violation.src;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				NormalEdge *edge = graph.getEdge(prev, *it);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow += minUpperBound;
// if (neti)					cout << "Flow of " << prev << " " << *it << " now " << edge->flow << endl;
					*graph.flowCost += edge->cost;
					if (edge->destNode < graph.totalVarNodes && li != NULL) {
						(*li)[edge->destNode] = (*graph.nodeToVal)[prev];
					}
					updateResidualGraph(prev, *it, *edge);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(*it, prev);
					edge->flow -= minUpperBound;
	// if (neti)				cout << "Flow of " << *it << " " << prev << " now " << edge->flow << endl;
					if (!edge->flow) {
						*graph.flowCost -= edge->cost;
					}
					updateResidualGraph(*it, prev, *edge);
				}
				prev = *it;
			}
		}

		void sendFlowCycle(const vector<int>& cycle, int minUpperBound, LI* li) {
			// Send flow through the path edges and update residual graph
			int prev = cycle[0];
			for(unsigned int i = 1; i < cycle.size(); i++) {			
				NormalEdge *edge = graph.getEdge(prev, cycle[i]);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow += minUpperBound;
					*graph.flowCost += edge->cost;
					if (edge->destNode < graph.totalVarNodes && li != NULL) {
						(*li)[edge->destNode] = (*graph.nodeToVal)[prev];
					}
					updateResidualGraph(prev, cycle[i], *edge);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(cycle[i], prev);
					edge->flow -= minUpperBound;
					if (!edge->flow) {
						*graph.flowCost -= edge->cost;
					}
					updateResidualGraph(cycle[i], prev, *edge);
				}
				prev = cycle[i];
			}
		}

		int findMinUpperBoundCycle(const vector<int>& cycle) const {
			int prev = cycle[0];
			int minUpperBound = INF_INT;
			for(unsigned int i = 1; i < cycle.size(); i++) {
				// Bellman returns the path in reverse, so traverse it in reverse
				ResidualEdge *edge = graph.getResidualEdge(prev, cycle[i]);
				minUpperBound = min(minUpperBound, edge->upperBound);
				prev = cycle[i];
			}
			return minUpperBound;
		}
		bool debug = false;
		int findMinUpperBound(const EdgeInfo& violation, vector<int>& shortestPath, int* flowCost) {
			// Find min upper bound along shortest path
			int prev = violation.src;
			int minUpperBound = INF_INT;
			// cout << "Violation " << prev << " " << violation.second << endl;
			*flowCost = *graph.flowCost;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				// Bellman returns the path in reverse, so traverse it in reverse
				ResidualEdge *edge = graph.getResidualEdge(prev, *it);
/*				if (debug) {
					graph.print();
					graph.printResidual();
					if (edge == NULL) {
						cout << "cant find res edge " << prev << "->" << *it << endl;
					}
				}*/
				if (edge == NULL) {
					cout << "can't find " << prev << "->" << *it << endl;
					graph.print();
					graph.printResidual();
					assert(false); 
					exit(1);
				}
				minUpperBound = min(minUpperBound, edge->upperBound);
				NormalEdge *e = graph.getEdge(prev, *it);
				if (e != NULL) {
					*flowCost += e->cost;
				} else {
					e = graph.getEdge(*it, prev);
					*flowCost -= e->cost; // here in costsym check if flow becomes zero by all substraction
				}

				prev = *it;
			}
//			cout << endl;
			return minUpperBound;
		}

		bool minCostFlowIteration(const EdgeInfo& violation, bool *isCycle, LI* li, Int::IntView costUpperBound) {		
			vector<int> shortestPath;
			vector<int> dist;
			int pathCost; 
			//  cout << "Violation " << violation.first << "->" << violation.second << endl;
			//graph.print();
			//graph.printResidual();
			if (!findShortestPathNegativeCosts(violation.dest, violation.src, 
																				 shortestPath, dist, pathCost, isCycle)) {
				// Constraint is not consistent
				return false;
			}

			if (isCycle != NULL && *isCycle) {
				int minUpperBound = findMinUpperBoundCycle(shortestPath);
				sendFlowCycle(shortestPath, minUpperBound, li);
/*					graph.dist->clear();
			for (auto d: dist) {
				graph.dist->push_back(d);
			}*/
				return true;		
			}

		//	updatePotentials(visitedNodes, dist, pathCost);
		//	updateCosts();

			int flowCost = 0;
			int minUpperBound = findMinUpperBound(violation, shortestPath, &flowCost);		
			if (flowCost > costUpperBound.max()) {
				return false;
			}
			sendFlow(violation, shortestPath, minUpperBound, li);
			
			// TODO: optimize
/*			graph.dist->clear();
			for (auto d: dist) {
				graph.dist->push_back(d);
			}*/

			return true;
		}

		// Shortest path from source node to dest node
		// Return the shortest path and cost through parameters,
		// or false as return value in case of no path
		bool findShortestPathNegativeCosts(int source, int dest, 
																			 vector<int>& path, vector<int>& dist,
																			int& cost, bool *isCycle) 
																			 const {
			vector<int> prev;
			
			if (isCycle == NULL) {
				bellmanFordShortestPaths(source, prev, dist, &dest);
			} else {
				bellmanFordShortestPathsCycles(source, prev, dist, path, isCycle, &dest);
				if (path.size()) {
					*isCycle = true;
/*					for (auto node: path) {
						cout << node << "->";
					}
					cout << endl;*/
					return true;
				} 
			}


			// No path exists
			if (dist[dest] == INF_INT) {
				return false;
			}
			// Trace back shortest path
			int node = dest;
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
		void findShortestPathsReducedCosts(int source, 
																			 unordered_set<int>& targetNodes, 
																			 vector<int>& dist, 
																			 int costLowerBound) const {
			vector<bool> visited;
			visited.assign(graph.nodeList.size(), false);
			dist.assign(graph.nodeList.size(), INF_INT);
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

				for (auto& edge: *graph.nodeList[node].residualEdgeList) {
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
		// Practical Improvements section of the research paper
		bool earlyPrune(int a, int b, int y, 
										int m) const {

			auto mFactor = [](int b, int y, FlowGraph& graph) {
				int min = INF_INT;
				for (int z = 0 ; z < graph.totalVarNodes ; z++) {
					ResidualEdge *edgeZB;
					if (z == y || ((edgeZB = graph.getResidualEdge(z, b)) == NULL)) {
						continue;
					}
					for (auto& edgeZC: *graph.nodeList[z].residualEdgeList) {
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
			return false; 
			int mB = mFactor(b, y, graph);
			if (mB != INF_INT && edgeSA->flow < edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mB > m) {
				//if (a == 10 && !y) return false;
				cout << "\tCONDITION 2 EARLY PRUNNING " << a << " (val "  << (*graph.nodeToVal)[a] << ") from " << y << endl;
				cout << edgeSA->flow << " " << edgeSA->upperBound << " " << edgeSB->flow << " " << 
				edgeSB->lowerBound << " " << edgeAY->reducedCost << " " << mB << " " << m << endl; 
				//graph.print();
				//graph.printResidual();
				//if (a == 12 && !y) exit(1);
				return true;
			}
			int mA = mFactor(a, y, graph);
			if (mA != INF_INT && edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow > edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA > m - edgeYB->reducedCost) {
		//		cout << "\tCONDITION 3 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			if (mA != INF_INT && mB != INF_INT && edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA + mB > m) {
		//		cout << "\tCONDITION 4 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			return false;
		}

	public:
		FlowGraphAlgorithms(FlowGraph& graph) : graph(graph) {}

		bool findMinCostFlow(LI* li, Int::IntView costUpperBound) {
			EdgeInfo violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
//				cout << "Violation " << violation.first << "->" << violation.second << endl;
				if (!minCostFlowIteration(violation, NULL, li, costUpperBound)) {
				//	cout << "incosistent" << endl;
					return false;
				}
			}
			isMinCost = true;
			//graph.calculateFlowCost(li);
			// graph.print();
			return graph.checkFlowCost(costUpperBound);
		}

		// Given updatedEdges contains the edges whose bounds have been tightened
		// since last execution, do the following:
		// - Update the residual graph to match the changes
		// - If the old flow is not still feasible, find a new one, using the 
		//   incremental algorithm from the publication
		bool updateMinCostFlow(vector<EdgeInfo>& updatedEdges, LI* li, Int::IntView costUpperBound) {
			  //  cout << "Propagate: update min cost flow" << endl;
		//	graph.print();
			buildResidualGraph(li);
			// if (*graph.oldFlowIsFeasible) {
			// 	isMinCost = false;
			// 	cout << "is already feasible" << endl;
			// 	return true;
			// }
			assert(updatedEdges.size());
			if (!updatedEdges.size()) {
				cout << "problem" << endl;
				graph.print();
				exit(1);
			}

			bool isCycle;
			do {
					isCycle = false;
					vector<int> prev, path;
					vector<int> dist;
					bellmanFordShortestPathsCycles(graph.tNode(), prev, dist, path, &isCycle, NULL);
					if (path.size()) {
						isCycle = true;
/*						for (auto node: path) {
							cout << node << "->";
						}
						cout << endl;*/
						int minUpperBound = findMinUpperBoundCycle(path);
						sendFlowCycle(path, minUpperBound, li);
						//assert(false);
						//exit(1);
					}
				} while (isCycle);

			// graph.print();
			for (auto& e: updatedEdges) {
				// Among the edges that changed, look for one violating bounds
				// Assume violating lower bound initially
				  // cout << e.src << "->" << e.dest << endl;
				//  if (e.src == 23 && e.dest == 11) {
				// 	 neti = true;
				//  } else {
				// 	 neti = false;
				//  }
				auto res = graph.getEdge(e.src, e.dest);
				if (res == NULL) {
					//  cout << "double entry" << endl;
					continue;
				}

				if (res->flow) {

				int src = e.src;
				int dest = e.dest;
				// Violating upper bound, swap direction of initial violating edge
					std::swap(src, dest);
					//src = e.dest;
					//dest = graph.tNode();
				
				bool isCycle = false;
				// graph.printResidual();
				if (!minCostFlowIteration({src, dest}, NULL, li, costUpperBound)) {
					return false;
				}
				// graph.printResidual();
				if (isCycle) {
					assert(false);
					cout << "i didn't expect a cycle here" << endl;
					exit(1);
					while (graph.getEdge(e.src, e.dest)->flow) {
						isCycle = false;
						debug = true;
						if (!minCostFlowIteration({src, dest}, &isCycle, li, costUpperBound)) {
							return false;
						}
					}
				}

				}
				
				graph.deleteEdge(e.src, e.dest);
				// We have already removed flow from this Val->Var edge and updated 
				// residual graph, so we know for sure that the only residual edge 
				// direction will be the same as the original edge.
				graph.deleteResidualEdge(e.src, e.dest);
				int val = (*graph.nodeToVal)[e.src];
				graph.varToVals[e.dest].deleteVal(val);
			
				while (isCycle) {
					isCycle = false;
					vector<int> prev, path;
					vector<int> dist;
					bellmanFordShortestPathsCycles(graph.sNode(), prev, dist, path, &isCycle, NULL);
					if (path.size()) {
						isCycle = true;
/*						for (auto node: path) {
							cout << node << "->";
						}
						cout << endl;*/
						int minUpperBound = findMinUpperBoundCycle(path);
						sendFlowCycle(path, minUpperBound, li);
						//assert(false);
						//exit(1);
					}					
					
					// TODO: alg that just send flow to cycle without violation at the same time
				}
			
			}
			if (!updatedEdges.size()) {
				cout << "this shouldn't happen" << endl;
				assert(false);
				exit(1);	
				bool isCycle = false;
				do {
					isCycle = false;
					vector<int> prev, path;
					vector<int> dist;
					bellmanFordShortestPathsCycles(graph.sNode(), prev, dist, path, &isCycle, NULL);
					if (path.size()) {
						isCycle = true;
/*						for (auto node: path) {
							cout << node << "->";
						}
						cout << endl;*/
						int minUpperBound = findMinUpperBoundCycle(path);
						sendFlowCycle(path, minUpperBound, li);
						//assert(false);
						//exit(1);
					}
				}	while (isCycle);				
			}
			isMinCost = true;
			return graph.checkFlowCost(costUpperBound);
		}
	
		// In addition to pruning, hold the affected Val->Var edges in updatedEdges
		// so we can update the residual graph later accordingly. We do not update
		// it here, because in case the home space fails due to another
		// constraint, or if we find a solution from this pruning, we would have 
		// updated it for no reason, as we wouldn't need to re-check the validity of 
		// costgcc, the search would backtrack to previous instances.
		// The reason why 

		ExecStatus performArcConsistency(Space& home, ViewArray<Int::IntView>& vars, 
															       LI* li, Int::IntView costUpperBound) {
		//	graph.addTResidualEdges(); // opt?
			vector<int> distances;
			vector<int> prev;
		//	cout << "in arc" << endl;
		isMinCost = true;
			if (isMinCost) {
		// 		using std::chrono::high_resolution_clock;
    // using std::chrono::duration_cast;
    // using std::chrono::duration;
    // using std::chrono::milliseconds;

    // auto t1 = high_resolution_clock::now();
				bellmanFordShortestPaths(graph.tNode(), prev, distances, NULL);
			// auto t2 = high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    // auto ms_int = duration_cast<milliseconds>(t2 - t1);

    // /* Getting number of milliseconds as a double. */
    // duration<double, std::milli> ms_double = t2 - t1;

    // std::cout << ms_int.count() << "ms\n";
			} else {
				bool isCycle = false;
				do {
					isCycle = false;
					// cout << "after cycle check" << endl;
					vector<int> path;
					bellmanFordShortestPathsCycles(graph.tNode(), prev, distances, path, &isCycle, NULL);
					// cout << "after cycle done" << endl;
					if (path.size()) {
						// cout << "after cycle FOUND ARC" << endl;
						isCycle = true;
						// for (auto node: path) {
							// cout << node << "->";
						// }
						// cout << endl;
						int minUpperBound = findMinUpperBoundCycle(path);
						sendFlowCycle(path, minUpperBound, li);
						
						//assert(false);
						//exit(1);
					}
				} while (isCycle);
			}
			
		  graph.calculateReducedCosts(distances);
			//graph.printResidual();
 
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
			auto& sNode = graph.nodeList[graph.sNode()];
			for (int i = 0; i < sNode.edgeList.listSize; i++) {
				auto& edge = (*sNode.edgeList.list)[i];
				if (edge.flow > 0) {
					vector<int> yList;
					int b = edge.destNode;
					unordered_set<int> targetNodes;
					vector< pair<int, int>> ayList;
					auto& bNode = graph.nodeList[b];
					for (int j = 0; j < bNode.edgeList.listSize; j++) {
						auto& edgeBY = (*bNode.edgeList.list)[j];
						if (edgeBY.flow == 1) {
							int y = edgeBY.destNode;
						//	minDist[y].bestVal = b;
							for (IntVarValues v(vars[y]); v(); ++v) {
								int a = (*graph.valToNode)[v.val()];
								if (a != b) {
									if (earlyPrune(a, b, y, costUpperBound.max() - 
															   *(graph.flowCost))) {
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
																				costUpperBound.max() - *(graph.flowCost));
					
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
						int costAY = residualEdge->reducedCost;
						int costYB = graph.getResidualEdge(y, b)->reducedCost;
						if ((int)reducedDistances[a] > (costUpperBound.max() - *(graph.flowCost) 
																		  - (int)costAY - (int)costYB)) {
							edgesToPrune.push_back(EdgeWithVal(a, y,
																							graph.nodeToVal->find(a)->second));
						} 
					}
				}
			}


			// Do the actual pruning and update data structures
			for (auto& edge: edgesToPrune) {
				// Push to updatedEdges so we can modify the residual graph accordingly
				// on the next min cost flow computation
	//			updatedEdges.push_back(EdgeInfo(edge.src, edge.dest, false, false, true));
				// Prune
				GECODE_ME_CHECK(vars[edge.dest].nq(home, edge.val));
	//			cout << "Prunning " << edge.src << " " << edge.dest << endl;
				// Also remove from varToVals
	//			auto& vals = graph.varToVals[edge.dest];
	//			vals.deleteVal(edge.val);
				// Update upper bound
	//			assert(!actualEdge->flow);
	//			graph.deleteEdge(edge.src, edge.dest);
				/*if (vars[edge.dest].assigned()) {
					// If a variable got assigned by pruning, set corresponding edge
					// lower bound to 1
					int assignedVal = vars[edge.dest].val();
					assert(*vals.begin() == assignedVal);
					auto valNode = graph.valToNode->find(assignedVal)->second;
					graph.getEdge(valNode, edge.dest)->lowerBound = 1;
				}*/
			}
			// cout << "done prunning" << endl;

   // 	graph.removeTResidualEdges();
			return ES_OK;
		}

ExecStatus performArcConsistencyBell(Space& home, ViewArray<Int::IntView>& vars, 
															       Int::IntView costUpperBound) {
		//	graph.addTResidualEdges(); // opt?
			vector<int> distances;
			vector<int> prev;
   // using std::chrono::milliseconds;

    // auto t1 = high_resolution_clock::now();
				bellmanFordShortestPaths(graph.tNode(), prev, distances, NULL);
 
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
			auto& sNode = graph.nodeList[graph.sNode()];
			for (int i = 0; i < sNode.edgeList.listSize; i++) {
				auto& edge = (*sNode.edgeList.list)[i];
				if (edge.flow > 0) {
					vector<int> yList;
					int b = edge.destNode;
					unordered_set<int> targetNodes;
					vector< pair<int, int>> ayList;
					auto& bNode = graph.nodeList[b];
					for (int j = 0; j < bNode.edgeList.listSize; j++) {
						auto& edgeBY = (*bNode.edgeList.list)[j];
						if (edgeBY.flow == 1) {
							int y = edgeBY.destNode;
						//	minDist[y].bestVal = b;
							for (IntVarValues v(vars[y]); v(); ++v) {
								int a = (*graph.valToNode)[v.val()];
								if (a != b) {
									targetNodes.insert(a);
									ayList.push_back({a, y});
								}
							}
						}
					}

					if (targetNodes.empty()) {
						continue;
					}
					prev.clear();
					distances.clear();
					bellmanFordShortestPaths(b, prev, distances, NULL);;
					
					for (const auto& ay: ayList) {
						const auto a = ay.first;
						const auto y = ay.second;
						ResidualEdge *residualEdge = graph.getResidualEdge(a, y);
						int costAY = residualEdge->cost;
						int costYB = graph.getResidualEdge(y, b)->cost;
						if ((int)distances[a] > (costUpperBound.max() - *(graph.flowCost) 
																		  - (int)costAY - (int)costYB)) {
							edgesToPrune.push_back(EdgeWithVal(a, y,
																							graph.nodeToVal->find(a)->second));
						} 
					}
				}
			}


			// Do the actual pruning and update data structures
			for (auto& edge: edgesToPrune) {
				// Push to updatedEdges so we can modify the residual graph accordingly
				// on the next min cost flow computation
	//			updatedEdges.push_back(EdgeInfo(edge.src, edge.dest, false, false, true));
				// Prune
				GECODE_ME_CHECK(vars[edge.dest].nq(home, edge.val));
	//			cout << "Prunning " << edge.src << " " << edge.dest << endl;
				// Also remove from varToVals
	//			auto& vals = graph.varToVals[edge.dest];
	//			vals.deleteVal(edge.val);
				// Update upper bound
	//			assert(!actualEdge->flow);
	//			graph.deleteEdge(edge.src, edge.dest);
				/*if (vars[edge.dest].assigned()) {
					// If a variable got assigned by pruning, set corresponding edge
					// lower bound to 1
					int assignedVal = vars[edge.dest].val();
					assert(*vals.begin() == assignedVal);
					auto valNode = graph.valToNode->find(assignedVal)->second;
					graph.getEdge(valNode, edge.dest)->lowerBound = 1;
				}*/
			}
			// cout << "done prunning" << endl;

   // 	graph.removeTResidualEdges();
			return ES_OK;
		}

};

#endif