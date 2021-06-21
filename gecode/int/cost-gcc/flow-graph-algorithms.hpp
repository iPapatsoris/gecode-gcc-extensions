#ifndef FLOW-GRAPH-ALGORITHMS
#define FLOW-GRAPH-ALGORITHMS 

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
unsigned int countt = 0;

class FlowGraphAlgorithms {
	private:
		FlowGraph& graph;

		// Add / update / delete residual edges related to the original graph edge 
		// source->dest. Updates / deletions are needed because in each iteration,
		// instead of building the residual graph from scratch, we modify the 
		// previous one only in the edges that change
		void updateResidualGraph(unsigned int source, unsigned int dest, 
														 NormalEdge edge, bool onPath) {
			unsigned int residualEdgeIndex;
			unsigned int residualBackwardsEdgeIndex;
			ResidualEdge* residualEdgeSearch = graph.getResidualEdge(source, dest, 
																														&residualEdgeIndex);
			ResidualEdge* residualBackwardsEdgeSearch = graph.getResidualEdge(dest, source,
																								 &residualBackwardsEdgeIndex);
			unsigned int reducedCost = 0;
			if (edge.flow < edge.upperBound) {
				// Add / update forward residual edge
				if (!onPath) {		
					reducedCost = (residualEdgeSearch != NULL ? residualEdgeSearch->reducedCost : residualBackwardsEdgeSearch->reducedCost);	
				}
				graph.setOrCreateResidualEdge(residualEdgeSearch, source, 
																		  ResidualEdge(dest, 
																								   edge.upperBound - edge.flow, 
																									 edge.cost, reducedCost, onPath));
			//	cout << "Adding res edge " << source << "->" << dest << endl;
			} else if (residualEdgeSearch != NULL) {
				// Delete forward residual edge that should no longer exist
				auto it = graph.nodeList[source].residualEdgeList.begin() + 
																				 residualEdgeIndex;
				graph.nodeList[source].residualEdgeList.erase(it);
			//	graph.orderGraph.removeEdge(source, dest);
			//	cout << "Deleting res edge " << source << "->" << dest << endl;
			}

			reducedCost = 0;
			if (edge.flow > edge.lowerBound) {
				// Add / update backward residual edge
				if (!onPath) {
					reducedCost = (residualBackwardsEdgeSearch != NULL ? residualBackwardsEdgeSearch->reducedCost : residualEdgeSearch->reducedCost);
				}
				graph.setOrCreateResidualEdge(residualBackwardsEdgeSearch, dest, 
																		  ResidualEdge(source, 
																									 edge.flow - edge.lowerBound, 
																									 edge.cost, reducedCost, onPath));
				//graph.orderGraph.addEdge(dest, source);
			//	cout << "Adding res edge " << dest << "->" << source << endl;
			} else {
				/*if (edge.flow == edge.upperBound && edge.flow == edge.lowerBound) {
					graph.completedResidualEdges.push_back({dest, ResidualEdge(source, 0, edge.cost, 0, true)});
				} else*/ if (residualBackwardsEdgeSearch != NULL) {
				// Delete backward residual edge that should no longer exist
				auto it = graph.nodeList[dest].residualEdgeList.begin() + 
																			 residualBackwardsEdgeIndex;
				graph.nodeList[dest].residualEdgeList.erase(it);
				//graph.orderGraph.removeEdge(dest, source);
			//	cout << "Deleting res edge " << dest << "->" << source << endl;
				}
			//graph.orderGraph.print();
			}
	}

		// Bellman-Ford algorithm for shortest paths with negative costs.
		// If dest is not NULL, ignore any direct source->dest edge.
		// This is needed when searching for shortest path to a specific 
		// destination, by the min cost flow algorithm.
		// TODO: consider testing the randomized variation improvement
     void bellmanFordShortestPaths(unsigned int source, 
																	vector<unsigned int>& prev, vector<unsigned int>& dist, 
																	unsigned int* dest = NULL) const {
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist.assign(graph.nodeList.size(), INF_UINT);
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
							dist[node] != INF_UINT && dist[node] + edge.cost < 
																			 dist[edge.destNode]) {
							// Ignore direct source->dest edge when looking for shortest path 
							// to specific requested destination
							dist[edge.destNode] = dist[node] + edge.cost;
							prev[edge.destNode] = node;
							foundUpdate = true;
							updatedNodesNew->push_back(edge.destNode);
							countt++;
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
			//cout << "Total bellman: " << countt << "\n";
		}

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
		}

			void sendFlow(pair<unsigned int, unsigned int>& violation, vector<unsigned int>& shortestPath, unsigned int minUpperBound) {
				// Send flow through the path edges and update residual graph
			unsigned int prev = violation.first;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				NormalEdge *edge = graph.getEdge(prev, *it);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow += minUpperBound;
					updateResidualGraph(prev, *it, *edge, prev != violation.first);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(*it, prev);
					edge->flow -= minUpperBound;
					updateResidualGraph(*it, prev, *edge, prev != violation.first);
				}
				prev = *it;
			}
		}

		unsigned int findMinUpperBound(pair<unsigned int, unsigned int>& violation, vector<unsigned int>& shortestPath) {
			// Find min upper bound along shortest path
			unsigned int prev = violation.first;
			unsigned int minUpperBound = INF_UINT;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				// Bellman returns the path in reverse, so traverse it in reverse
				ResidualEdge *edge = graph.getResidualEdge(prev, *it);
				minUpperBound = min(minUpperBound, edge->upperBound);
				prev = *it;
				/*if (violation.first == 49 && violation.second == 7) {
				cout << *it << (it != shortestPath.rend()-1 ? "->" : "\n");
				}*/
			}
			return minUpperBound;
		}

		bool minCostFlowIteration(pair<unsigned int, unsigned int> violation) {
		//	cout << "Violation " << violation.first << "->" << violation.second << "\n";
		/*if (violation.first == 49 && violation.second == 7) 
		 { 
				cout << "Violation " << violation.first << "->" << violation.second << "\n";
				graph.print();
				graph.printResidual();
			}*/
			vector<unsigned int> shortestPath, dist, visitedNodes;
			int pathCost; 
			if (!findShortestPathNegativeCosts(violation.second, violation.first, 
																				 shortestPath, dist, visitedNodes, pathCost)) {
				// Constraint is not consistent
				return false;
			}


			updatePotentials(visitedNodes, dist, pathCost);

			/*for (unsigned int n = 0; n < visitedNodes.size(); n++) {
				if (visitedNodes[n]) {
					assert(dist[n] != INF_UINT);
					graph.nodeList[n].potential += (pathCost - dist[n]);
				}
			}*/

			updateCosts();

			/*for (unsigned int i = 0; i < graph.nodeList.size(); i++) { // can opt to not iterate all
				for (auto& e: graph.nodeList[i].residualEdgeList) {
												// SIGOURA cost edw kai OXI reduced. alliws sto paper bgianoun arnitika
					e.reducedCost = (e.isBackwards ? 0 : e.cost - graph.nodeList[i].potential + graph.nodeList[e.destNode].potential);
				}
			}*/

			unsigned int minUpperBound = findMinUpperBound(violation, shortestPath);

			// Find min upper bound along shortest path
			/*unsigned int prev = violation.first;
			unsigned int minUpperBound = INF_UINT;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				// Bellman returns the path in reverse, so traverse it in reverse
				ResidualEdge *edge = graph.getResidualEdge(prev, *it);
				minUpperBound = min(minUpperBound, edge->upperBound);
				prev = *it;
			}*/
		//	cout << "Graph after reduced costs\n";
		//	graph.printResidual();
			//assert(minUpperBound != INF_UINT);

			sendFlow(violation, shortestPath, minUpperBound);

			// Send flow through the path edges and update residual graph
			/*prev = violation.first;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				NormalEdge *edge = graph.getEdge(prev, *it);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow += minUpperBound;
					updateResidualGraph(prev, *it, *edge, prev != violation.first);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(*it, prev);
					edge->flow -= minUpperBound;
					updateResidualGraph(*it, prev, *edge, prev != violation.first);
				}
				prev = *it;
			}*/
//			cout << "New residual graph\n";
	//		graph.printResidual();
		/*	if (violation.first == 49 && violation.second == 7) {
				graph.print();
			}*/ 
			return true;
		}

		// Shortest path from source node to dest node
		// Return the shortest path and cost through parameters,
		// or false as return value in case of no path
		bool findShortestPathNegativeCosts(unsigned int source, unsigned int dest, 
																			 vector<unsigned int>& path, vector<unsigned int>& dist,
																			vector<unsigned int>& visited, int& cost) 
																			 const {
			vector<unsigned int> prev;
			
			// bellmanFordShortestPaths(source, prev, dist, &dest);
			disjjtra2(source, dest, prev, dist, visited);

			/*cout << "SPs from " << source << " to " << dest << endl;
			for (auto n: dist) {
				cout << n << "\n";
			}
			cout << "\n";
			for (auto n: prev) {
				cout << n << endl;
			}*/

			// No path exists
			if (dist[dest] == INF_UINT) {
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
		void disjjtra2(unsigned int source, unsigned int dest, vector<unsigned int>& prev, 
																			 vector<unsigned int>& dist, vector<unsigned int>& visited
																			 ) const {
			dist.assign(graph.nodeList.size(), INF_UINT);
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist[source] = 0;
			visited.assign(graph.nodeList.size(), false);

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

				if (dist[node] < curItem.dist) {
					// We have already found a better path than the one this heap item 
					// suggests
					continue;
				}

				for (auto& edge: graph.nodeList[node].residualEdgeList) {
					if (visited[edge.destNode] || (node == source && edge.destNode == dest)) {
						continue;
					}
					unsigned int newDist = dist[node] + edge.reducedCost; //graph.getReducedCost(edge, node, edge.destNode);
					if (newDist < dist[edge.destNode]) {
						// Found path of lower cost
						dist[edge.destNode] = newDist;
						prev[edge.destNode] = node;
						heap.push(HeapItem(edge.destNode, newDist));
						countt++;
					}
				}
				if (node == dest) {
					break;
				}
			}
		//	cout << "Total disktra: " << countt << "\n";
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
				cout << "\tCONDITION 1 EARLY PRUNNING " << a << " from " << y << endl; 
				return true;
			} 
			unsigned int mB = mFactor(b, y, graph);
			if (edgeSA->flow < edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mB > m) {
				cout << "\tCONDITION 2 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			unsigned int mA = mFactor(a, y, graph);
			if (edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow > edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA > m - edgeYB->reducedCost) {
				cout << "\tCONDITION 3 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			if (edgeSA->flow == edgeSA->upperBound 
					&& edgeSB->flow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA + mB > m) {
				cout << "\tCONDITION 4 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			return false;
		}

	public:
		FlowGraphAlgorithms(FlowGraph& graph) : graph(graph) {}
		
		void benchmark() const {
			vector<unsigned int> dist;
			vector<unsigned int> visitedNodes;
			vector<unsigned int> prev;
			for (unsigned int it = 0; it < 5000; it++) {
				for (unsigned int i = 0; i < graph.nodeList.size(); i++) {
					for (unsigned int j = 0; j < graph.nodeList.size(); j++) {
						disjjtra2(i, j, prev, dist, visitedNodes);
						/*for (unsigned int n = 0; n < visitedNodes.size(); n++) {
							if (visitedNodes[n]) {
								assert(dist[n] != INF_UINT);
								graph.nodeList[n].potential += (dist[j] - dist[n]);
							}
						}
						bellmanFordShortestPaths(i, prev, dist, &j);*/
					}
				}
			}

			cout << "Iterations: " << countt << "\n";
			return;
		}

		bool findMinCostFlow() {
			pair<unsigned int, unsigned int> violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
				if (!minCostFlowIteration(violation)) {
					return false;
				}
			}
			graph.oldFlowIsFeasible = true;
			graph.calculateFlowCost();
			return graph.checkFlowCost();
		}

		// Given updatedEdges contains the edges whose bounds have been tightened
		// since last execution, do the following:
		// - Update the residual graph to match the changes
		// - If the old flow is not still feasible, find a new one, using the 
		//   incremental algorithm from the publication
		bool updateMinCostFlow(vector<EdgeNodes>& updatedEdges) {
			//if (graph.printDebug) {
		//		cout << "Advisor on this (updated) graph\n";
		//		graph.print();
		//		graph.printResidual();
		//	}
			vector<NormalEdge*> edgeReference;
			for (auto& edge: updatedEdges) {
			//	cout << "Advising " << edge.first << "->" << edge.second.destNode << " feasible " << graph.oldFlowIsFeasible <<  endl;
				edgeReference.push_back(graph.getEdge(edge.first, edge.second));
				updateResidualGraph(edge.first, edge.second, *edgeReference.back(), false);
			}
			/*if (graph.printDebug) {
				cout << "Updated residual\n";
				graph.printResidual();
			}*/
			if (graph.oldFlowIsFeasible) {
				return true;
			}
			for (unsigned int i = 0 ; i < updatedEdges.size(); i++) {
				// Among the edges that changed, look for one violating bounds
				// Assume we are violating lower bound on init
				auto e = edgeReference[i];
				unsigned int src = updatedEdges[i].first;
				unsigned int dest = updatedEdges[i].second;
				/*	if (src == 49 && dest == 7) {
					cout << "UPDATEZ\n";
					cout << e.flow << " " << e.lowerBound << " " << e.upperBound << endl;
				}*/ 
				if (e->flow >= e->lowerBound && e->flow <= e->upperBound) {
					// All bounds satisfied, try another edge
					continue;
				}
				if (e->flow > e->upperBound) {
					// Violating upper bound, swap direction of initial violating edge
					std::swap(src, dest);
				}
				if (!minCostFlowIteration({src, dest})) {
				//	cout << "rip violation " << edge.first << " " << e->destNode << endl;
					return false;
				}
			}
			/*if (!foundFeasibleFlow) {
				return false;
			}*/
			graph.oldFlowIsFeasible = true;
			graph.calculateFlowCost();
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
			/*	cout << "Arc consistency on this graph\n";
				graph.print();
				graph.printResidual();
			*/
//			graph.addTResidualEdges();
//			cout << "after addTResidualEdges\n";
		//	graph.printResidual();
	//		vector<unsigned int> distances;
	//		cout << "distances\n";
		//	findShortestPathsNegativeCosts(graph.tNode(), distances);
		//	for (auto d: distances) {
			//	cout << d << "\n";
	//		}
		//	graph.calculateReducedCosts(distances);
  	//	graph.printResidual();

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

			// Gather the targetNodes we want to find shortests paths to from B,
			// and check early prune conditions to skip finding some
			for (auto& edge: graph.nodeList[graph.sNode()].edgeList) {
				if (edge.flow > 0) {
					unsigned int b = edge.destNode;
					unordered_set<unsigned int> targetNodes;
					vector< pair<unsigned int, unsigned int>> ayList;
					for (auto& edgeBY: graph.nodeList[b].edgeList) {
						if (edgeBY.flow == 1) {
							unsigned int y = edgeBY.destNode;
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
							cout << "\tEARLY PRUNING " << a << " FROM " << y << endl;
							edgesToPrune.push_back(EdgeWithVal(a, y, 
																						  graph.nodeToVal->find(a)->second));
							continue;
						}
						ResidualEdge *residualEdge = graph.getResidualEdge(a, y);
						unsigned int costAY = residualEdge->reducedCost;
						unsigned int costYB = graph.getResidualEdge(y, b)->reducedCost;
					//	cout << reducedDistances[a] << "\n" << graph.costUpperBound << "\n" << graph.flowCost << "\n" << costAY << "\n" << costYB << "\n" << endl;
						//cout << graph.costUpperBound - graph.flowCost 
						//												  - (int)costAY - (int)costYB << endl;
						if ((int)reducedDistances[a] > (graph.costUpperBound - graph.flowCost 
																		  - (int)costAY - (int)costYB)) {
					/*		cout << "\tPRUNING " << a << " FROM " << y << endl;
							if (a == 49 && y == 7) {
								graph.print();
								graph.printResidual();
								for (auto d: reducedDistances) {
									cout << d << "\n";
								}
								//exit(1);
							}*/
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

    //	graph.removeTResidualEdges();
			return ES_OK;
		}

};

#endif