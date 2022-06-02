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
#define MINUS_INF_INT INT_MIN
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
			unsigned int flow = graph.getEdgeFlow(source, dest);
			if (flow < edge.upperBound) {
				// Add / update forward residual edge
				graph.setOrCreateResidualEdge(residualEdgeSearch, source, 
																		  ResidualEdge(dest, 
																								   edge.upperBound - flow, 
																									 edge.cost, graph.getReducedCost(source, dest, edge.cost)));
			} else if (residualEdgeSearch != NULL) {
				// Delete forward residual edge that should no longer exist
				auto it = graph.nodeList[source].residualEdgeList.begin() + 
																				 residualEdgeIndex;
				graph.nodeList[source].residualEdgeList.erase(it);
			//	graph.orderGraph.removeEdge(source, dest);
			}

			if (flow > edge.lowerBound) {
				// Add / update backward residual edge
				graph.setOrCreateResidualEdge(residualBackwardsEdgeSearch, dest, 
																		  ResidualEdge(source, 
																									 flow - edge.lowerBound, 
																									 -edge.cost, graph.getReducedCost(dest, source, -edge.cost)));
				//graph.orderGraph.addEdge(dest, source);
			} else if (residualBackwardsEdgeSearch != NULL) {
				// Delete backward residual edge that should no longer exist
				auto it = graph.nodeList[dest].residualEdgeList.begin() + 
																			 residualBackwardsEdgeIndex;
				graph.nodeList[dest].residualEdgeList.erase(it);
			}
		}

		void sendFlow(pair<unsigned int, unsigned int>& violation, vector<unsigned int>& shortestPath, unsigned int minUpperBound, LI* li) {
			struct EdgeInfo {
				unsigned int src;
				NormalEdge *edge;
				bool isForward;

				EdgeInfo(unsigned int src, NormalEdge *edge, bool isForward) :
					src(src), edge(edge), isForward(isForward) {}
			};
			
			// Send flow through the path edges and update residual graph
			unsigned int prev = violation.first;
			vector<EdgeInfo> todo;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				NormalEdge *edge = graph.getEdge(prev, *it);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					if (edge->destNode < graph.totalVarNodes) {
						graph.flow.addValVarFlow(prev, *it);
						graph.flowCost += edge->cost;
						if (edge->destNode < graph.totalVarNodes && li != NULL) {
							(*li)[edge->destNode] = (*graph.nodeToVal)[prev];
						}
					updateResidualGraph(prev, *it, *edge);	
					} else {
						todo.push_back(EdgeInfo(prev, edge, true));
					}
				
					//edge->flow += minUpperBound;
//					cout << "Flow of " << prev << " " << *it << " now " << graph.getEdgeFlow(prev, *it) << endl;
				
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(*it, prev);
					if (edge->destNode < graph.totalVarNodes) {
						graph.flow.removeValVarFlow(*it, prev);
							graph.flowCost -= edge->cost;
					updateResidualGraph(*it, prev, *edge);
					} else {
						todo.push_back(EdgeInfo(*it, edge, false));
					}
					//edge->flow -= minUpperBound;
//					cout << "Flow of " << *it << " " << prev << " now " << graph.getEdgeFlow(*it, prev) << endl;
				
				}
				prev = *it;
			}

			for (auto p: todo) {
				if (graph.flow.isFirstFlow() && p.edge->destNode == graph.tNode()) {
					if (p.isForward)
						graph.flow.addVarTFlow(p.src);
					else 
						graph.flow.removeVarTFlow(p.src);
			}
				updateResidualGraph(p.src, p.edge->destNode, *p.edge);

		}

		}

		unsigned int findMinUpperBound(pair<unsigned int, unsigned int>& violation, vector<unsigned int>& shortestPath, int* flowCost) {
			// Find min upper bound along shortest path
			unsigned int prev = violation.first;
			unsigned int minUpperBound = INF_UINT;
			// cout << "SP: ";
			// for (auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
			// 	cout << *it << "->";
			// }
			*flowCost = graph.flowCost;
			for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				// Bellman returns the path in reverse, so traverse it in reverse
				ResidualEdge *edge = graph.getResidualEdge(prev, *it);
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
			// cout << endl;
			return minUpperBound;
		}
		unsigned int debugCount = 0;


		bool minCostFlowIteration(pair<unsigned int, unsigned int> violation, LI* li) {		
			vector<unsigned int> shortestPath, dist;
			// cout << debugCount++ << endl; 
			// cout << "Violation " << violation.first << "->" << violation.second << endl;
/*			if (violation.first == 19 && violation.second == 23) {
				cout << "neti" << endl;
				cout << graph.flowCost << endl;
				graph.print();
			}*/
			//graph.print();
			//graph.printResidual();
			if (!findShortestPathNegativeCosts(violation.second, violation.first, 
																				 shortestPath, dist)) {
				// Constraint is not consistent
				return false;
			}

		//	updatePotentials(visitedNodes, dist, pathCost);
		//	updateCosts();

			int flowCost = 0;
			unsigned int minUpperBound = findMinUpperBound(violation, shortestPath, &flowCost);		
			//assert(flowCost == 11);
				// TODO: OPTIMIZATION FOR DISJKTRA ONLY OR???
			if (flowCost > graph.costUpperBound) {
				return false;
			}
				graph.updatePotentials(dist);
			graph.calculateReducedCosts();
			sendFlow(violation, shortestPath, minUpperBound, li);
//				graph.print();
			return true;
		}

		// Shortest path from source node to dest node
		// Return the shortest path and cost through parameters,
		// or false as return value in case of no path
		bool findShortestPathNegativeCosts(unsigned int source, unsigned int dest, 
																			 vector<unsigned int>& path, vector<unsigned int>& dist) 
																			 const {
			vector<unsigned int> prev;
			
			findShortestPathsReducedCosts2(source, dest, prev, dist);
	
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
			return true;
		}

		void findShortestPathsReducedCosts2(unsigned int source, unsigned int dest, vector<unsigned int>& prev,
																			 vector<unsigned int>& dist) const {
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist.assign(graph.nodeList.size(), INF_UINT);
			dist[source] = 0;
			vector<bool> visited;
			visited.assign(graph.nodeList.size(), false);
		//	cout << "looking for " << source << " paths" << endl; 
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
					assert(dist[node] != INF_UINT);
					if (visited[edge.destNode]) {
						continue;
					}
					unsigned int newDist = dist[node] + edge.reducedCost;
				//	cout << "newDist of " << node << " " << newDist << endl;
					if (!(node == source && edge.destNode == dest) && 
						newDist < dist[edge.destNode]) {
						// Found path of lower cost
						dist[edge.destNode] = newDist;
						prev[edge.destNode] = node;
						heap.push(HeapItem(edge.destNode, newDist));
					}
				}
			}
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
					unsigned int newDist = dist[node] + edge.reducedCost;
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
			unsigned int saFlow = graph.getEdgeFlow(graph.sNode(), a);
			unsigned int sbFlow = graph.getEdgeFlow(graph.sNode(), b);  

			if (saFlow < edgeSA->upperBound 
					&& sbFlow > edgeSB->lowerBound 
					&& edgeAY->reducedCost > m - edgeYB->reducedCost) {
		//		cout << "\tCONDITION 1 EARLY PRUNNING " << a << " from " << y << endl; 
				return true;
			}
			return false; 
			unsigned int mB = mFactor(b, y, graph);
			if (mB != INF_UINT && saFlow < edgeSA->upperBound 
					&& sbFlow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mB > m) {
				//if (a == 10 && !y) return false;
				cout << "\tCONDITION 2 EARLY PRUNNING " << a << " (val "  << (*graph.nodeToVal)[a] << ") from " << y << endl;
				cout << saFlow << " " << edgeSA->upperBound << " " << sbFlow << " " << 
				edgeSB->lowerBound << " " << edgeAY->reducedCost << " " << mB << " " << m << endl; 
				//graph.print();
				//graph.printResidual();
				//if (a == 12 && !y) exit(1);
				return true;
			}
			unsigned int mA = mFactor(a, y, graph);
			if (mA != INF_UINT && saFlow == edgeSA->upperBound 
					&& sbFlow > edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA > m - edgeYB->reducedCost) {
		//		cout << "\tCONDITION 3 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			if (mA != INF_UINT && mB != INF_UINT && saFlow == edgeSA->upperBound 
					&& sbFlow == edgeSB->lowerBound 
					&& edgeAY->reducedCost + mA + mB > m) {
		//		cout << "\tCONDITION 4 EARLY PRUNNING " << a << " from " << y << endl;
				return true;
			}
			return false;
		}

	public:
		FlowGraphAlgorithms(FlowGraph& graph) : graph(graph) {}

		bool findMinCostFlow(LI* li) {
			pair<unsigned int, unsigned int> violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
				// cout << "Violation " << violation.first << "->" << violation.second << endl;
				// graph.print();
				// graph.printResidual();
				if (!minCostFlowIteration(violation, li)) {
				//	cout << "incosistent" << endl;
					return false;
				}
			}
			graph.oldFlowIsFeasible = true;
			//graph.calculateFlowCost(li);
			graph.flow.setFirstFlow(false);
		//	graph.addTResidualEdges();
			return graph.checkFlowCost();
		}

		// Given updatedEdges contains the edges whose bounds have been tightened
		// since last execution, do the following:
		// - Update the residual graph to match the changes
		// - If the old flow is not still feasible, find a new one, using the 
		//   incremental algorithm from the publication
		bool updateMinCostFlow(vector<EdgeUpdate>& updatedEdges, LI* li) {
			//    cout << "Propagate: update min cost flow" << endl;
			// graph.print();
			// graph.printResidual();
			//buildResidualGraph();
			// graph.printResidual();
			if (!updatedEdges.size()) {
				cout << "this shouldn't happen" << endl;
				assert(false);
				exit(1);
			}
			for (auto& e: updatedEdges) {
				// Among the edges that changed, look for one violating bounds
				// Assume violating lower bound initially
				// cout << e.src << "->" << e.dest << endl;
				unsigned int src = e.src;
				unsigned int dest = e.dest;
					auto res = graph.getEdge(e.src, e.dest);
					if (res == NULL) {
						continue;
					}
				if (graph.getEdgeFlow(e.src, e.dest)) {
					// Violating upper bound, swap direction of initial violating edge
						std::swap(src, dest);
						//src = e.dest;
						//dest = graph.tNode();
					if (!minCostFlowIteration({src, dest}, li)) {
						return false;
					}
				}
				graph.deleteEdge(e.src, e.dest);
				graph.deleteEitherResidualEdge(e.src, e.dest);
				int val = (*graph.nodeToVal)[e.src];
				graph.varToVals[e.dest].deleteVal(val);
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
															       vector<EdgeUpdate>& updatedEdges) {
			// vector<int> distances;
			// vector<unsigned int> prev;
			// bellmanFordShortestPaths(graph.tNode(), prev, distances, NULL);
		  // graph.calculateReducedCosts(distances);
			//graph.printResidual();
 
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
			auto& sNode = graph.nodeList[graph.sNode()];
			for (unsigned int i = 0; i < sNode.edgeListSize; i++) {
				auto& edge = (*sNode.edgeList)[i];
				if (graph.getEdgeFlow(graph.sNode(), edge.destNode) > 0) {
					vector<unsigned int> yList;
					unsigned int b = edge.destNode;
					unordered_set<unsigned int> targetNodes;
					vector< pair<unsigned int, unsigned int>> ayList;
					auto& bNode = graph.nodeList[b];
					for (unsigned int j = 0; j < bNode.edgeListSize; j++) {
						auto& edgeBY = (*bNode.edgeList)[j];
						if (graph.getEdgeFlow(b, edgeBY.destNode) == 1) {
							unsigned int y = edgeBY.destNode;
						//	minDist[y].bestVal = b;
							for (IntVarValues v(vars[y]); v(); ++v) {
								unsigned int a = (*graph.valToNode)[v.val()];
								if (a != b) {
									if (earlyPrune(a, b, y, graph.costUpperBound - 
															   (graph.flowCost))) {
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
																				graph.costUpperBound - (graph.flowCost));
					
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
						if ((int)reducedDistances[a] > (graph.costUpperBound - (graph.flowCost) 
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
				updatedEdges.push_back(EdgeUpdate(edge.src, edge.dest));
				// Prune
				GECODE_ME_CHECK(vars[edge.dest].nq(home, edge.val));
				// Also remove from varToVals
				// auto& vals = graph.varToVals[edge.dest];
				// vals.deleteVal(edge.val);
				// Update upper bound
				// graph.deleteEdge(edge.src, edge.dest);
				// assert(!graph.getEdgeFlow(edge.src, edge.dest));
				/*if (vars[edge.dest].assigned()) {
					// If a variable got assigned by pruning, set corresponding edge
					// lower bound to 1
					int assignedVal = vars[edge.dest].val();
					assert(*vals.begin() == assignedVal);
					auto valNode = graph.valToNode->find(assignedVal)->second;
					graph.getEdge(valNode, edge.dest)->lowerBound = 1;
				}*/
			}
			return ES_OK;
		}

};

#endif