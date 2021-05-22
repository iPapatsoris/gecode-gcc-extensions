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
		void updateResidualGraph(unsigned int source, unsigned int dest, NormalEdge* edge) {
			unsigned int residualEdgeIndex;
			ResidualEdge *residualEdgeSearch = graph.getResidualEdge(source, dest, &residualEdgeIndex);
			if (edge->flow < edge->upperBound) {
				// Add / update forward residual edge
				graph.setOrCreateResidualEdge(residualEdgeSearch, source, ResidualEdge(dest, edge->upperBound - edge->flow, edge->cost));
			} else if (residualEdgeSearch != NULL) {
				// Delete forward residual edge that should no longer exist
				graph.nodeList[source].residualEdgeList.erase(graph.nodeList[source].residualEdgeList.begin() + residualEdgeIndex);
			}

			residualEdgeSearch = graph.getResidualEdge(dest, source, &residualEdgeIndex);
			if (edge->flow > edge->lowerBound) {
				// Add / update backward residual edge
				graph.setOrCreateResidualEdge(residualEdgeSearch, dest, ResidualEdge(source, edge->flow - edge->lowerBound, -edge->cost));
			} else if (residualEdgeSearch != NULL) {
				// Delete backward residual edge that should no longer exist
				graph.nodeList[dest].residualEdgeList.erase(graph.nodeList[dest].residualEdgeList.begin() + residualEdgeIndex);
			}
		}

		// Bellman-Ford algorithm for shortest paths with negative costs.
		// If dest is not NULL, ignore any direct source->dest edge.
		// This is needed when searching for shortest path to a specific 
		// destination, by the min cost flow algorithm.
		// TODO: consider testing the randomized variation improvement
		void bellmanFordShortestPaths(unsigned int source, 
																	vector<unsigned int>& prev, vector<int>& dist, 
																	unsigned int* dest = NULL) const {
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist.assign(graph.nodeList.size(), INF_INT);
			dist[source] = 0;

			for (unsigned int i = 0; i < graph.nodeList.size() - 1; i++) {
				bool foundUpdate = false;
				for (unsigned int node = 0; node < graph.nodeList.size(); node++) {
					for (auto &edge : graph.nodeList[node].residualEdgeList) {
						if (dest == NULL || (
							!(node == source && edge.destNode == *dest) && 
							dist[node] != INF_INT && dist[node] + edge.cost < 
																			 dist[edge.destNode])) {
							// Ignore direct source->dest edge when looking for shortest path 
							// to specific requested destination
							dist[edge.destNode] = dist[node] + edge.cost;
							prev[edge.destNode] = node;
							foundUpdate = true;
						}
					}
				}
				if (!foundUpdate) {
					// No updates between two iterations; early termination
					break;
				}
			}
		}

	public:
		FlowGraphAlgorithms(FlowGraph& graph) : graph(graph) {}

		// Shortest paths from source node to all other nodes
		// Return the costs for each shortest path in distances
		// To be used when graph can contain negative costs, thus for the residual
		void findShortestPathsNegativeCosts(unsigned int source, 
																			  vector<int> &distances) const {
			vector<unsigned int> prev;
			bellmanFordShortestPaths(source, prev, distances);
		}

		// Shortest path from source node to dest node
		// Return the shortest path and cost through parameters,
		// or false as return value in case of no path
		bool findShortestPathNegativeCosts(unsigned int source, unsigned int dest, 
																			 vector<unsigned int>& path, int& cost) 
																			 const {
			vector<unsigned int> prev;
			vector<int> dist;
			
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

		bool findMinCostFlow() {
			pair<unsigned int, NormalEdge> violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
				/*cout << "Edge violation " << violation.first << "->" 
						 << violation.second.destNode << " flow " << violation.second.flow 
						 << " lower bound " << violation.second.lowerBound << "\n";*/
				vector<unsigned int> shortestPath;
				int pathCost; 
				if (!findShortestPathNegativeCosts(violation.second.destNode, 
																					 violation.first, shortestPath, 
																					 pathCost)) {
					// Constraint is not consistent
					cout << "No path for edge violation " << violation.first << "->" 
							 << violation.second.destNode << endl;
					graph.print();
					graph.printResidual();
					return false;
				}

				//cout << "Path cost " << pathCost << endl;

				// Find min upper bound along shortest path
				unsigned int prev = violation.first;
				unsigned int minUpperBound = INF_UINT;
				for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
					// Bellman returns the path in reverse, so traverse it in reverse
					ResidualEdge *edge = graph.getResidualEdge(prev, *it);
					minUpperBound = min(minUpperBound, edge->upperBound);
					prev = *it;
					cout << *it << (it != shortestPath.rend()-1 ? "->" : "\n");
				}
				//cout << "\nk = " << minUpperBound << endl;
				assert(minUpperBound != INF_UINT);

				// Send flow through the path edges and update residual graph
				prev = violation.first;
				for(auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
					NormalEdge *edge = graph.getEdge(prev, *it);
					if (edge != NULL) {
						// Path residual edge is a forward edge in the original graph
						edge->flow += minUpperBound;
						updateResidualGraph(prev, *it, edge);
					}	else {
						// Path residual edge is a backward edge in the original graph
						cout << "\tBackwards edge flow!!!!!!!" << endl;
						edge = graph.getEdge(*it, prev);
						edge->flow -= minUpperBound;
						updateResidualGraph(*it, prev, edge);
					}
					prev = *it;
				}
			}
			graph.calculateFlowCost();
			cout << "Min cost flow " << graph.getFlowCost() << endl;
			return true;
		}
};