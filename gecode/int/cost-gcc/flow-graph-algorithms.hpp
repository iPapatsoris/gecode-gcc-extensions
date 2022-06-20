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
#define MINUS_INF_INT INT_MIN
#define INF_UINT UINT_MAX
#define NONE_UINT INF_UINT - 1

using namespace std;
int debugCounter = 0;
bool neti = false;
class FlowGraphAlgorithms {

	public: 

	void debugPath() {
			unsigned int n;

			n = 0;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 1;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 2;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 3;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 4;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 5;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 6;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 7;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 8;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 9;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 10;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 11;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 12;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 13;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 14;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 15;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 16;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 17;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 18;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 19;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(26, 1, 1, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 20;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(5, 0, 1, 4));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(13, 0, 1, 4));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(1, 0, 1, 3));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(0, 0, 1, 10));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(15, 0, 1, 7));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(17, 0, 1, 5));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(19, 0, 1, 6));
			graph.nodeList[n].edgeList->back().flow = 0;

graph.nodeList[n].edgeList->push_back(NormalEdge(14, 0, 1, 4));
			graph.nodeList[n].edgeList->back().flow = 0;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 21;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(18, 0, 1, 8));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(19, 0, 1, 7));
			graph.nodeList[n].edgeList->back().flow = 0;

			graph.nodeList[n].edgeList->push_back(NormalEdge(14, 0, 1, 9));
			graph.nodeList[n].edgeList->back().flow = 0;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 22;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(10, 0, 1, 5));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(16, 0, 1, 7));
			graph.nodeList[n].edgeList->back().flow = 1;
			
			graph.nodeList[n].edgeList->push_back(NormalEdge(9, 0, 1, 7));
			graph.nodeList[n].edgeList->back().flow = 1;
			
			graph.nodeList[n].edgeList->push_back(NormalEdge(12, 0, 1, 8));
			graph.nodeList[n].edgeList->back().flow = 1;
			

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 23;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(7, 0, 1, 7));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(19, 0, 1, 4));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(14, 0, 1, 7));
			graph.nodeList[n].edgeList->back().flow = 1;



			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 24;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(8, 0, 1, 9));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(6, 0, 1, 6));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(4, 0, 1, 5));
			graph.nodeList[n].edgeList->back().flow = 1;

graph.nodeList[n].edgeList->push_back(NormalEdge(3, 0, 1, 2));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(2, 0, 1, 2));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(19, 0, 1, 6));
			graph.nodeList[n].edgeList->back().flow = 0;

			graph.nodeList[n].edgeList->push_back(NormalEdge(14, 0, 1, 2));
			graph.nodeList[n].edgeList->back().flow = 0;

						graph.nodeList[n].edgeList->push_back(NormalEdge(11, 0, 1, 8));
			graph.nodeList[n].edgeList->back().flow = 0;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 25;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(20, 2, 12, 0));
			graph.nodeList[n].edgeList->back().flow = 6;

			graph.nodeList[n].edgeList->push_back(NormalEdge(21, 1, 11, 0));
			graph.nodeList[n].edgeList->back().flow = 1;

			graph.nodeList[n].edgeList->push_back(NormalEdge(22, 4, 10, 0));
			graph.nodeList[n].edgeList->back().flow = 4;

			graph.nodeList[n].edgeList->push_back(NormalEdge(23, 3, 10, 0));
			graph.nodeList[n].edgeList->back().flow = 4;

			graph.nodeList[n].edgeList->push_back(NormalEdge(24, 0, 11, 0));
			graph.nodeList[n].edgeList->back().flow = 5;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			n = 26;
			graph.nodeList[n].edgeList->clear();
			graph.nodeList[n].edgeToPos->clear();

			graph.nodeList[n].edgeList->push_back(NormalEdge(25, 20, 20, 0));
			graph.nodeList[n].edgeList->back().flow = 20;

			graph.nodeList[n].edgeListSize = graph.nodeList[n].edgeList->size();
			for (unsigned int i = 0; i < graph.nodeList[n].edgeListSize; i++) {
				graph.nodeList[n].edgeToPos->insert({(*graph.nodeList[n].edgeList)[i].destNode, i});
			}

			buildResidualGraph(NULL);
			graph.print();
			graph.printResidual();

			vector<int> distances;
			vector<unsigned int> prev, cycle;
			bool isCycle = false;
				using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
		//exit(1);
		//unsigned int src = 11;
		//unsigned int dest = 21;

    auto t1 = high_resolution_clock::now();
			//bellmanFordShortestPathsCycles(graph.tNode(), prev, distances, cycle, &isCycle, NULL);
			 // bellmanFordShortestPathsCycles(src, prev, distances, cycle, &isCycle, &dest);
			auto t2 = high_resolution_clock::now();

    /* Getting number of milliseconds as an integer. */
    auto ms_int = duration_cast<milliseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_int.count() << "ms\n";

		if (cycle.size()) {
//						cout << "after cycle FOUND" << endl;
						isCycle = true;
						for (auto node: cycle) {
							cout << node << "->";
						}
						cout << endl;
		}

		}



	private:
		FlowGraph& graph;
		bool isMinCost;

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
			for (unsigned int i = 0; i < graph.tNode(); i++) {
				graph.nodeList[i].residualEdgeList->clear();
			}

			for (unsigned int i = graph.totalVarNodes; i < graph.tNode(); i++) {
				auto& node = graph.nodeList[i];
				for (unsigned int j = 0; j < node.edgeListSize; j++) {
					auto& edge = (*node.edgeList)[j];
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
     void bellmanFordShortestPaths(unsigned int source, 
																	vector<unsigned int>& prev, vector<int>& dist, 
																	unsigned int* dest = NULL) const {
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist.assign(graph.nodeList.size(), INF_INT);
			dist[source] = 0;
			bool debug = false;
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
			vector<unsigned int> updatedNodes1 = {source};
			vector<unsigned int> updatedNodes2;
			vector<unsigned int>* updatedNodesOld = &updatedNodes1;
			vector<unsigned int> *updatedNodesNew = &updatedNodes2;
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
			if (dest == NULL)
				;// cout << "done " << endl;
		}

		void traceCycle(const vector<unsigned int>& prev, unsigned int node, vector<unsigned int>& cycle) const {
				unordered_set<unsigned int> stackContent;
				stack<unsigned int> stack;

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

		void bellmanFordShortestPathsCycles(unsigned int source, 
																	vector<unsigned int>& prev, vector<int>& dist,
																	vector<unsigned int>& cycle, bool *isCycle,
																	unsigned int* dest = NULL) const {
			// if (dest != NULL)
			// cout << "look for " << source << "->" << *dest << endl;
			// else 
			// cout << "look for " << source << "->" << "all" << endl;
				
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist.assign(graph.nodeList.size(), INF_INT);
			dist[source] = 0;
			vector<unsigned int> len;
			len.assign(graph.nodeList.size(), 0);
		
			vector<unsigned int> updatedNodes1 = {source};
			vector<unsigned int> updatedNodes2;
			vector<unsigned int>* updatedNodesOld = &updatedNodes1;
			vector<unsigned int> *updatedNodesNew = &updatedNodes2;
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
							if (++len[edge.destNode] == graph.nodeList.size() - 1) {
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

	/*	void bellmanCycles(unsigned int source, 
																	vector<unsigned int>& prev, vector<int>& dist, 
																	unsigned int* dest = NULL) const {
			prev.assign(graph.nodeList.size(), NONE_UINT);
			dist.assign(graph.nodeList.size(), INF_INT);
			dist[source] = 0;
			bool debug = false;
			if (source == 21 && *dest == 5) {
				cout << "look for " << source << "->" << *dest << endl;
				graph.print();
				graph.printResidual();
				//exit(1);
				debug = true;
			}	

			for (unsigned int iterations = 0; iterations < graph.nodeList.size() - 1; iterations++) {
				bool foundUpdate = false;
				for (unsigned int node = 0; node < graph.nodeList.size(); node++) {
					for (auto& e: *graph.nodeList[node].residualEdgeList)
						if ((dest == NULL || 
							!(node == source && e.destNode == *dest)) && 
							dist[node] != INF_INT && dist[node] + e.cost < 
																			 dist[e.destNode]) {
							// Ignore direct source->dest edge when looking for shortest path 
							// to specific requested destination
							dist[e.destNode] = dist[node] + e.cost;
							prev[e.destNode] = node;
							foundUpdate = true;
					}
				}
				if (!foundUpdate) {
					break;
				}
			}

			for (unsigned int iterations = 0; iterations < graph.nodeList.size() - 1; iterations++) {
				bool foundUpdate = false;
				for (unsigned int node = 0; node < graph.nodeList.size(); node++) {
					for (auto& e: *graph.nodeList[node].residualEdgeList)
						if ((dest == NULL || 
							!(node == source && e.destNode == *dest)) && 
							dist[node] != INF_INT && dist[node] + e.cost < 
																			 dist[e.destNode]) {
							// Ignore direct source->dest edge when looking for shortest path 
							// to specific requested destination
							dist[e.destNode] = MINUS_INF_INT;
							prev[e.destNode] = node;
							foundUpdate = true;
					}
				}
				if (!foundUpdate) {
					break;
				}
			}

			//cout << "done " << endl;
		}*/

		void sendFlow(pair<unsigned int, unsigned int>& violation, const vector<unsigned int>& shortestPath, unsigned int minUpperBound, LI* li) {
			// Send flow through the path edges and update residual graph
			unsigned int prev = violation.first;
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

		void sendFlowCycle(const vector<unsigned int>& cycle, unsigned int minUpperBound, LI* li) {
			// Send flow through the path edges and update residual graph
			unsigned int prev = cycle[0];
			for(unsigned int i = 1; i < cycle.size(); i++) {			
				NormalEdge *edge = graph.getEdge(prev, cycle[i]);
				if (edge != NULL) {
					// Path residual edge is a forward edge in the original graph
					edge->flow += minUpperBound;
if (neti)					cout << "Flow of " << prev << " " << cycle[i] << " now " << edge->flow << endl;
					*graph.flowCost += edge->cost;
					if (edge->destNode < graph.totalVarNodes && li != NULL) {
						(*li)[edge->destNode] = (*graph.nodeToVal)[prev];
					}
					updateResidualGraph(prev, cycle[i], *edge);
				}	else {
					// Path residual edge is a backward edge in the original graph
					edge = graph.getEdge(cycle[i], prev);
					edge->flow -= minUpperBound;
if (neti)					cout << "Flow of " << cycle[i] << " " << prev << " now " << edge->flow << endl;
					if (!edge->flow) {
						*graph.flowCost -= edge->cost;
					}
					updateResidualGraph(cycle[i], prev, *edge);
				}
				prev = cycle[i];
			}
		}

		unsigned int findMinUpperBoundCycle(const vector<unsigned int>& cycle) const {
			unsigned int prev = cycle[0];
			unsigned int minUpperBound = INF_UINT;
if (neti)			cout << "SP cycle: ";
			for(unsigned int i = 1; i < cycle.size(); i++) {
				// Bellman returns the path in reverse, so traverse it in reverse
if (neti)				cout << prev << "->" << cycle[i]; 
				ResidualEdge *edge = graph.getResidualEdge(prev, cycle[i]);
				minUpperBound = min(minUpperBound, edge->upperBound);
				prev = cycle[i];
			}
if (neti)			cout << endl;
			return minUpperBound;
		}
		bool debug = false;
		unsigned int findMinUpperBound(pair<unsigned int, unsigned int>& violation, vector<unsigned int>& shortestPath, int* flowCost) {
			// Find min upper bound along shortest path
			unsigned int prev = violation.first;
			unsigned int minUpperBound = INF_UINT;
			// cout << "Violation " << prev << " " << violation.second << endl;
			if (neti) {
			cout << "SP: ";
			for (auto it = shortestPath.rbegin(); it != shortestPath.rend(); it++) {
				cout << *it << "->" << flush;
			}
			cout << endl;
			}
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

		bool minCostFlowIteration(pair<unsigned int, unsigned int> violation, bool *isCycle, LI* li, Int::IntView costUpperBound) {		
			vector<unsigned int> shortestPath;
			vector<int> dist;
			int pathCost; 
			//  cout << "Violation " << violation.first << "->" << violation.second << endl;
			//graph.print();
			//graph.printResidual();
			if (!findShortestPathNegativeCosts(violation.second, violation.first, 
																				 shortestPath, dist, pathCost, isCycle)) {
				// Constraint is not consistent
				return false;
			}

			if (isCycle != NULL && *isCycle) {
				unsigned int minUpperBound = findMinUpperBoundCycle(shortestPath);
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
			unsigned int minUpperBound = findMinUpperBound(violation, shortestPath, &flowCost);		
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
		bool findShortestPathNegativeCosts(unsigned int source, unsigned int dest, 
																			 vector<unsigned int>& path, vector<int>& dist,
																			int& cost, bool *isCycle) 
																			 const {
			vector<unsigned int> prev;
			
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

// if (debugCounter == 4569161 || debugCounter == 4569171 || debugCounter == 4569439 || debugCounter == 4569449) {
// 					vector<unsigned int> prev2, path2;
// 					vector<int> dist2;
// 					bool isCycle2 = false;
// 					cout << "check for cycle bug" << endl;
// 					bellmanFordShortestPathsCycles(graph.tNode(), prev2, dist2, path2, &isCycle2, NULL);
// 					if (path2.size()) {
// 						cout << "found cycle from T but not from normal path lol" << endl;
// 						exit(0);
// 					}
// 					exit(1);
// 				}

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

				for (auto& edge: *graph.nodeList[node].residualEdgeList) {
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
			unsigned int mB = mFactor(b, y, graph);
			if (mB != INF_UINT && edgeSA->flow < edgeSA->upperBound 
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

		bool findMinCostFlow(LI* li, Int::IntView costUpperBound) {
			pair<unsigned int, unsigned int> violation;
			while (graph.getLowerBoundViolatingEdge(violation)) {
//				cout << "Violation " << violation.first << "->" << violation.second << endl;
				if (!minCostFlowIteration(violation, NULL, li, costUpperBound)) {
				//	cout << "incosistent" << endl;
					return false;
				}
			}
			*graph.oldFlowIsFeasible = true;
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
		bool updateMinCostFlow(vector<EdgeUpdate>& updatedEdges, LI* li, Int::IntView costUpperBound) {
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
					if (neti) cout << "still cycle" << endl;
					vector<unsigned int> prev, path;
					vector<int> dist;
					bellmanFordShortestPathsCycles(graph.tNode(), prev, dist, path, &isCycle, NULL);
					if (path.size()) {
if (neti)						cout << "after cycle FOUND" << endl;
						isCycle = true;
/*						for (auto node: path) {
							cout << node << "->";
						}
						cout << endl;*/
						unsigned int minUpperBound = findMinUpperBoundCycle(path);
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

				unsigned int src = e.src;
				unsigned int dest = e.dest;
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
						if (neti) cout << "cycle that didn't fix violation, fixing now" << endl;
						debug = true;
						if (!minCostFlowIteration({src, dest}, &isCycle, li, costUpperBound)) {
							return false;
						}
					}
				}

				}
				
				graph.deleteEdge(e.src, e.dest);
				graph.deleteResidualEdge(e.src, e.dest);
				int val = (*graph.nodeToVal)[e.src];
				graph.varToVals[e.dest].deleteVal(val);
			
				while (isCycle) {
					isCycle = false;
					if (neti) cout << "still cycle" << endl;
					vector<unsigned int> prev, path;
					vector<int> dist;
					bellmanFordShortestPathsCycles(graph.sNode(), prev, dist, path, &isCycle, NULL);
					if (path.size()) {
if (neti)						cout << "after cycle FOUND" << endl;
						isCycle = true;
/*						for (auto node: path) {
							cout << node << "->";
						}
						cout << endl;*/
						unsigned int minUpperBound = findMinUpperBoundCycle(path);
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
					vector<unsigned int> prev, path;
					vector<int> dist;
					bellmanFordShortestPathsCycles(graph.sNode(), prev, dist, path, &isCycle, NULL);
					if (path.size()) {
						isCycle = true;
/*						for (auto node: path) {
							cout << node << "->";
						}
						cout << endl;*/
						unsigned int minUpperBound = findMinUpperBoundCycle(path);
						sendFlowCycle(path, minUpperBound, li);
						//assert(false);
						//exit(1);
					}
				}	while (isCycle);				
			}
			*(graph.oldFlowIsFeasible) = true;
			isMinCost = true;
			// cout << isMinCost << endl;

// 			bool isCycle = false;
// //					cout << "after cycle check" << endl;
// 					vector<unsigned int> prev, path;
// 					vector<int> dist;
// 					bellmanFordShortestPathsCycles(graph.sNode(), prev, dist, path, &isCycle, NULL);
// 					if (path.size()) {
// 						cout << "bug found" << endl;
// 					}
// 			graph.printResidual();
			// cout << "Update end" << endl;
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
															       vector<EdgeUpdate>& updatedEdges, LI* li, Int::IntView costUpperBound) {
		//	graph.addTResidualEdges(); // opt?
			vector<int> distances;
			vector<unsigned int> prev;
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
					vector<unsigned int> path;
					bellmanFordShortestPathsCycles(graph.tNode(), prev, distances, path, &isCycle, NULL);
					// cout << "after cycle done" << endl;
					if (path.size()) {
						// cout << "after cycle FOUND ARC" << endl;
						isCycle = true;
						// for (auto node: path) {
							// cout << node << "->";
						// }
						// cout << endl;
						unsigned int minUpperBound = findMinUpperBoundCycle(path);
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
				if (edge.flow > 0) {
					vector<unsigned int> yList;
					unsigned int b = edge.destNode;
					unordered_set<unsigned int> targetNodes;
					vector< pair<unsigned int, unsigned int>> ayList;
					auto& bNode = graph.nodeList[b];
					for (unsigned int j = 0; j < bNode.edgeListSize; j++) {
						auto& edgeBY = (*bNode.edgeList)[j];
						if (edgeBY.flow == 1) {
							unsigned int y = edgeBY.destNode;
						//	minDist[y].bestVal = b;
							for (IntVarValues v(vars[y]); v(); ++v) {
								unsigned int a = (*graph.valToNode)[v.val()];
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

					vector<unsigned int> reducedDistances;
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
						unsigned int costAY = residualEdge->reducedCost;
						unsigned int costYB = graph.getResidualEdge(y, b)->reducedCost;
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
				NormalEdge* actualEdge = graph.getEdge(edge.src, edge.dest);
				assert(actualEdge != NULL);
				// Push to updatedEdges so we can modify the residual graph accordingly
				// on the next min cost flow computation
	//			updatedEdges.push_back(EdgeUpdate(edge.src, edge.dest, false, false, true));
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
															       vector<EdgeUpdate>& updatedEdges, Int::IntView costUpperBound) {
		//	graph.addTResidualEdges(); // opt?
			vector<int> distances;
			vector<unsigned int> prev;
   // using std::chrono::milliseconds;

    // auto t1 = high_resolution_clock::now();
				bellmanFordShortestPaths(graph.tNode(), prev, distances, NULL);
 
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
			auto& sNode = graph.nodeList[graph.sNode()];
			for (unsigned int i = 0; i < sNode.edgeListSize; i++) {
				auto& edge = (*sNode.edgeList)[i];
				if (edge.flow > 0) {
					vector<unsigned int> yList;
					unsigned int b = edge.destNode;
					unordered_set<unsigned int> targetNodes;
					vector< pair<unsigned int, unsigned int>> ayList;
					auto& bNode = graph.nodeList[b];
					for (unsigned int j = 0; j < bNode.edgeListSize; j++) {
						auto& edgeBY = (*bNode.edgeList)[j];
						if (edgeBY.flow == 1) {
							unsigned int y = edgeBY.destNode;
						//	minDist[y].bestVal = b;
							for (IntVarValues v(vars[y]); v(); ++v) {
								unsigned int a = (*graph.valToNode)[v.val()];
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
						unsigned int costAY = residualEdge->cost;
						unsigned int costYB = graph.getResidualEdge(y, b)->cost;
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
				NormalEdge* actualEdge = graph.getEdge(edge.src, edge.dest);
				assert(actualEdge != NULL);
				// Push to updatedEdges so we can modify the residual graph accordingly
				// on the next min cost flow computation
	//			updatedEdges.push_back(EdgeUpdate(edge.src, edge.dest, false, false, true));
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