#ifndef H_NODE
#define H_NODE

#include <vector>
#include "bt-vector.hpp"
#include "edge.hpp"

using namespace std;

class Node {
	// Efficient backtracking for edges
	BtVector<NormalEdge> edgeList;
	// No backtracking for residual graph; we build it from scratch each time
	shared_ptr<vector<ResidualEdge>> residualEdgeList;

	Node(int totalEdges) : edgeList(totalEdges) {
		residualEdgeList = make_shared<vector<ResidualEdge>>();
	}
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif