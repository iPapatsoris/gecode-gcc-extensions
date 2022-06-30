#ifndef H_NODE
#define H_NODE

#include <vector>
#include "bt-vector.hpp"
#include "edge.hpp"
#include "residual-graph.hpp"

using namespace std;

class Node {
	// Efficient backtracking for edges
	BtVector<NormalEdge> edgeList;
	// No backtracking for residual graph; we build it from scratch each time

	Node(int totalEdges) : edgeList(totalEdges) {}
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif