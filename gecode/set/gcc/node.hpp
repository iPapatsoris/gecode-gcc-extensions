#ifndef H_NODE
#define H_NODE

#include <vector>
#include "bt-vector.hpp"
#include "edge.hpp"

using namespace std;

class Node {
	BtVector<NormalEdge> edgeList;
	vector<ResidualEdge> residualEdgeList;

	Node(int totalEdges) : edgeList(totalEdges) {}
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif