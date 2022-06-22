#ifndef H_NODE
#define H_NODE

#include <iostream>
#include <vector>
#include "bt-vector.hpp"
#include "edge.hpp"

using namespace std;

class Node {
	BtVector<NormalEdge, unsigned int> edgeList;
	vector<ResidualEdge> *residualEdgeList;

	Node(unsigned int totalEdges) : edgeList(totalEdges) {
		residualEdgeList = new vector<ResidualEdge>();
	}
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif