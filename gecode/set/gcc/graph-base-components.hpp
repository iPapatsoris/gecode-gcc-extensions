#ifndef GRAPH_BASE_COMPONENTS
#define GRAPH_BASE_COMPONENTS

#include <iostream>
#include <vector>
#include <climits>

using namespace std;

/**
 * Base Edge class
 */
class Edge {
protected:
	unsigned int destNode;

	Edge(unsigned int destNode) : destNode(destNode) {}
	Edge() {}

	friend class FlowGraphAlgorithms;
};

/**
 * Regular Edge for the original flow graph
 */
class NormalEdge : public Edge {
	unsigned int lowerBound;
	unsigned int upperBound;
	unsigned int flow;

public:
	NormalEdge(unsigned int destNode, unsigned int lowerBound,
						 unsigned int upperBound)
			: Edge(destNode), lowerBound(lowerBound), upperBound(upperBound),
				flow(0) {}
	NormalEdge() {}
	unsigned int getUpperBound() const { return upperBound; }
	unsigned int getDestNode() const { return destNode; }
	void print() const {
		cout << destNode << " (" << lowerBound << "," << upperBound
				 << ") flow: " << flow << "\n";
	}

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

/**
 * Edge for residual graph
 */
class ResidualEdge : public Edge {
	unsigned int upperBound;

public:
	ResidualEdge(unsigned int destNode, unsigned int upperBound)
			: Edge(destNode), upperBound(upperBound) {}
	ResidualEdge(const NormalEdge &edge)
			: Edge(edge.getDestNode()), upperBound(edge.getUpperBound()) {}
	ResidualEdge() {}
	void print() const {
		cout << destNode << " upper " << upperBound <<  "\n";
	}

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

class Node {
	vector<NormalEdge> edgeList;
	vector<ResidualEdge> residualEdgeList;

	Node(unsigned int totalEdges) {
		edgeList.reserve(totalEdges);
	}
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif