#ifndef H_EDGE
#define H_EDGE

#include <iostream>

using namespace std;

/**
 * Base Edge class
 */
class Edge
{
protected:
	int destNode;

	Edge(int destNode) : destNode(destNode) {}
	Edge() {}	

	friend class FlowGraphAlgorithms;
};

/**
 * Regular Edge for the original flow graph
 */
class NormalEdge : public Edge
{
	int lowerBound;
	int upperBound;
	int flow;

public:
	NormalEdge(int destNode, int lowerBound, int upperBound)
			: Edge(destNode), lowerBound(lowerBound), upperBound(upperBound),
				flow(0) {}
	NormalEdge() {}
	int getUpperBound() const { return upperBound; }
	int getDestNode() const { return destNode; }
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
class ResidualEdge : public Edge
{
	int upperBound;

public:
	ResidualEdge(int destNode, int upperBound)
			: Edge(destNode), upperBound(upperBound) {}
	ResidualEdge(const NormalEdge &edge)
			: Edge(edge.getDestNode()), upperBound(edge.getUpperBound()) {}
	ResidualEdge() {}
	void print() const {
		cout << destNode << " upper " << upperBound;
	}

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif