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
	int cost;

	Edge(int destNode, int cost) : destNode(destNode), cost(cost) {}
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
	NormalEdge(int destNode, int lowerBound, int upperBound, int cost)
			: Edge(destNode, cost), lowerBound(lowerBound), upperBound(upperBound),
				flow(0) {}
	NormalEdge() {}
	int getUpperBound() const { return upperBound; }
	int getCost() const { return cost; }
	int getDestNode() const { return destNode; }
	void print() const {
		cout << destNode << " (" << lowerBound << "," << upperBound << "," << cost
				 << ") flow: " << flow << "\n";
	}

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif