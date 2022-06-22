#ifndef H_EDGE
#define H_EDGE

#include <iostream>
#include <vector>
#include "bt-vector.hpp"

using namespace std;

/**
 * Base Edge class
 */
class Edge
{
protected:
	unsigned int destNode;
	int cost;

	Edge(unsigned int destNode, unsigned int cost) : destNode(destNode), 
																									 cost(cost) {}
	Edge() {}	

	friend class FlowGraphAlgorithms;
};

/**
 * Regular Edge for the original flow graph
 */
class NormalEdge : public Edge
{
	unsigned int lowerBound;
	unsigned int upperBound;
	unsigned int flow;

public:
	NormalEdge(unsigned int destNode, unsigned int lowerBound,
						 unsigned int upperBound, int cost)
			: Edge(destNode, cost), lowerBound(lowerBound), upperBound(upperBound),
				flow(0) {}
	NormalEdge() {}
	unsigned int getUpperBound() const { return upperBound; }
	int getCost() const { return cost; }
	unsigned int getDestNode() const { return destNode; }
	void print() const
	{
		cout << destNode << " (" << lowerBound << "," << upperBound << "," << cost
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
	unsigned int upperBound;
	unsigned int reducedCost;

public:
	ResidualEdge(unsigned int destNode, unsigned int upperBound, int cost)
			: Edge(destNode, cost), upperBound(upperBound), reducedCost(0) {}
	ResidualEdge(const NormalEdge &edge)
			: Edge(edge.getDestNode(), edge.getCost()), upperBound(edge.getUpperBound()) {}
	ResidualEdge() {}
	void print() const
	{
		cout << destNode << " upper " << upperBound
				 << " reduced cost " << reducedCost << " cost " << cost << "\n"; 
	}

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

#endif