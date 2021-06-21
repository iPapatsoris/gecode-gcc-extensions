#include <iostream>
#include <vector>

using namespace std;

/**
 * Base Edge class
 */
class Edge
{
protected:
	unsigned int destNode;
	unsigned int cost;

	Edge(unsigned int destNode, unsigned int cost) : destNode(destNode), cost(cost) {}
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
						 unsigned int upperBound, unsigned int cost)
			: Edge(destNode, cost), lowerBound(lowerBound), upperBound(upperBound),
				flow(0) {}
	NormalEdge() {}
	unsigned int getUpperBound() const { return upperBound; }
	unsigned int getCost() const { return cost; }
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
	bool isBackwards;

public:
	ResidualEdge(unsigned int destNode, unsigned int upperBound, int cost, int reducedCost, bool isBackwards)
			: Edge(destNode, cost), upperBound(upperBound), reducedCost(reducedCost), isBackwards(isBackwards) {}
	ResidualEdge(const NormalEdge &edge)
			: Edge(edge.getDestNode(), edge.getCost()), upperBound(edge.getUpperBound()),
				reducedCost(edge.getCost()), isBackwards(false) {}
	ResidualEdge() {}
	void print() const
	{
		cout << destNode << " upper " << upperBound << " cost " << cost
				 << " reduced cost " << reducedCost <<  "\n";
	}

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
	friend class OrderGraph;
};

// Edge containing both source node and destination info
// Normally source is not included in Edge class, because node ID N corresponds
// to the N-th position in the node list arrays
// It is typically used in updatedEdges vectors, which hold which edges got
// updated due to some pruning or assignment, and need to be checked on when
// we next update the residual graph. It is important to have NormalEdge object
// and not NormalEdge*, because when search finds a solution or fails, it will
// clone the graph and destroy the original, thus invalidating the pointer.
//typedef pair<unsigned int, NormalEdge> FullEdge;

class Node {
	vector<NormalEdge> edgeList;
	vector<ResidualEdge> residualEdgeList;
	int potential;

	Node(unsigned int totalEdges) : potential(0) {
		edgeList.reserve(totalEdges);
	}
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
	friend class OrderGraph;

};