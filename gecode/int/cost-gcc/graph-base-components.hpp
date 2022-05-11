#ifndef GRAPH_BASE_COMPONENTS
#define GRAPH_BASE_COMPONENTS

#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

/**
 * Base Edge class
 */
class Edge
{
protected:
	unsigned int destNode;
	int cost;

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

class Node {
	vector<NormalEdge> *edgeList;
	unordered_map<unsigned int, unsigned int> *edgeToPos;
	unsigned int edgeListSize;

	vector<ResidualEdge> *residualEdgeList;

	Node(unsigned int totalEdges) : edgeListSize(totalEdges) {
		edgeList = new vector<NormalEdge>();
		edgeToPos = new unordered_map<unsigned int, unsigned int>();
		edgeToPos->reserve(totalEdges);
		edgeList->reserve(totalEdges);
		residualEdgeList = new vector<ResidualEdge>();
	}
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;

	void print() const {
		cout << "edgeListSize: " << edgeListSize << "\n";
		cout << "edgeToPos:";
		for (auto& p: *edgeToPos) {
			cout << p.first << ": " << p.second << " ";
		}
		cout << "edgeList:\n";
		for (auto& e: *edgeList) {
			cout << e.getDestNode() << " ";
		}
		cout << endl;

	}
};

#endif