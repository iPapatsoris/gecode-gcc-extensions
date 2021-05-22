#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

using namespace Gecode;
using namespace std;

typedef unordered_map<int, unordered_set<int>> MapToSet;

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
	unsigned int cost;
	unsigned int flow;

	public:
		NormalEdge(unsigned int destNode, unsigned int lowerBound, 
						   unsigned int upperBound, unsigned int cost) 
				: Edge(destNode), lowerBound(lowerBound), upperBound(upperBound), 
					cost(cost), flow(0) {}
		NormalEdge() {}
		unsigned int getUpperBound() const { return upperBound; }
		unsigned int getCost() const { return cost; }
		unsigned int getDestNode() const { return destNode; }
		void print() const {
			cout << destNode << " (" << lowerBound << "," << upperBound << "," << cost 
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
	int cost;

	public: 
		ResidualEdge(unsigned int destNode, unsigned int upperBound, int cost) 
				: Edge(destNode), upperBound(upperBound), cost(cost) {}
		ResidualEdge(const NormalEdge& edge) 
				: Edge(edge.getDestNode()), upperBound(edge.getUpperBound()), 
					cost(edge.getCost()) {}
		ResidualEdge() {}
		void print() const {
			cout << destNode << " upper " << upperBound << " cost " << cost << "\n";
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

/**
 * Graph used to solve min cost flow problem, to prove consistency of
 * costgcc, and also achieve arc consistency later on.
 * There is a node for each variable, for each value, and an S and T node.
 * There are value->variable edges according to the respective domains,
 * variable->T edges, S->value edges and a T->S edge.
 * The lower/upper bounds on the edges are chosen in a way to respect 
 * the original contraints, while also making sure that each variable gets
 * eventually assigned to a value.
 * 
 * The graph is built once, and later on any changes in domain such as value 
 * prunes or variable assignments, are expressed by setting value->variable 
 * upper bounds to 0 (prune) or lower bounds to 1 (assignment)
 * 
 * The graph's residual graph exists through Node's residualEdgeList member.
 * It is calculated once on construction, and later on when we have a new flow,
 * it is only modified in the edges that differ from the previous flow.
 */
class FlowGraph {
	friend class FlowGraphAlgorithms;
	private:
		// The node order in nodeList is variables,values,S,T.
		vector<Node> nodeList;
		unsigned int totalVarNodes; 
		// Fast lookup of what node a value corresponds to
		unordered_map<int, unsigned int> valToNode;

		// Total flow through the graph, starts at 0
		int flowCost;

		// Position of S node
		unsigned int sNode() const { return nodeList.size() - 2; }
		// Position of T node
		unsigned int tNode() const { return nodeList.size() - 1; }

		// Search for an edge flow violating lower bounds
		// Return false if none exists
		bool getLowerBoundViolatingEdge(pair<unsigned int, NormalEdge>& violation) 
			const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				for (auto& edge: nodeList[i].edgeList) {
					if (edge.flow < edge.lowerBound) {
						violation = {i, edge};
						return true;
					}
				}
			}
			return false;
		}

		// To compute reduced costs on the residual graph, we need a starting node
		// from which we can reach every other node. Use the T node for this.
		// The final residual graph of the graph where we find a feasible min cost 
		// flow, does not contain T edges, so we restore them with this function
		void addTResidualEdges() {
			for (unsigned int var = 0; var < totalVarNodes; var++) {
				nodeList[tNode()].residualEdgeList.push_back(ResidualEdge(var, 1, 0));
			}
		}

		// Search for source->dest edge, return pointer to it or NULL if it 
		// doesn't exist
		NormalEdge* getEdge(unsigned int source, unsigned int dest) {
			for (auto& edge: nodeList[source].edgeList) {
				if (edge.destNode == dest) {
					return &edge;
				}
			}
			return NULL;
		}
		
		// Search for source->dest residual edge, return pointer to it or NULL if it
		// doesn't exists. If index is not NULL, also return its position in that
		// node's residual edges list 
		ResidualEdge* getResidualEdge(unsigned int source, unsigned int dest, 
																  unsigned int *index = NULL) {
			for (unsigned int i=0; i<nodeList[source].residualEdgeList.size(); i++) {
				ResidualEdge& edge = nodeList[source].residualEdgeList[i];
				if (edge.destNode == dest) {
					if (index != NULL) {
						*index = i;
					}
					return &edge;
				}
			}
			return NULL;
		}

		// If existingEdge is NULL, create edge from source to newEdge.dest
		// If existingEdge is not NULL, change it to newEdge
		void setOrCreateResidualEdge(ResidualEdge* existingEdge, 
															   unsigned int source, 
																 const ResidualEdge& newEdge) {
			if (existingEdge != NULL) {
				*existingEdge = newEdge;
			} else {
				nodeList[source].residualEdgeList.push_back(newEdge);
			}
		}

		// Iterate through each edge that has flow, to find its total cost
		int calculateFlowCost() {
			for (unsigned int i = totalVarNodes; i < sNode(); i++) {
				for (auto& edge: nodeList[i].edgeList) {
					if (edge.flow > 0) {
						flowCost += edge.cost;
					}
				}
			}

			return flowCost;
		}

		// Given 'distances' contains the shortest paths from one node to every 
		// other node in the graph, transform the costs of the residual graph to 
		// reduced costs, making sure that they all become non-negative
		void calculateReducedCosts(const vector<int>& distances) {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				for (auto& edge : nodeList[i].residualEdgeList) {
					edge.cost = distances[i] + edge.cost - distances[edge.destNode];
				}
			}
		}

	public:
		// Create FlowGraph. Parameter includePruned controls whether early pruned
		// values without any matching variables, should be included as edge-less
		// nodes or not. If propagator post function calls checkLowerBounds,
		// includePruned should be false, since the lower bounds of the pruned
		// values have already been checked for validity.
		// If checkLowerBounds isn't called, it should be true, to include them 
		// and let the min cost flow algorithm check the lower bounds restrictions.
		// This is only for testing purposes, to see how much checkLowerBounds helps
		// performance. 
		// TODO: Eventually, decide on one way and remove this parameter.
		FlowGraph(const ViewArray<Int::IntView>& vars, const MapToSet& valToVars,
					const IntArgs& inputVals, const IntArgs& lowerBounds, 
					const IntArgs& upperBounds, const IntArgs& costs, bool includePruned) 
				: flowCost(0) {

			totalVarNodes = vars.size();
			unsigned int totalValNodes = (includePruned ? inputVals.size() 
																									: valToVars.size());
			// Nodes are variable nodes, values nodes, S and T nodes
			int totalNodes = totalVarNodes + totalValNodes + 2;
			// S node position
			int sNode = totalNodes - 2;
			// T node position
			int tNode = totalNodes - 1;
			nodeList.reserve(totalNodes);

			// Insert variable nodes and var->T edges
			for (unsigned int x = 0; x < totalVarNodes; x++) {
				nodeList.push_back(Node(1));
				nodeList.back().edgeList.push_back(NormalEdge(tNode, 1, 1, 0));
			}

			Matrix<IntArgs> c(costs, inputVals.size(), vars.size());

			// Insert Value nodes and Val->Var edges
			// It is important to iterate through inputVals and not through the
			// domains or valToVars, because values might have been early pruned from
			// the latter. We still need to include pruned values, to respect
			// their lower bound restriction
			for (int i = 0; i < inputVals.size(); i++) {
				int val = inputVals[i];
				auto it = valToVars.find(val);
				if (it != valToVars.end() || includePruned) {
					valToNode.insert({val, nodeList.size()});
					cout << "node " << nodeList.size() << " corresponds to val " << val
							 << "\n";
					if (it != valToVars.end()) {
						nodeList.push_back(Node(it->second.size()));
						for (auto& var : it->second) {
							int lowerBound = (vars[var].assigned() ? 1 : 0);
							nodeList.back().edgeList.push_back(NormalEdge(var, lowerBound, 1, 
																								c(i, var)));
						}
					} else if (includePruned) {
						// This value has been pruned early from all possible variables
						// Insert it in the graph to respect its lower bound restriction,
						// but do not add any edges to it
						nodeList.push_back(Node(0));
					}
				}
			}

			// Insert S node and S->Val edges
			nodeList.push_back(Node(inputVals.size()));
			for (int i = 0; i < inputVals.size(); i++) {
				int val = inputVals[i];
				auto valNode = valToNode.find(val);
				if (valNode != valToNode.end()) {
					nodeList.back().edgeList.push_back(NormalEdge(valNode->second, 
																				    lowerBounds[i], upperBounds[i], 0));
				}
			}

			// Insert T node and T->S edge
			nodeList.push_back(Node(1));
			nodeList.back().edgeList.push_back(NormalEdge(sNode, totalVarNodes, 
																				 totalVarNodes, 0));

			// Create residual graph
			for (auto &node : nodeList) {
				copy(node.edgeList.begin(), node.edgeList.end(), 
						 back_inserter(node.residualEdgeList));
			}
		}

		int getFlowCost() const {
			return flowCost;
		}

		void print() const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				auto& node = nodeList[i];
				for (unsigned int j = 0; j < node.edgeList.size(); j++) {
					auto& edge = node.edgeList[j];
					cout << i << " -> ";
					edge.print();
				}
			}
			cout << endl;
		}

		void printResidual() const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				auto& node = nodeList[i];
				for (unsigned int j = 0; j < node.residualEdgeList.size(); j++) {
					auto& edge = node.residualEdgeList[j];
					cout << i << " -> ";
					edge.print();
				}
			}
			cout << endl;
		}
};