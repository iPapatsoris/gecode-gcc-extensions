#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

using namespace Gecode;
using namespace std;

template <class T1, class T2>
struct MapToSet {
	unordered_map<T1, unordered_set<T2>> map;
	MapToSet() {}
	MapToSet(const MapToSet &c) {
		map = c.map;
	}
};

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
	unsigned int reducedCost;

	public: 
		ResidualEdge(unsigned int destNode, unsigned int upperBound, int cost) 
				: Edge(destNode), upperBound(upperBound), cost(cost), reducedCost(0) {}
		ResidualEdge(const NormalEdge& edge) 
				: Edge(edge.getDestNode()), upperBound(edge.getUpperBound()), 
					cost(edge.getCost()), reducedCost(0) {}
		ResidualEdge() {}
		void print() const {
			cout << destNode << " upper " << upperBound << " cost " << cost 
					 << " reduced cost " << reducedCost << "\n";
		}
	
	friend class FlowGraph;
	friend class FlowGraphAlgorithms;
};

// Edge containing both source node and destination info
// Normally source is not included in Edge class, because node ID N corresponds
// to the N-th position in the node list arrays
typedef pair<unsigned int, NormalEdge *> FullEdge;

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
		// Fast lookup of what value a node corresponds to
		unordered_map<unsigned int, int> nodeToVal;
		// Holds each variable's domain. Is needed to find out which values got
		// pruned during an iteration, by comparing with Gecode's variable domains.
		// We know which variables got changed by using advisors, so domain
		// comparison is only done among those.
		// Should be kept up to date with assignments and pruning.
		MapToSet<unsigned int, int> varToVals;

		// Total flow through the graph, starts at 0. Is calculated at once using
		// appropriate function, not gradually
		int flowCost;

		// Cost upper bound as defined by the constraint input
		const int costUpperBound;

		// When values are pruned or variables are assigned, we update the 
		// bounds of the corresponding edges. At that point, set this variable
		// to note whether our old flow is still feasible, or if we need to find a 
		// new one. It is much more efficient to check this when we update 
		// the bounds of specific edges, than to scan the whole graph later to see 
		// if the old flow still stands
		bool oldFlowIsFeasible;

		// Position of S node
		unsigned int sNode() const { return nodeList.size() - 2; }
		// Position of T node
		unsigned int tNode() const { return nodeList.size() - 1; }

		// Search for an edge flow violating lower bounds
		// Return false if none exists
		bool getLowerBoundViolatingEdge(pair<unsigned int, unsigned int>& violation) 
			const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				for (auto& edge: nodeList[i].edgeList) {
					if (edge.flow < edge.lowerBound) {
						violation = {i, edge.destNode};
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
			flowCost = 0;
			for (unsigned int i = totalVarNodes; i < sNode(); i++) {
				for (auto& edge: nodeList[i].edgeList) {
					if (edge.flow > 0) {
						flowCost += edge.cost;
					}
				}
			}
			return flowCost;
		}

		bool checkFlowCost() const {
			cout << "Min cost flow " << flowCost << endl;
			return flowCost <= costUpperBound;
		}

		// Given 'distances' contains the shortest paths from one node to every 
		// other node in the graph, transform the costs of the residual graph to 
		// reduced costs, making sure that they all become non-negative
		void calculateReducedCosts(const vector<int>& distances) {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				for (auto& edge : nodeList[i].residualEdgeList) {
					edge.reducedCost = distances[i] + edge.cost - distances[edge.destNode];
				}
			}
		}

		// Assert varToVals is synchronized with Gecode variable X domain
		void assertVarToValsInSync(Int::IntView x, int xIndex) const {
			auto vals = varToVals.map.find(xIndex)->second;
			assert(vals.size() == x.size());
			for (IntVarValues v(x); v(); ++v) {
				assert(vals.find(v.val()) != vals.end());
			}
			for (auto val: vals) {
				assert(x.in(val));
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
		FlowGraph(
			const ViewArray<Int::IntView>& vars, 
 			const MapToSet<unsigned int, int>& varToVals,
			const MapToSet<int, unsigned int>& valToVars,
			const IntArgs& inputVals, const IntArgs& lowerBounds, 
			const IntArgs& upperBounds, const IntArgs& costs, int costUpperBound, 
			bool includePruned) 
				: flowCost(0), costUpperBound(costUpperBound) {
			
			this->varToVals = varToVals;
			totalVarNodes = vars.size();
			unsigned int totalValNodes = (includePruned ? inputVals.size() 
																									: valToVars.map.size());
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
				auto it = valToVars.map.find(val);
				if (it != valToVars.map.end() || includePruned) {
					valToNode.insert({val, nodeList.size()});
					nodeToVal.insert({nodeList.size(), val});
					//cout << "node " << nodeList.size() << " corresponds to val " << val
					//		 << "\n";
					if (it != valToVars.map.end()) {
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

		// Update graph state to match variable X domain pruning/assignment.
		// Update is made by tightening the bounds of edge V->X as follows:
		// - If X got assigned to value V, set the lower bound to 1.
		// - For every value V that has been pruned off X, set the upper bound 
		//   to 0. 
	  // If we prune a value that is used by current flow, or assign a value 
		// that is not used by it, set oldFlowIsFeasible to false.
		void updatePrunedValues(Int::IntView x, unsigned int xIndex, 
													  vector<FullEdge>& updatedEdges) {
			oldFlowIsFeasible = true;
			// Hold iterators to the values that we end up prunning, so we can also
			// remove them from valToVars
			vector< unordered_set<int>::iterator > prunedValues; 
			// Iterate all values for which there is an edge to X
			auto& values = varToVals.map.find(xIndex)->second;
			for (auto valueIt = values.begin(); valueIt != values.end(); valueIt++) {
				auto value = *valueIt;
				auto valueNode = valToNode.find(value)->second;
				NormalEdge* edge = getEdge(valueNode, xIndex);
				if (!edge->lowerBound && edge->upperBound == 1) {
					// If edge hasn't already been pruned or assigned
					if (!x.in(value)) {
						// Value has been pruned from variable X's domain, update graph
						edge->upperBound = 0;
						if (edge->flow == 1) {
							oldFlowIsFeasible = false;
						}
						updatedEdges.push_back({valueNode, edge});
						prunedValues.push_back(valueIt);
					}
					if (x.assigned() && x.val() == value && !edge->lowerBound) {
						// Variable has been assigned with a value, update graph
						edge->lowerBound = 1;
						if (edge->flow == 0) {
							oldFlowIsFeasible = false;
						}
						updatedEdges.push_back({valueNode, edge});
					}	
				}
			}

			for (auto it: prunedValues) {
				values.erase(it);
			}

			#ifndef NDEBUG
				assertVarToValsInSync(x, xIndex);
			#endif
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