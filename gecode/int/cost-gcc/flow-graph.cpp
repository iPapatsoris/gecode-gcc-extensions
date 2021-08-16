#include "flow-graph.hpp"
#include "example/cost-gcc-example.hpp"

		FlowGraph::FlowGraph(
			const ViewArray<Int::IntView>& vars, 
 			const MapToSet<unsigned int, int>& varToVals,
			const MapToSet<int, unsigned int>& valToVars,
			const IntArgs& inputVals, const IntArgs& lowerBounds, 
			const IntArgs& upperBounds, const IntArgs& costs, int costUpperBound) 
				: flowCost(0), costUpperBound(costUpperBound), oldFlowIsFeasible(true), 
				firstTime(true) {
			
			this->varToVals = varToVals;
			totalVarNodes = vars.size();
			unsigned int totalValNodes = inputVals.size();
			// Nodes are variable nodes, values nodes, S and T nodes
			int totalNodes = totalVarNodes + totalValNodes + 2;
			// S node position
			int sNode = totalNodes - 2;
			// T node position
			int tNode = totalNodes - 1;
			//nodeList = new vector<Node>();
			nodeToVal = new unordered_map<unsigned int, int>();
			valToNode = new unordered_map<int, unsigned int>();
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
				if (it != valToVars.map.end()) {
					valToNode->insert({val, nodeList.size()});
					nodeToVal->insert({nodeList.size(), val});
					//cout << "node " << nodeList.size() << " corresponds to val " << val
					//		 << "\n";
					nodeList.push_back(Node(it->second.size()));
					// Add edges
					for (auto &var : it->second) {
						int lowerBound = (vars[var].assigned() ? 1 : 0);
						nodeList.back().edgeList.push_back(NormalEdge(var, lowerBound, 1,
																													c(i, var)));
					}
				}
			}

			// Insert S node and S->Val edges
			nodeList.push_back(Node(inputVals.size()));
			for (int i = 0; i < inputVals.size(); i++) {
				int val = inputVals[i];
				auto valNode = valToNode->find(val);
				if (valNode != valToNode->end()) {
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
		// Populate updatedEdges, so we know where we should update the old residual
		// graph later on
		void FlowGraph::updatePrunedValues(Int::IntView x, unsigned int xIndex, 
													  vector<EdgeNodes>& updatedEdges) {
			// Hold iterators to the values that we end up prunning, so we can also
			// remove them from valToVars
			vector< unordered_set<int>::iterator > prunedValues; 
			// Iterate all values for which there is an edge to X
			auto& values = varToVals.map.find(xIndex)->second;
			for (auto valueIt = values.begin(); valueIt != values.end(); valueIt++) {
				auto value = *valueIt;
				auto valueNode = valToNode->find(value)->second;
				NormalEdge* edge = getEdge(valueNode, xIndex);
				if (!edge->lowerBound && edge->upperBound == 1) {
					// If edge hasn't already been pruned or assigned
					if (!x.in(value)) {
						// Value has been pruned from variable X's domain, update graph
						edge->upperBound = 0;
						if (edge->flow == 1) {
							oldFlowIsFeasible = false;
						}
						updatedEdges.push_back({valueNode, xIndex});
						prunedValues.push_back(valueIt);
					}
					if (x.assigned() && x.val() == value) {
						// Variable has been assigned with a value, update graph
						edge->lowerBound = 1;
						if (edge->flow == 0) {
							oldFlowIsFeasible = false;
						}
						updatedEdges.push_back({valueNode, xIndex});
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

// Iterate through each edge that has flow, to find its total cost
		int FlowGraph::calculateFlowCost(LI& lii) {
			flowCost = 0;
			for (unsigned int i = totalVarNodes; i < sNode(); i++) {
				for (auto& edge: nodeList[i].edgeList) {
					if (edge.flow > 0) {
						flowCost += edge.cost;
						lii[edge.destNode] = (*nodeToVal)[i];
					}
				}
			}
			return flowCost;
		}

		void FlowGraph::print() const {
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

		void FlowGraph::printResidual() const {
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