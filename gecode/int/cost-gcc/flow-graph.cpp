#include "flow-graph.hpp"
#include "example/cost-gcc-example.hpp"
unsigned long countNi = 0;
		FlowGraph::FlowGraph(
			const ViewArray<Int::IntView>& vars, 
 			const vector<unordered_set<int> >& varToVals,
			const MapToSet<int, unsigned int>& valToVars,
			const IntArgs& inputVals, const IntArgs& lowerBounds, 
			const IntArgs& upperBounds, const IntArgs& costs, int costUpperBound) 
				: costUpperBound(costUpperBound), firstTimeValidCost(true) {
			
			for (unsigned int var = 0; var < varToVals.size(); var++) {
				this->varToVals.push_back(BtVector(varToVals[var].size()));
				for (auto val: varToVals[var]) {
					this->varToVals[var].pushVal(val);
				}
			}

			flowCost = 0;
			oldFlowIsFeasible = new bool(true);
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
			//debug.reserve(totalNodes);

			// Insert variable nodes and var->T edges
			for (unsigned int x = 0; x < totalVarNodes; x++) {
				nodeList.push_back(Node(1));
				//debug.push_back(1);
				nodeList.back().edgeList->push_back(NormalEdge(tNode, 1, 1, 0));
				nodeList.back().edgeToPos->insert(pair<unsigned int, unsigned int>
																			    (tNode, 0));
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
				valToNode->insert({val, nodeList.size()});
				nodeToVal->insert({nodeList.size(), val});
				//cout << "node " << nodeList.size() << " corresponds to val " << val
				//		 << "\n";
				bool valIsPruned = (it == valToVars.map.end());
				nodeList.push_back(Node(valIsPruned ? 0 : it->second.size()));
				//debug.push_back(1);
				if (!valIsPruned) {
					// Add edges
					for (auto &var : it->second) {
						// int lowerBound = (vars[var].assigned() ? 1 : 0);
						nodeList.back().edgeList->push_back(NormalEdge(var, 0, 1,
																													c(i, var)));
						nodeList.back().edgeToPos->insert(pair<unsigned int, unsigned int>
																			    (var, nodeList.back().edgeList->size() - 1));
					}
				}
			}

			// Insert S node and S->Val edges
			nodeList.push_back(Node(inputVals.size()));
		//	debug.push_back(1);
			for (int i = 0; i < inputVals.size(); i++) {
				int val = inputVals[i];
				auto valNode = valToNode->find(val);
				nodeList.back().edgeList->push_back(NormalEdge(valNode->second, 
																				   lowerBounds[i], upperBounds[i], 0));
				nodeList.back().edgeToPos->insert(pair<unsigned int, unsigned int>
																			    (valNode->second, nodeList.back().edgeList->size() - 1));
			}

			// Insert T node and T->S edge
			nodeList.push_back(Node(1));
		//	debug.push_back(1);
			nodeList.back().edgeList->push_back(NormalEdge(sNode, totalVarNodes, 
																				 totalVarNodes, 0));
			nodeList.back().edgeToPos->insert(pair<unsigned int, unsigned int>
																			    (sNode, nodeList.back().edgeList->size() - 1));


			// Create residual graph
			/*for (auto &node : nodeList) {
				copy(node.edgeList->begin(), node.edgeList->end(), 
						 back_inserter(*node.residualEdgeList));
			}*/
			for (auto &node : nodeList) {
				for (unsigned int e = 0; e < node.edgeListSize; e++) {
					auto& edge = (*node.edgeList)[e];
					node.residualEdgeList->push_back(ResidualEdge(edge));
				}
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
		bool FlowGraph::updatePrunedValues(Int::IntView x, unsigned int xIndex, 
													  vector<EdgeUpdate>& updatedEdges) {
			// Hold iterators to the values that we end up prunning, so we can also
			// remove them from valToVars
			vector<int> prunedValues; 
			bool isFeasible = true;
			// Iterate all values for which there is an edge to X
			auto& values = varToVals[xIndex];
			for (unsigned int i = 0; i < values.listSize; i++) {
				auto value = (*values.list)[i];
				auto valueNode = valToNode->find(value)->second;
				NormalEdge* edge = getEdge(valueNode, xIndex);
				if (edge == NULL) {
				//	cout << "Looking for edge " << valueNode << "->" <<  xIndex << endl; 
				}
				if (!edge->lowerBound && edge->upperBound == 1) { // TODO: not needed
					// If edge hasn't already been pruned or assigned
					if (!x.in(value)) {
						// Value has been pruned from variable X's domain, update graph
						//edge->upperBound = 0;
						if (flow.getValVarFlow(valueNode, xIndex) == 1) {
							oldFlowIsFeasible = false;
							isFeasible = false;
							updatedEdges.push_back(EdgeUpdate(valueNode, xIndex));
//  cout << "upper bound violation " << valueNode << " " << xIndex << endl;							//*flowCost -= edge->cost;
							//edge->flow -= 1;
							//getEdge(xIndex, tNode())->flow -= 1;
							//getEdge(tNode(), sNode())->flow -= 1;
							//getEdge(sNode(), valueNode)->flow -= 1;
						} else {
														//  cout << "deleting " << valueNode << " " << xIndex << endl;

							deleteEdge(valueNode, xIndex);
							prunedValues.push_back(value);
							//deleteResidualEdge()
						}
//						cout << "ADDING " << valueNode << "->" << xIndex << endl;
						//prunedValues.push_back(value);
					}
					if (x.assigned() && x.val() == value) {
						// Variable has been assigned with a value, update graph
						//edge->lowerBound = 1;
						if (flow.getValVarFlow(valueNode, xIndex) == 0) {
							oldFlowIsFeasible = false;
							isFeasible = false;
														//  cout << "lower bound violation " << valueNode << " " << xIndex << endl;

						}
//						cout << "ADDING " << valueNode << "->" << xIndex << endl;
					}	
				}
			}

			for (auto val: prunedValues) {
			//	cout << "Deleting val " << val << " (node )" << (*valToNode)[val] << " from varToVals " << xIndex << endl;
				values.deleteVal(val);
			}
/*			cout << "Updated Edges:\n";
			for (auto& e: updatedEdges) {
				cout << e.src << "->" << e.dest << " " << e.lowerBoundViolation << " " << e.upperBoundViolation << " " << e.deleted << endl;
			}*/
			
			#ifndef NDEBUG
				//assertVarToValsInSync(x, xIndex);
			#endif
			return isFeasible;
		}

// Iterate through each edge that has flow, to find its total cost
	/*	int FlowGraph::calculateFlowCost(LI& lii) {
			*flowCost = 0;
			for (unsigned int i = totalVarNodes; i < sNode(); i++) {
				for (unsigned int j = 0; j < nodeList[i].edgeListSize; j++) {
					auto& edge = (*(nodeList[i].edgeList))[j];
					if (edge.flow > 0) {
						*flowCost += edge.cost;
						lii[edge.destNode] = (*nodeToVal)[i];
					}
				}
			}
			return *flowCost;
		}*/

		void FlowGraph::print() const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				auto& node = nodeList[i];
				//node.print();
				//if (i < totalVarNodes)
					//varToVals[i].print();
				for (unsigned int j = 0; j < nodeList[i].edgeListSize; j++) {
					auto& edge = (*(nodeList[i].edgeList))[j];
					cout << i << " -> ";
					edge.print(i, sNode(), tNode(), flow);
				}
			}
			cout << endl;
		}

		void FlowGraph::printResidual() const {
			for (unsigned int i = 0; i < nodeList.size(); i++) {
				auto& node = nodeList[i];
				//cout << "address " << node.residualEdgeList << endl;
				for (auto& edge: *node.residualEdgeList) {
					//cout << "capacity " << (*node.residualEdgeList) << " ";
					cout << i << " -> ";
					edge.print();
				}
			}
			cout << endl;
		}