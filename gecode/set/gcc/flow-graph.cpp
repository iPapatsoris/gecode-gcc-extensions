#include "flow-graph.hpp"
#include "example/sym-gcc-example.hpp"
#include <climits>

#define INF_UINT UINT_MAX

FlowGraph::FlowGraph( 
	const ViewArray<Set::SetView>& vars, 
	const vector<unordered_set<int> >& varToVals,
	const MapToSet& valToVars,
	const IntArgs& inputVals, const IntArgs& lowerValBounds, 
	const IntArgs& upperValBounds, const IntArgs& lowerVarBounds, 
	const IntArgs& upperVarBounds) {
	
	totalVarNodes = vars.size();
	int totalValNodes = inputVals.size();
	// Nodes are variable nodes, values nodes, S and T nodes
	int totalNodes = totalVarNodes + totalValNodes + 2;
	// S node position
	int sNode = totalNodes - 2;
	// T node position
	int tNode = totalNodes - 1;

	edgeListSize.assign(totalNodes, 0);
	varToValsSize.assign(totalVarNodes, 0);
	backtrackStable = make_shared<BacktrackStableContent>();
	auto& nodeList = backtrackStable->nodeList;
	auto& valToNode = backtrackStable->valToNode;
	auto& nodeToVal = backtrackStable->nodeToVal;
	nodeList.reserve(totalNodes);

	for (unsigned int var = 0; var < varToVals.size(); var++) {
		this->backtrackStable->varToVals.push_back(BtVector<int>(
																							 varToVals[var].size()));
		for (auto val: varToVals[var]) {
			this->backtrackStable->varToVals[var].pushVal(val, &varToValsSize[var]);
		}
	}

	// Insert variable nodes and var->T edges
	for (int x = 0; x < totalVarNodes; x++) {
		nodeList.push_back(Node(1));
		NormalEdge e(tNode, lowerVarBounds[x], upperVarBounds[x]);
		nodeList.back().edgeList.pushVal(e, &edgeListSize[x]);
	}

	// Insert Value nodes and Val->Var edges
	// It is important to iterate through inputVals and not through the
	// domains or valToVars, because values might have been early pruned from
	// the latter because of another constraint. We still need to include 
	// pruned values, to respect their lower bound restriction.
	for (int i = 0; i < inputVals.size(); i++) {
		int val = inputVals[i];
		auto it = valToVars.find(val);
		int valNode = nodeList.size();
		valToNode.insert({val, valNode});
		nodeToVal.insert({valNode, val});
		bool valIsPruned = (it == valToVars.end());
		nodeList.push_back(Node(valIsPruned ? 0 : it->second.size()));
		if (!valIsPruned) {
			// Add edges only if value is not already pruned
			for (auto &var : it->second) {
				NormalEdge e(var, 0, 1);
				nodeList.back().edgeList.pushVal(e, &edgeListSize[valNode]);
			}
		}
	}

	// Insert S node and S->Val edges
	nodeList.push_back(Node(inputVals.size()));
	for (int i = 0; i < inputVals.size(); i++) {
		int val = inputVals[i];
		auto valNode = valToNode.find(val);
		NormalEdge e(valNode->second, lowerValBounds[i], upperValBounds[i]);
		nodeList.back().edgeList.pushVal(e, &edgeListSize[sNode]);
	}

	// Insert T node and T->S edge
	nodeList.push_back(Node(1));
	NormalEdge e(sNode, 0, INT_MAX);
	nodeList.back().edgeList.pushVal(e, &edgeListSize[tNode]);

																		
	// Create residual graph
	for (unsigned int i = 0; i < nodeList.size(); i++) {
		auto& node = nodeList[i];
		auto& edges = node.edgeList;
		for (int e = 0; e < edgeListSize[i]; e++) {
			auto& edge = (edges.list)[e];
			node.residualEdgeList.push_back(ResidualEdge(edge));
		}
	}
}


//	TODO: update this comment

// Update graph state to match variable X domain pruning/assignment.
// Update is made by tightening the bounds of edge V->X as follows:
// - If X got assigned to value V, set the lower bound to 1.
// - For every value V that has been pruned off X, set the upper bound 
//   to 0. If we prune a value that is used by current flow, or assign a value 
// that is not used by it, set oldFlowIsFeasible to false.
// Populate updatedEdges, so we know where we should update the old residual
// graph later on
bool FlowGraph::updatePrunedValues(Set::SetView x, int xIndex, 
																   vector<EdgeInfo>& updatedEdges) {
	// Hold the values that we end up prunning, so we can remove them from 
	// valToVars after iteration is done
	vector<int> prunedValues; 
	// Iterate all values for which there is an edge to X
	auto& values = backtrackStable->varToVals[xIndex];
	bool isFeasible = true;
	// cout << "var " << xIndex << " lower bound: ";
	// for (SetVarGlbValues i(x); i(); ++i)
	// 	std::cout << i.val() << " ";
	// cout << "\nvar " << xIndex << " upper bound: ";
	// for (SetVarLubValues i(x); i(); ++i)
	// 	std::cout << i.val() << " ";
	// cout << "cardinalities " << x.cardMin() << " " << x.cardMax() << endl;
	// print();
	for (int i = 0; i < varToValsSize[xIndex]; i++) {
		auto value = (values.list)[i];
		auto valueNode = backtrackStable->valToNode.find(value)->second;
		NormalEdge* edge = getEdge(valueNode, xIndex);
		assert(edge != NULL);
		if (x.notContains(value)) {
			// cout << "val " << value << " (node " << valueNode << " no longer in var " << xIndex << " upper bound" << endl;
			// Value has been pruned from variable X's domain
			if (edge->flow == 1) {
				// cout << "upper bound violation" << endl;
				// Mark infeasible flow and edge to be repaired
				isFeasible = false;
				updatedEdges.push_back(EdgeInfo(valueNode, xIndex, false, 0));
			} else {
				// cout << "delete on the spot" << endl;
				// No flow through the edge, can delete on the spot
				deleteEdge(valueNode, xIndex);
				prunedValues.push_back(value);
			}
		}
		if (x.contains(value) && !edge->flow) {
			// cout << "val " << value << " (node " << valueNode << " in var " << xIndex << " lower bound but with no flow" << endl;
			// cout << "lower bound violation" << endl;
			isFeasible = false;
			updatedEdges.push_back(EdgeInfo(valueNode, xIndex, true, 1));
		}
	}

	for (auto val: prunedValues) {
		values.deleteVal(val, &varToValsSize[xIndex]);
	}

	// Update Var->T edge bounds
	auto& nodeList = backtrackStable->nodeList;
	auto edge = nodeList[xIndex].edgeList.getVal(tNode(), edgeListSize[xIndex]);
	if (edge->flow < (int) x.cardMin()) { 
		// cout << "flow " << edge->flow << " vs cardMin " << x.cardMin() << endl;
		isFeasible = false;
		updatedEdges.push_back(EdgeInfo(xIndex, tNode(), true, x.cardMin()));
	} else if (edge->flow > (int) x.cardMax()) {
		// cout << "flow " << edge->flow << " vs cardMax " << x.cardMax() << endl;
		isFeasible = false;
		updatedEdges.push_back(EdgeInfo(xIndex, tNode(), false, x.cardMax()));
	}

	return isFeasible;																		 
	

	// for (auto valIt = varToLub[xIndex].begin(); valIt != varToLub[xIndex].end(); valIt++) {
	// 	auto val = *valIt;
	// 	if (x.notContains(*valIt)) {
	// 		// Value has been pruned from variable's X's domain
	// 		auto valNode = valToNode->find(val)->second;
	// 		NormalEdge* edge = getEdge(valNode, xIndex);
	// 		if (!edge->lowerBound && edge->upperBound == 1) {
	// 			// If edge hasn't already been pruned or assigned, update graph
	// 			edge->upperBound = 0;
	// 			if (edge->flow == 1) {
	// 				oldFlowIsFeasible = false;
	// 				// if (bestBranch != NULL) {
	// 				// 	(*bestBranch)[xIndex].erase(val);
	// 				// }
	// 			}
	// 			updatedEdges.push_back({valNode, xIndex});
	// 			prunedValues.push_back(valIt);
	// 		}
	// 	}
	// }	
	
	// for (auto it: prunedValues) {
	// 	//cout << "\npruning from lub " << *it;
	// 	varToLub[xIndex].erase(it);
	// }

	// vector<int> valuesToAdd;
	// for (SetVarGlbValues i(x); i(); ++i) {
	// 	auto val = i.val();
	// 	if (varToGlb[xIndex].find(val) == varToGlb[xIndex].end()) {
	// 		// Value has been included in variable X's domain
	// 		auto valNode = valToNode->find(val)->second;
	// 		NormalEdge* edge = getEdge(valNode, xIndex);
	// 		if (!edge->lowerBound && edge->upperBound == 1) {
	// 			// If edge hasn't already been pruned or assigned, update graph
	// 			edge->lowerBound = 1;
	// 			if (!edge->flow) {
	// 				oldFlowIsFeasible = false;
	// 			} /*else if (bestBranch != NULL) {
	// 				for (SetVarUnknownValues nextVal(x); nextVal(); ++nextVal) {
	// 					auto nextValNode = (*valToNode)[nextVal.val()];
	// 					auto nextValEdge = getEdge(nextValNode, xIndex);
	// 					cout << "checking val " << nextVal.val(); 
	// 					if (nextValEdge->flow) {
	// 						cout << "ok";
	// 						(*bestBranch)[xIndex] = val;
	// 					}
	// 				}
	// 			}*/
	// 			updatedEdges.push_back({valNode, xIndex});
	// 			valuesToAdd.push_back(val);
	// 		}
	// 	}
	// }	
	
	// for (auto val: valuesToAdd) {
	// 	//cout << "\nadding to glb " << val;
	// 	varToGlb[xIndex].insert(val);
	// }

	// Update Var->T edge bounds
	// auto& edge = nodeList[xIndex].edgeList[0];
	// edge.lowerBound = x.cardMin();
	// edge.upperBound = x.cardMax();
	// if (edge.flow < edge.lowerBound || edge.flow > edge.upperBound) {
	// 	oldFlowIsFeasible = false;
	// 	updatedEdges.push_back({xIndex, tNode()});
	// }
	//cout << "\ncard " << edge.lowerBound << " " << edge.upperBound << "\n";

	//cout << "end of update: feasible flow " << oldFlowIsFeasible << endl; 
}

void FlowGraph::print() const {
	for (unsigned int i = 0; i < backtrackStable->nodeList.size(); i++) {
		auto& node = backtrackStable->nodeList[i];
		auto& edges = node.edgeList;
		for (int j = 0; j < edgeListSize[i]; j++) {
			auto& edge = (edges.list)[j];
			cout << i << " -> ";
			edge.print();
		}
	}
	cout << endl;
}

void FlowGraph::printResidual() const {
	for (unsigned int i = 0; i < backtrackStable->nodeList.size(); i++) {
		auto& node = backtrackStable->nodeList[i];
		for (auto& edge: node.residualEdgeList) {
			cout << i << " -> ";
			edge.print();
		}
	}
	cout << endl;
}

// 	cout << "glb\n";
// 	for (auto v: varToGlb[x]) {
// 		cout << v << " ";
// 	}
// 	cout << "\nlub\n";
// 	for (auto v: varToLub[x]) {
// 		cout << v << " ";
// 	}
// 	cout << "\ncard " << nodeList[x].edgeList.front().lowerBound
// 			 << " " << nodeList[x].edgeList.front().upperBound << "\n";
// }
