#include "flow-graph.hpp"
FlowGraph::FlowGraph(
	const ViewArray<Int::IntView>& vars, 
	const vector<unordered_set<int> >& varToVals,
	const MapToSet& valToVars,
	const IntArgs& inputVals, const IntArgs& lowerBounds, 
	const IntArgs& upperBounds, const IntArgs& costs) 
		: firstTimeValidCost(true) {
	
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
		NormalEdge e(tNode, 1, 1, 0);
		nodeList.back().edgeList.pushVal(e, &edgeListSize[x]);
	}

	Matrix<IntArgs> c(costs, inputVals.size(), vars.size());

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
				NormalEdge e(var, 0, 1, c(i, var));
				nodeList.back().edgeList.pushVal(e, &edgeListSize[valNode]);
			}
		}
	}

	// Insert S node and S->Val edges
	nodeList.push_back(Node(inputVals.size()));
	for (int i = 0; i < inputVals.size(); i++) {
		int val = inputVals[i];
		auto valNode = valToNode.find(val);
		NormalEdge e(valNode->second, lowerBounds[i], upperBounds[i], 0);
		nodeList.back().edgeList.pushVal(e, &edgeListSize[sNode]);
	}

	// Insert T node and T->S edge
	nodeList.push_back(Node(1));
	NormalEdge e(sNode, totalVarNodes, totalVarNodes, 0);
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

// Update graph state to match variable X domain prunings.
// If a variable-value pair is pruned that has no flow, delete it on the 
// spot. If it has flow, insert it on updatedEdges, for the flow repair
// algorithm to fix it on the next propagation, and mark flow as infeasible.
// If a variable is assigned but has no flow, mark flow as infeasible. 
// Returns whether the current flow is still feasible or not.
bool FlowGraph::updatePrunedValues(Int::IntView x, int xIndex, 
												vector<EdgeInfo>& updatedEdges) {
	// Hold the values that we end up prunning, so we can remove them from 
	// valToVars after iteration is done
	vector<int> prunedValues; 
	// Iterate all values for which there is an edge to X
	auto& values = backtrackStable->varToVals[xIndex];
	bool isFeasible = true;
	for (int i = 0; i < varToValsSize[xIndex]; i++) {
		auto value = (values.list)[i];
		auto valueNode = backtrackStable->valToNode.find(value)->second;
		NormalEdge* edge = getEdge(valueNode, xIndex);
		assert(edge != NULL);
		if (!x.in(value)) {
			// Value has been pruned from variable X's domain
			if (edge->flow == 1) {
				// Mark infeasible flow and edge to be repaired
				isFeasible = false;
				updatedEdges.push_back(EdgeInfo(valueNode, xIndex));
			} else {
				// No flow through the edge, can delete on the spot
				deleteEdge(valueNode, xIndex);
				prunedValues.push_back(value);
			}
		}
		if (x.assigned() && x.val() == value) {
			if (edge->flow == 0) {
				// Lower bound violation
				isFeasible = false;
			}
		}	
	}

	for (auto val: prunedValues) {
		values.deleteVal(val, &varToValsSize[xIndex]);
	}
	return isFeasible;
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