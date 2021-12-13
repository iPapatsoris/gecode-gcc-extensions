#include "flow-graph.hpp"
#include "example/sym-gcc-example.hpp"
#include <climits>

#define INF_UINT UINT_MAX

FlowGraph::FlowGraph( 
	const ViewArray<Int::BoolView>& vars, 
	const MapToSet<int, unsigned int>& valToVars, const vector<unordered_set<int> >& varToVals,
	const IntArgs& inputVals, const IntArgs& lowerValBounds, 
	const IntArgs& upperValBounds, const IntArgs& lowerVarBounds, 
	const IntArgs& upperVarBounds, const VarUtil& varUtil) 
		: oldFlowIsFeasible(true), varToVals(varToVals), varUtil(varUtil) {



	totalVarNodes = varToVals.size();
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
		nodeList.back().edgeList.push_back(NormalEdge(tNode, lowerVarBounds[x], 
																			 upperVarBounds[x]));
	}

	// Insert Value nodes and Val->Var edges
	// It is important to iterate through inputVals and not through the
	// domains or valToVars, because values might have been early pruned from
	// the latter (from a different contraint). We still need to include pruned 
	// values, to respect their lower bound restriction
	for (int i = 0; i < inputVals.size(); i++) {
		int val = inputVals[i];
		auto it = valToVars.map.find(val);
			valToNode->insert({val, nodeList.size()});
			nodeToVal->insert({nodeList.size(), val});
			//cout << "node " << nodeList.size() << " corresponds to val " << val
			//		 << "\n";
			bool valIsPruned = (it == valToVars.map.end());
			nodeList.push_back(Node(valIsPruned ? 0 : it->second.size()));
		if (!valIsPruned) {
			// Add edges
			for (auto &var : it->second) {				
				int lowerBound = (varUtil.inputVarIsAssigned(var, vars) ? 1 : 0);
				nodeList.back().edgeList.push_back(NormalEdge(var, lowerBound, 1));
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
																				 lowerValBounds[i], upperValBounds[i]));
		}
	}

	// Insert T node and T->S edge
	nodeList.push_back(Node(1));
	nodeList.back().edgeList.push_back(NormalEdge(sNode, 0, INF_UINT));

	// Create residual graph
	for (auto &node : nodeList) {
		copy(node.edgeList.begin(), node.edgeList.end(), 
				 back_inserter(node.residualEdgeList));
	}
}

//	TODO: update this comment

// Update graph state to match variable X domain pruning/assignment.
// Update is made by tightening the bounds of edge V->X as follows:
// - If X got assigned to value V, set the lower bound to 1.
// - For every value V that has been pruned off X, set the upper bound 
//   to 0. 
// If we prune a value that is used by current flow, or assign a value 
// that is not used by it, set oldFlowIsFeasible to false.
// Populate updatedEdges, so we know where we should update the old residual
// graph later on
void FlowGraph::updatePrunedValues(Int::BoolView x, unsigned int xIndex, 
																   vector<EdgeNodes>& updatedEdges, LI* li) {
	unsigned int xNodeIndex = varUtil.getInputVarFromXIndex(xIndex);		
	int val = varUtil.getValFromXIndex(xIndex);						
	
	auto valNode = valToNode->find(val)->second;
	NormalEdge* edge = getEdge(valNode, xNodeIndex);
	if (edge->lowerBound == 1 || !edge->upperBound) {
		// Edge has already been updated
		//cout << "end of update: feasible flow " << oldFlowIsFeasible << endl; 
		return;
	}

	updatedEdges.push_back({valNode, xNodeIndex});
	if (!x.val()) {
		// Value has been pruned from variable X domain
		edge->upperBound = 0;
		if (edge->flow == 1) {
			oldFlowIsFeasible = false;
		}
	} else if (x.val() == 1) {
		// Value has been assigned to variable X domain
		edge->lowerBound = 1;
		if (!edge->flow) {
			oldFlowIsFeasible = false;
		}
	}

	// Update Var->T edge bounds
	/*auto& edge = nodeList[xNodeIndex].edgeList[0];
	edge.lowerBound = x.cardMin();
	edge.upperBound = x.cardMax();
	if (edge.flow < edge.lowerBound || edge.flow > edge.upperBound) {
		oldFlowIsFeasible = false;
		updatedEdges.push_back({xNodeIndex, tNode()});
	}
	cout << "\ncard " << edge.lowerBound << " " << edge.upperBound << "\n";*/

	//cout << "end of update: feasible flow " << oldFlowIsFeasible << endl; 
}

void FlowGraph::print() const {

	for (unsigned int inputVarIndex = 0; inputVarIndex < varToVals.size(); inputVarIndex++) {
		cout << "Input var " << inputVarIndex << " has vals:\n";
		for (auto val: varToVals[inputVarIndex]) {
			cout << val << " which maps to node " << (*valToNode)[val] << "\n";
		}
		cout << "\n";
	}

	for (auto& n: *nodeToVal) {
		cout << "Node " << n.first << " maps to val " << n.second << "\n";
	}
	cout << "\n";

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

/*void FlowGraph::printBounds(int x) const {
	cout << "glb\n";
	for (auto v: varToGlb[x]) {
		cout << v << " ";
	}
	cout << "\nlub\n";
	for (auto v: varToLub[x]) {
		cout << v << " ";
	}
	cout << "\ncard " << nodeList[x].edgeList.front().lowerBound
			 << " " << nodeList[x].edgeList.front().upperBound << "\n";
}
*/