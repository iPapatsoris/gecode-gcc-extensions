#include "flow-graph.hpp"
#include "example/sym-gcc-example.hpp"
#include <climits>

#define INF_UINT UINT_MAX

FlowGraph::FlowGraph( 
	const ViewArray<Set::SetView>& vars, 
	const MapToSet<int, unsigned int>& valToVars,
	const IntArgs& inputVals, const IntArgs& lowerValBounds, 
	const IntArgs& upperValBounds, const IntArgs& lowerVarBounds, 
	const IntArgs& upperVarBounds) 
		: oldFlowIsFeasible(true) {
	
	varToGlb.assign(vars.size(), unordered_set<int>());
	varToLub.assign(vars.size(), unordered_set<int>());
	for (int x = 0; x < vars.size(); x++) {
		for (SetVarGlbValues i(vars[x]); i(); ++i) {
			varToGlb[x].insert(i.val());
		}
		for (SetVarLubValues i(vars[x]); i(); ++i) {
			varToLub[x].insert(i.val());
		}
	}

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
				int lowerBound = (vars[var].assigned() ? 1 : 0);
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
void FlowGraph::updatePrunedValues(Set::SetView x, unsigned int xIndex, 
																   vector<EdgeNodes>& updatedEdges, LI* li) {
	// Hold iterators to the values that we end up prunning, so we can also
	// remove them from valToVars
	vector< unordered_set<int>::iterator > prunedValues; 
	for (auto valIt = varToLub[xIndex].begin(); valIt != varToLub[xIndex].end(); valIt++) {
		auto val = *valIt;
		if (x.notContains(*valIt)) {
			// Value has been pruned from variable's X's domain
			auto valNode = valToNode->find(val)->second;
			NormalEdge* edge = getEdge(valNode, xIndex);
			if (!edge->lowerBound && edge->upperBound == 1) {
				// If edge hasn't already been pruned or assigned, update graph
				edge->upperBound = 0;
				if (edge->flow == 1) {
					oldFlowIsFeasible = false;
					if (li != NULL) {
						(*li)[xIndex].erase(val);
					}
				}
				updatedEdges.push_back({valNode, xIndex});
				prunedValues.push_back(valIt);
			}
		}
	}	
	
	for (auto it: prunedValues) {
		//cout << "\npruning from lub " << *it;
		varToLub[xIndex].erase(it);
	}

	vector<int> valuesToAdd;
	for (SetVarGlbValues i(x); i(); ++i) {
		auto val = i.val();
		if (varToGlb[xIndex].find(val) == varToGlb[xIndex].end()) {
			// Value has been included in variable X's domain
			auto valNode = valToNode->find(val)->second;
			NormalEdge* edge = getEdge(valNode, xIndex);
			if (!edge->lowerBound && edge->upperBound == 1) {
				// If edge hasn't already been pruned or assigned, update graph
				edge->lowerBound = 1;
				if (!edge->flow) {
					oldFlowIsFeasible = false;
				} /*else if (li != NULL) {
					for (SetVarUnknownValues nextVal(x); nextVal(); ++nextVal) {
						auto nextValNode = (*valToNode)[nextVal.val()];
						auto nextValEdge = getEdge(nextValNode, xIndex);
						cout << "checking val " << nextVal.val(); 
						if (nextValEdge->flow) {
							cout << "ok";
							(*li)[xIndex] = val;
						}
					}
				}*/
				updatedEdges.push_back({valNode, xIndex});
				valuesToAdd.push_back(val);
			}
		}
	}	
	
	for (auto val: valuesToAdd) {
		//cout << "\nadding to glb " << val;
		varToGlb[xIndex].insert(val);
	}

	// Update Var->T edge bounds
	auto& edge = nodeList[xIndex].edgeList[0];
	edge.lowerBound = x.cardMin();
	edge.upperBound = x.cardMax();
	if (edge.flow < edge.lowerBound || edge.flow > edge.upperBound) {
		oldFlowIsFeasible = false;
		updatedEdges.push_back({xIndex, tNode()});
	}
	//cout << "\ncard " << edge.lowerBound << " " << edge.upperBound << "\n";

	//cout << "end of update: feasible flow " << oldFlowIsFeasible << endl; 
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

void FlowGraph::printBounds(int x) const {
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
