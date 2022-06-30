#ifndef H_RESIDUAL_GRAPH
#define H_RESIDUAL_GRAPH

#include <iostream>
#include <vector>
#include <memory>
#include <climits>

using namespace std;

class ResidualEdge2 {
	int active;
	int next;
	int prev;
	
	int cost;
	int reducedCost;
	int upperBound;

	friend class ResidualGraph;

	public:
	ResidualEdge2() : active(INT_MIN), next(-1), prev(-1), cost(0), 
							      reducedCost(0), upperBound(0) {}
};

class ResidualGraph {
	struct Bounds {
		int active;
		int start;
		int end;

		Bounds() : active(INT_MIN), start(-1), end(-1) {}
	};

	vector<ResidualEdge2> list;
	vector<Bounds> bounds;
	int activeFlag;
	int dimensionSize;

	public:
	ResidualGraph(int dimensionSize) : activeFlag(INT_MIN + 1), dimensionSize(dimensionSize) {
		list.assign(dimensionSize * dimensionSize, ResidualEdge2());
		bounds.assign(dimensionSize, Bounds());
	}

	void clear() {
		activeFlag++;
	}

	ResidualEdge2 *getResidualEdge(int src, int dest) {
		auto& res = list[dimensionSize * src + dest];
		return res.active >= activeFlag ? &res : NULL;
	}

	void addResidualEdge(int src, int dest, int cost, int reducedCost, int upperBound) {
		// cout << "adding " << src << "->" << dest << endl;
		int pos = dimensionSize * src + dest;
		auto& edge = list[pos];
		edge.active = activeFlag;
		edge.next = -1;
		if (bounds[src].active >= activeFlag) {
			int lastPos = bounds[src].end;
			bounds[src].end = pos;
			list[lastPos].next = pos; 
			edge.prev = lastPos;
		} else {
			bounds[src].active = activeFlag;
			bounds[src].start = bounds[src].end = pos;
			edge.prev = -1;
		}
		
		edge.cost = cost;
		edge.reducedCost = reducedCost;
		edge.upperBound = upperBound;
	}

	void deleteResidualEdge(int src, int dest) {
		// cout << "removing " << src << "->" << dest << endl;
		int pos = dimensionSize * src + dest;
		auto& edge = list[pos];
		edge.active--;
		if (bounds[src].start == pos && bounds[src].end == pos) {
			bounds[src].active--;
		} else if (bounds[src].start == pos) {
			bounds[src].start = edge.next;
			list[edge.next].prev = -1;
		} else if (bounds[src].end == pos) {
			bounds[src].end = edge.prev;
			list[edge.prev].next = -1;
		} else {
			list[edge.next].prev = edge.prev;
			list[edge.prev].next = edge.next;
		}
	}

	void addOrUpdate(int src, int dest, int cost, int reducedCost, int upperBound) {
		auto res = getResidualEdge(src, dest);
		if (res == NULL) {
			addResidualEdge(src, dest, cost, reducedCost, upperBound);
		} else {
			res->cost = cost;
			res->reducedCost = reducedCost;
			res->upperBound = upperBound;
		}
	}

	void markEdgesAlwaysActive(int src) {
		bounds[src].active = INT_MAX;
		int pos =  bounds[src].start;
		ResidualEdge2 *edge;
		do {
			edge = &list[pos];
			edge->active = INT_MAX;
			pos = edge->next;
		} while (pos != -1);
	}

	void calculateReducedCosts(const vector<int>& distances) {
		for (int i = 0; i < dimensionSize; i++) {
			if (bounds[i].active < activeFlag) {
				continue;
			}
			int pos = bounds[i].start;
			ResidualEdge2 *edge;
			do {
				edge = &list[pos];
				edge->reducedCost = distances[i] + edge->cost - distances[pos % dimensionSize];
				pos = edge->next;
			} while (pos != -1);
		}
	}


	void print() {
		cout << "New:\n";
		for (int i = 0; i < dimensionSize; i++) {
			if (bounds[i].active < activeFlag) {
				continue;
			}
			int z = bounds[i].start;
			ResidualEdge2 *edge;
			do {
				edge = &list[z];
				cout << i << " -> " << z % dimensionSize << " upper " << edge->upperBound << " reduced cost " << edge->reducedCost << " cost " << edge->cost << endl;
				z = edge->next;
			} while (z != -1);
		}
	}
};

#endif