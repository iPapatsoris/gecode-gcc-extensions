#ifndef H_BT_VECTOR
#define H_BT_VECTOR

#include <vector>
#include <unordered_map>
#include "edge.hpp"

using namespace std;

/**
 * Specialized vector for backtracking. Suitable only for deletions.
 * List and valToPos are pointers so that they will not be copied on branching.
 * Only listSize is backtracked. 
 * On deletion, the element to be deleted is swapped with the last one,
 * and listSize is decremented by 1. The element is still within the vector.
 * When backtracking occurs, the previous value of listSize is restored,
 * and thus the deleted elements can be accessed again with minimal overhead.
 * 
 * During iteration, client MUST iterate taking into account listSize and not 
 * the vector's actual internal size. The class is exposed to its friend classes 
 * for easy and performant iteration. Alternatively for more encapsulation, 
 * custom iterators can be implemented.
 * 
 * typename Index is the id used to identify T in valToPos. Alternatively 
 * we could have specified a hash function for T, but in that case we 
 * would have to create a dummy T object every time we want to look for its
 * existence. Instead we use directly the ID for performance, thus the reason 
 * why we need to specify it as a separate type (FlowGraph::varToVals uses
 * int as Index, since a value can be negative, Node::edgeList uses 
 * unsigned int, since nodes are only positive).  
 */
template <typename T, typename Index>
class BtVector {
	vector<T>* list;
	unordered_map<Index, unsigned int> *valToPos;
	unsigned int listSize; 

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;

	public:
		BtVector(unsigned int totalEdges) : listSize(totalEdges) {
			list = new vector<T>();
			valToPos = new unordered_map<Index, unsigned int>();
			valToPos->reserve(totalEdges);
			list->reserve(totalEdges);
		}

		T* getVal(Index val) {
				unsigned int pos;
				if (!getValPos(val, &pos)) {
					return NULL;
				}
				return &(*list)[pos];
			}

		// Deletion for BTVector<NormalEdge, unsinged int> (Node::edgeList)
		void deleteVal(unsigned int val) {
			unsigned int pos;
			if (!getValPos(val, &pos)) {
				return;
			}
			auto lastElement = ((NormalEdge)(*list)[listSize - 1]).getDestNode();
			deleteValAtPos(val, pos, lastElement);
		}

		// Deletion for BTVector<int, int> (FlowGraph::varToVals)
		void deleteVal(int val) {
			unsigned int pos;
			if (!getValPos(val, &pos)) {
				return;
			}
			auto lastElement = (unsigned int)(*list)[listSize - 1];
			deleteValAtPos(val, pos, lastElement);
		}

		// Values should NOT be added after search begins
		void pushVal(NormalEdge val) {
			list->push_back(val);
			valToPos->insert(pair<Index, unsigned int>(val.getDestNode(), list->size() - 1));
		}

		// Values should NOT be added after search begins
		void pushVal(int val) {
			list->push_back(val);
			valToPos->insert(pair<Index, unsigned int>(val, list->size() - 1));
		}

		void print() const {
			cout << "size: " << listSize << "\n";
			cout << "valToPos:";
			for (auto& p: *valToPos) {
				cout << p.first << ": " << p.second << " ";
			}
			cout << "list:\n";
			for (auto& e: *list) {
				cout << e << " ";
			}
			cout << endl;
		}
	
	private:
		// Return true and the postion of val if found, false otherwise 
		bool getValPos(Index val, unsigned int *pos) const {
			auto res = valToPos->find(val);
			if (res == valToPos->end() || res->second >= listSize) {
				return false;
			}

			*pos = res->second;
			return true;
		}

		// Swap element to be deleted with the last one, decrement listSize
		void deleteValAtPos(Index val, unsigned int pos, 
												unsigned int lastElement)
		{
			swap((*list)[pos], (*list)[listSize - 1]);
			
			(*valToPos)[val] = listSize - 1;
			(*valToPos)[lastElement] = pos;
			listSize -= 1;
		}

};

#endif