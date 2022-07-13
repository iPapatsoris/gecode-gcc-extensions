#ifndef H_BT_VECTOR
#define H_BT_VECTOR

#include <vector>
#include <unordered_map>
#include "edge.hpp"

using namespace std;

/**
 * Specialized vector for backtracking. Suitable only for deletions.
 * Place structure within a shared_ptr so that it will not be copied on cloning.
 * For each instance created, there should be a "listSize" field that is placed 
 * outside the class, because it is the only one that needs to be copied. 
 * See "FlowGraph::edgeListSize" and "FlowGraph::varToValsSize".
 * On deletion, the element to be deleted is swapped with the last one,
 * and listSize is decremented by 1. The element is still within the vector.
 * When backtracking occurs, the previous value of listSize is restored,
 * and thus the deleted elements can be accessed again with minimal overhead.
 * 
 * The data we hold can be any type T, but the map keys for the elements are 
 * type int and not T, because in the latter case we would need to create a 
 * dummy object of type T every time we want to search for the existence of an 
 * element, and we would also have more overhead for the map keys for no reason. 
 * Instead, implement BtVector::getValID() for each different possible type, to 
 * convert them to a suitable int key.
 * 
 * During iteration, client MUST iterate taking into account listSize and not 
 * the vector's actual internal size. The class is exposed to its friend classes 
 * for convenient and performant iteration. Alternatively for more encapsulation
 * custom iterators can be implemented.
 */
template <typename T>
class BtVector {
	vector<T> list;
	unordered_map<int, int> valToPos;

	friend class FlowGraph;
	friend class FlowGraphAlgorithms;

	public:
		BtVector(int totalEdges) {
			valToPos.reserve(totalEdges);
			list.reserve(totalEdges);
		}

		// Search val in map, return the actual element
		T* getVal(int val, int listSize) {
				int pos;
				if (!getValPos(val, &pos, listSize)) {
					return NULL;
				}
				return &(list)[pos];
			}
			
		// Swap element to be deleted with the last one, decrement listSize
		void deleteVal(int val, int *listSize) {
			int pos;
			if (!getValPos(val, &pos, *listSize)) {
				return;
			}
			auto lastElement = getValID((list)[*listSize - 1]);
			deleteValAtPos(val, pos, lastElement, listSize);
		}

		// Values should NOT be added after search begins
		void pushVal(T val, int *listSize) {
			list.push_back(val);
			valToPos.insert(pair<int, int>(getValID(val), list.size() - 1));
			(*listSize)++;
		}

		void print(int listSize) const {
			cout << "size: " << listSize << "\n";
			cout << "valToPos:";
			for (auto& p: valToPos) {
				cout << p.first << ": " << p.second << " ";
			}
			cout << "list:\n";
			for (auto& e: list) {
				cout << e << " ";
			}
			cout << endl;
		}
	
	private:
		// Search val in map, return true and the postion in list if found,
		// false otherwise 
		bool getValPos(int val, int *pos, int listSize) const {
			auto res = valToPos.find(val);
			if (res == valToPos.end() || res->second >= listSize) {
				return false;
			}

			*pos = res->second;
			return true;
		}

		// Swap element to be deleted with the last one, decrement listSize
		void deleteValAtPos(int val, int pos, 
												int lastElement, int *listSize)
		{
			swap((list)[pos], (list)[*listSize - 1]);
			
			(valToPos)[val] = *listSize - 1;
			(valToPos)[lastElement] = pos;
			*listSize -= 1;
		}

		// Get identification for map
		int getValID(NormalEdge val) const {
			return val.getDestNode();
		}

		// Get identification for int
		int getValID(int val) const {
			return val;
		}

};

#endif