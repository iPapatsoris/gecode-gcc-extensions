#include <vector>
#include <unordered_map>

using namespace std;

class BtVector {
	vector<int>* list;
	unordered_map<int, unsigned int> *valToPos;
	unsigned int listSize; 

	public:

		BtVector(unsigned int totalEdges) : listSize(totalEdges) {
			list = new vector<int>();
			valToPos = new unordered_map<int, unsigned int>();
			valToPos->reserve(totalEdges);
			list->reserve(totalEdges);
		}

		friend class FlowGraph;
		friend class FlowGraphAlgorithms;

		int* getVal(int val) {
				auto res = valToPos->find(val);
				if (res == valToPos->end() || res->second >= listSize) {
					return NULL;
				}
				return &(*list)[res->second];
			}

		void deleteVal(int val) {
			auto res = valToPos->find(val);
			if (res == valToPos->end() || res->second >= listSize) {
				return;
			}

			auto pos = res->second;
			auto lastElement = (*list)[listSize - 1];
			swap((*list)[pos], (*list)[listSize - 1]);
			
			(*valToPos)[val] = listSize - 1;
			(*valToPos)[lastElement] = pos;
			listSize -= 1;
		}

		void pushVal(int val) {
			list->push_back(val);
			valToPos->insert(pair<unsigned int, unsigned int>(val, list->size() - 1));
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
};
