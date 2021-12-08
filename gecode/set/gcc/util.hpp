#ifndef H_SYM_UTIL
#define H_SYM_UTIL

#include <unordered_map>
#include <unordered_set>
#include <gecode/int.hh>

using namespace std;

template <class T1, class T2>
struct MapToSet {
	unordered_map<T1, unordered_set<T2>> map;
	MapToSet() {}
	MapToSet(const MapToSet &c) {
		map = c.map;
	}
};

struct VarUtil {
	struct InputVarInfo {
		unsigned int xIndex;
		unordered_map<int, unsigned int> valToValIndex; // in xToInputVar
	};

	struct XInfo {
		unsigned int varIndex;
		int val;
	};

	vector<InputVarInfo> inputVarToXIndex;
	vector<XInfo> xToInputVar;

	unsigned int getInputVarFromXIndex(unsigned int xIndex) const {
		return xToInputVar[xIndex].varIndex;
	}

	int getValFromXIndex(unsigned int xIndex) const {
		return xToInputVar[xIndex].val;
	}

	unsigned int getXFromInputVarVal(unsigned int inputVar, int val) const {
		unsigned int xIndex = inputVarToXIndex[inputVar].xIndex; 
		return xIndex + inputVarToXIndex[inputVar].valToValIndex.find(val)->second;
	}

	bool inputVarIsAssigned(unsigned int inputVarIndex, const ViewArray<Int::BoolView>& x) const {
		auto& inputVar = inputVarToXIndex[inputVarIndex];
		auto xIndex = inputVar.xIndex;
		for (; xIndex <= xIndex + inputVar.valToValIndex.size(); xIndex++) {
			if (!x[xIndex].assigned()) {
				return false;
			}
		}
		return true;
	}

	/*bool inputVarIsAssignedToVal(unsigned int inputVarIndex, int val, const ViewArray<Int::BoolView>& x) const {
		auto& inputVar = inputVarToXIndex[inputVarIndex];
		auto xIndex = inputVar.xIndex;
		for (; xIndex <= xIndex + inputVar.valToValIndex.size(); xIndex++) {
			if (!x[xIndex].assigned()) {
				return false;
			}
		}
		return true;
	}*/

	void print() const {
		cout << "inputVarToXIndex";
		for (unsigned int inputVar = 0; inputVar < inputVarToXIndex.size(); inputVar++) {
			auto& inputVarInfo = inputVarToXIndex[inputVar];
			cout << "\n\nInput var " << inputVar << " at index " << inputVarInfo.xIndex << endl;
			for (auto& valInfo: inputVarInfo.valToValIndex) {
				cout << "val " << valInfo.first << " at local pos " << valInfo.second << "\n";
			}
		}

		cout << "\nxToInputVar";
		for (unsigned int x = 0; x < xToInputVar.size(); x++) {
			auto& xInfo = xToInputVar[x];
			cout << "x[" << x << "] maps to inputVar " << xInfo.varIndex << ", val " << xInfo.val << endl;
		}
	}
};

#endif