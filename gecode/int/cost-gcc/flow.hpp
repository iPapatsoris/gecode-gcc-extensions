#ifndef H_FLOW
#define H_FLOW

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>

using namespace std;

class Flow {
	typedef unordered_set<unsigned int> IntSet;

	unordered_map<unsigned int, IntSet> valVarFlow;
	IntSet *varTFlow; // TODO: put outside, add/remove within initial sendFlow function and not in here 
										// (we don't want to update it from advisors)
	bool firstFlow;

	public:

	Flow() : varTFlow(new IntSet()), firstFlow(true) {}

	void addValVarFlow(unsigned int val, unsigned int var) {
		auto res = valVarFlow.find(val);
		if (res == valVarFlow.end()) {
			auto inserted = valVarFlow.insert(pair<unsigned int, IntSet>(val, unordered_set<unsigned int>()));
			inserted.first->second.insert(var);
		} else {
			assert(res->second.find(var) == res->second.end());
			res->second.insert(var);
		}
	}

	void removeValVarFlow(unsigned int val, unsigned int var) {
		auto res = valVarFlow.find(val);
		assert(res != valVarFlow.end());
		int erased = res->second.erase(var);
		assert(erased);	
	}

	void addVarTFlow(unsigned var) {
		assert(firstFlow);
		assert(varTFlow->find(var) == varTFlow->end());
		varTFlow->insert(var);
	}

	void removeVarTFlow(unsigned var) {
		assert(firstFlow);
		int erased = varTFlow->erase(var);
		assert(erased);
	}

	unsigned int getValVarFlow(unsigned int val, unsigned int var) const {
		auto res = valVarFlow.find(val);
		if (res == valVarFlow.end()) {
			return 0;
		}
		return res->second.find(var) != res->second.end() ? 1 : 0; 
	}

	unsigned int getSValFlow(unsigned int val) const {
		auto res = valVarFlow.find(val);
		if (res == valVarFlow.end()) {
			return 0;
		}
		return res->second.size();
	}

	unsigned int getVarTFlow(unsigned int var) const {
		if (firstFlow) {
			return varTFlow->find(var) != varTFlow->end() ? 1 : 0;
		}

		for (auto& val: valVarFlow) {
			if (val.second.find(var) != val.second.end()) {
				return 1;
			}
		}
		return 0;
	}

	unsigned int getTSFlow() const {
		if (firstFlow) {
			return varTFlow->size();
		}
		unsigned int total = 0;
		for (auto& val: valVarFlow) {
			total += val.second.size();
		}
		return total;
	}

	void setFirstFlow(bool flag) {
		firstFlow = flag;
	}
};

#endif
