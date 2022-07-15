#ifndef H_GCC_EXAMPLE
#define H_GCC_EXAMPLE

#include <gecode/driver.hh>
#include <gecode/set.hh>
#include <gecode/minimodel.hh>
//#include "BestBranch.hpp"
#include "read-input.hpp"

using namespace Gecode;
using namespace std;


class SymmetricGccExample : public Script {
protected:
	SetVarArray x;
	IntVarArray y;
	BoolVarArray z;
	int varsCount;
	int valsCount;


public:
	enum Model {
		MODEL_SINGLE, MODEL_COUNT, MODEL_LINEAR 
	};

	enum Branching {
		BRANCHING_SIMPLE, BRANCHING_CUSTOM
	};

	Model model;
	SymmetricGccExample(const InstanceOptions& opt);
	SymmetricGccExample(SymmetricGccExample &s) : Script(s) {
		model = s.model;
		varsCount = s.varsCount;
		valsCount = s.valsCount;
		switch (model) {
			case MODEL_SINGLE:
				x.update(*this, s.x);
				break;
			case MODEL_COUNT:
				y.update(*this, s.y);
				break;
			case MODEL_LINEAR:
				z.update(*this, s.z);
				break;
			}
	}
	virtual Space *copy(void) {
		return new SymmetricGccExample(*this);
	}
	void print(ostream& os) const {
		return;
		switch (model) {
			case MODEL_SINGLE:
				os << "\tSolution: " << x << "\n";
				break;
			case MODEL_COUNT:
			case MODEL_LINEAR:
				os << "\tSolution: \n";
				for (int i = 0; i < varsCount; i++) {
					for (int j = 0; j < valsCount; j++) {
						if (model == MODEL_COUNT) {
							cout << y[valsCount*i + j] << " ";
						} else if (model == MODEL_LINEAR) {
							cout << z[valsCount*i + j] << " ";
						}
					}
					cout << "\n";
				}
				break;
		}
	}
};

#endif