#ifndef H_GCC_EXAMPLE
#define H_GCC_EXAMPLE

#include <gecode/driver.hh>
#include <gecode/set.hh>
#include <gecode/minimodel.hh>

using namespace Gecode;
using namespace std;

class SymmetricGccExample : public Script {
protected:
	SetVarArray x;
	IntVarArray y;
	int weeks;
	int periods;

public:
	enum Model {
		MODEL_SINGLE, MODEL_COUNT
	};

	Model model;
	SymmetricGccExample(const InstanceOptions& opt);
	SymmetricGccExample(SymmetricGccExample &s) : Script(s) {
		model = s.model;
		weeks = s.weeks;
		periods = s.periods;
		switch (model) {
			case MODEL_SINGLE:
				x.update(*this, s.x);
				break;
			case MODEL_COUNT:
				y.update(*this, s.y);
				break;
			}
	}
	virtual Space *copy(void) {
		return new SymmetricGccExample(*this);
	}
	void print(ostream& os) const {
		os << "\tSolution: \n";
		if (model == MODEL_SINGLE) {
			for (int i = 0; i < periods; i++) {
				for (int j = 0; j < weeks; j++) {
					cout << x[weeks*i + j] << " ";
				}
				cout << "\n";
			}
		} else if (model == MODEL_COUNT) {
			for (int i = 0; i < periods; i++) {
				for (int j = 0; j < weeks * 2; j += 2) {
					cout << y[weeks*2*i + j] << "v" << y[weeks*2*i + j + 1] << " ";
				}
				cout << "\n";
			}
		}
	}	
};

#endif