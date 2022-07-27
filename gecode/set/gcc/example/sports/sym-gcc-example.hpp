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
		for (int i = 0; i < periods; i++) {
			for (int j = 0; j < weeks; j++) {
				if (model == MODEL_COUNT) {
					cout << y[weeks*i + j] << " ";
				} else if (model == MODEL_SINGLE) {
					cout << x[weeks*i + j] << " ";
				}
			}
			cout << "\n";
		}
	}
};

#endif