#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include "cost-gcc-post.hpp"

using namespace Gecode;

class CountCostsExample : public Space {
protected:
	IntVarArray x;

public:
	CountCostsExample(void) : x(*this, 3) {
		x[0] = IntVar(*this, IntSet({6, 9, 777, 888}));
		x[1] = IntVar(*this, IntSet({10, 7, 6, 777, 9}));
		x[2] = IntVar(*this, IntSet({777, 9}));

		IntArgs lowerBounds = {0, 1, 0, 1};
		IntArgs upperBounds = {2, 1, 1, 4};

		IntArgs vals = {9, 7, 10, 777};

		IntArgs costs = {
				14, 15, 16, 17,
				11, 12, 13, 14,
				11, 12, 13, 14
		};

		countCosts(*this, x, vals, lowerBounds, upperBounds, costs, 42, IntPropLevel::IPL_DOM);
		//branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
	}

	CountCostsExample(CountCostsExample &s) : Space(s) {
		x.update(*this, s.x);
	}
	virtual Space *copy(void) {
		return new CountCostsExample(*this);
	}
	void print(void) const {
		std::cout << "Solution: " << x << std::endl;
	}
};

int main(int argc, char *argv[]) {
	CountCostsExample *m = new CountCostsExample;
	return 0;

	DFS<CountCostsExample> e(m);
	delete m;

	// search and print all solutions
	while (CountCostsExample *s = e.next()) {
		s->print();
		delete s;
	}

	return 0;
}
