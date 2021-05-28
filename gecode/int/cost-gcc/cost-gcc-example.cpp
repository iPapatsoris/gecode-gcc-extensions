#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include "cost-gcc-post.hpp"

using namespace Gecode;
using namespace std;

const int vars = 3;
const int cost = 41;

IntSetArgs domain = {
	IntSet({6, 9, 777, 888}),
	IntSet({10, 7, 6, 777, 9}),
	IntSet({777, 9})
};

IntArgs lowerBounds = {0, 0, 0, 1};
IntArgs upperBounds = {2, 1, 1, 4};

IntArgs vals = {9, 7, 10, 777};

IntArgs costs = {
		14, 15, 16, 17,
		11, 12, 13, 14,
		11, 12, 13, 14};

class CountCostsExample : public Script {
protected:
	IntVarArray x;

public:
	enum {
		MODEL_SINGLE, MODEL_MULTI
	};

	CountCostsExample(const Options& opt) : Script(opt) {
		x = IntVarArray(*this, vars);
		for (int i = 0; i < vars; i++) {
			x[i] = IntVar(*this, domain[i]);
		}

		switch(opt.model()) {
			case MODEL_SINGLE:
				countCosts(*this, x, vals, lowerBounds, upperBounds, costs, cost, 
									 opt.ipl());

				branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
				break;

			case MODEL_MULTI: 
				// Boolean flags for each variable-value assignment. Needed for cost
				BoolVarArgs varValue(*this, costs.size(), 0, 1);
				Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
				for (int i = 0; i < m.height(); i++) {
					for (int j = 0; j < m.width(); j++) {
						rel(*this, x[i], IRT_EQ, vals[j], m(j, i));
					}
				}
				// Constrain the cost
				linear(*this, costs, varValue, IRT_LQ, cost);

				// Constrain the value bounds
				IntSetArgs bounds;
				for (int i = 0; i < vals.size(); i++) {
					bounds << IntSet(lowerBounds[i], upperBounds[i]);
				}
				count(*this, x, bounds, vals, opt.ipl());

				branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
				break;
		}
	}

	CountCostsExample(CountCostsExample &s) : Script(s) {
		x.update(*this, s.x);
	}
	virtual Space *copy(void) {
		return new CountCostsExample(*this);
	}
	void print(ostream& os) const {
		os << "\tSolution: " << x << endl;
	}
};

int main(int argc, char *argv[]) {
	Options opt("Cost GCC");
	opt.model(CountCostsExample::MODEL_SINGLE,
						"single", "use single costgcc");
	opt.model(CountCostsExample::MODEL_MULTI,
						"multi", "use multiple constraints");
	opt.model(CountCostsExample::MODEL_SINGLE);
	opt.ipl(IPL_DOM);
	opt.solutions(0);
	opt.parse(argc, argv);

	Script::run<CountCostsExample, DFS, Options>(opt);

	return 0;
}
