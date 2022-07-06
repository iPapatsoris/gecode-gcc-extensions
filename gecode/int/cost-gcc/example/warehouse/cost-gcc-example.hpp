#ifndef H_COST_GCC_EXAMPLE
#define H_COST_GCC_EXAMPLE

#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include "../BestBranch.hpp"
#include "read-input.hpp"

using namespace Gecode;
using namespace std;

class CountCostsExample : public Script {
protected:
	IntVarArray x;
	IntVar total;
	IntVar minCostFlowCost;

public:
	enum {
		MODEL_SINGLE, MODEL_MULTI, BRANCHING_SIMPLE, BRANCHING_CUSTOM
	};

	CountCostsExample(const InstanceOptions& opt);
	CountCostsExample(CountCostsExample &s) : Script(s) {
		x.update(*this, s.x);
		total.update(*this, s.total);
		minCostFlowCost.update(*this, s.minCostFlowCost);
	}
	virtual Space *copy(void) {
		return new CountCostsExample(*this);
	}	
	
	void print(ostream& os) const {
		os << "\tSolution: " << x << "\n" << cost() << endl;
	}

	virtual IntVar getMinCostFlowCost() const {
		return minCostFlowCost;
	}

	void constrain(const Space& _best);

	virtual IntVar cost(void) const {
		return total;
	}
	
};

#endif