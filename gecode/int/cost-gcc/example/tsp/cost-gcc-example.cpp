#include "cost-gcc-example.hpp"
#include "../../cost-gcc-post.hpp"
#include "../LI.hpp"
#include "../brancher.hpp"

 void CountCostsExample::constrain(const Space& _best) {
	const CountCostsExample* best =
		dynamic_cast<const CountCostsExample*>(&_best);
	if (best == nullptr)
		throw DynamicCastFailed("Space::constrain");
	rel(*this, cost(), IRT_LE, best->cost().val());
}

CountCostsExample::CountCostsExample(const InstanceOptions& opt) : Script(opt) {
	int vars;
	IntSetArgs domain;
	IntArgs lowerBounds, upperBounds, vals, costs;
	readInput(opt.instance(), vars, domain, vals, lowerBounds, upperBounds, 
						costs);
	
	// Init domain
	x = IntVarArray(*this, vars);
	for (int i = 0; i < vars; i++) {
		x[i] = IntVar(*this, domain[i]);
	}

	// Local object handle, to branch using heuristic information provided by 
	// the propagator 
	LI li(*this, x.size());
	
	// Boolean matrix of a row for each variable and a column for each value
	// Hold which variable is assigned to which value
	BoolVarArgs varValue(*this, costs.size(), 0, 1);
	Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
	for (int i = 0; i < m.height(); i++) {
		for (int j = 0; j < m.width(); j++) {
			rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
		}
	}
	
	// Hold the cost
	total = IntVar(*this, 0, Int::Limits::max);
	linear(*this, costs, varValue, IRT_EQ, total);
	
	auto simpleBranchVar = INT_VAR_SIZE_MIN();
	auto simpleBranchVal = INT_VAL_MIN();

	switch(opt.model()) {
		case MODEL_SINGLE:{
			circuit(*this, x);
			countCosts(*this, x, vals, lowerBounds, upperBounds, costs, total, 
									(opt.branching() == BRANCHING_CUSTOM ? &li : NULL),
									opt.ipl());

			if (opt.branching() == BRANCHING_CUSTOM) {
				branchBestVal(*this, x, li);
			} else if (opt.branching() == BRANCHING_SIMPLE) {
				branch(*this, x, simpleBranchVar, simpleBranchVal);
			}
			break;
		}
		case MODEL_MULTI: {
			circuit(*this, x);
			
			// Constrain the value bounds
			IntSetArgs bounds;
			for (int i = 0; i < vals.size(); i++) {
				bounds << IntSet(lowerBounds[i], upperBounds[i]);
			}
			count(*this, x, bounds, vals, opt.ipl());
			branch(*this, x, simpleBranchVar, simpleBranchVal);
		}
	}
}

int main(int argc, char *argv[]) {
	InstanceOptions opt("Cost GCC");
	opt.model(CountCostsExample::MODEL_SINGLE,
						"single", "use single costgcc");
	opt.model(CountCostsExample::MODEL_MULTI,
						"multi", "use multiple constraints");
	opt.model(CountCostsExample::MODEL_SINGLE);
	opt.ipl(IPL_DOM);
	opt.solutions(0);
	opt.branching(CountCostsExample::BRANCHING_SIMPLE,
	          "0", "use simple branching");
	opt.branching(CountCostsExample::BRANCHING_CUSTOM,
	          "1", "use custom branching");
	opt.branching(CountCostsExample::BRANCHING_CUSTOM);
	opt.parse(argc, argv);

	Script::run<CountCostsExample, BAB, InstanceOptions>(opt);
	return 0;
}
