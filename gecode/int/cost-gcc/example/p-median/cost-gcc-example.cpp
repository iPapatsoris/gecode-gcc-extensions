#include "cost-gcc-example.hpp"
#include "../../cost-gcc-post.hpp"
#include "../BestBranch.hpp"
#include "../brancher.hpp"

void CountCostsExample::constrain(const Space& _best) {
	const CountCostsExample* best =
		dynamic_cast<const CountCostsExample*>(&_best);
	if (best == nullptr)
		throw DynamicCastFailed("Space::constrain");
	rel(*this, cost(), IRT_LE, best->cost().val());
}


CountCostsExample::CountCostsExample(const InstanceOptions& opt) : Script(opt) {
	int vars, totalOpen;
	IntSetArgs domain;
	IntArgs lowerBounds, upperBounds, vals, costs, demands;
	readInput(opt.instance(), vars, domain, vals, lowerBounds, upperBounds, costs, 
						totalOpen);
	
	x = IntVarArray(*this, vars);
	for (int i = 0; i < vars; i++) {
		x[i] = IntVar(*this, domain[i]);
	}

	BoolVarArgs varValue(*this, costs.size(), 0, 1);
	Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
	for (int i = 0; i < m.height(); i++) {
		for (int j = 0; j < m.width(); j++) {
			rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
		}
	}

	nvalues(*this, x, IRT_LQ, totalOpen, opt.ipl());

	minCostFlowCost = IntVar(*this, 0, Int::Limits::max);
	linear(*this, costs, varValue, IRT_EQ, minCostFlowCost, opt.ipl());

	auto simpleBranchVar = INT_VAR_SIZE_MIN();
	auto simpleBranchVal = INT_VAL_MIN();

	BestBranch bestBranch(*this, x.size());
	switch(opt.model()) {
		case MODEL_SINGLE:{
			countCosts(*this, x, vals, lowerBounds, upperBounds, costs, minCostFlowCost, 
									(opt.branching() == BRANCHING_CUSTOM ? &bestBranch : NULL),
									opt.ipl());
			if (opt.branching() == BRANCHING_CUSTOM) {
				branchBestVal(*this, x, bestBranch);
			} else {
				branch(*this, x, simpleBranchVar, simpleBranchVal);
			}
			break;
		}
		case MODEL_MULTI:{
			branch(*this, x, simpleBranchVar, simpleBranchVal);
			break;
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
