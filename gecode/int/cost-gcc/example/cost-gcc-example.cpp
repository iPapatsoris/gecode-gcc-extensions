#include "cost-gcc-example.hpp"
#include "../cost-gcc-post.hpp"
#include "BestBranch.hpp"
#include "brancher.hpp"

CountCostsExample::CountCostsExample(const InstanceOptions& opt) : Script(opt) {
	int vars, cost;
	IntSetArgs domain;
	IntArgs lowerBounds, upperBounds, vals, costs;
	readInput(opt.instance(), vars, domain, vals, lowerBounds, upperBounds, costs, 
						cost);
	x = IntVarArray(*this, vars);
	for (int i = 0; i < vars; i++) {
		x[i] = IntVar(*this, domain[i]);
	}
	IntVar costVar(*this, cost, cost);

	// Local object handle, to branch using heuristic information provided by 
	// the propagator 
	BestBranch bestBranch(*this, x.size());

	auto simpleBranchVar = INT_VAR_SIZE_MIN();
	auto simpleBranchVal = INT_VAL_MIN();

	switch(opt.model()) {
		case MODEL_SINGLE:
			countCosts(*this, x, vals, lowerBounds, upperBounds, costs, costVar, 
									(opt.branching() == BRANCHING_CUSTOM ? &bestBranch : NULL),
									opt.ipl());

			if (opt.branching() == BRANCHING_CUSTOM) {
				branchBestVal(*this, x, bestBranch);
			} else {
				branch(*this, x, simpleBranchVar, simpleBranchVal);
			}
			break;

		case MODEL_MULTI: 
			// Boolean flags for each variable-value assignment. Needed for cost
			BoolVarArgs varValue(*this, costs.size(), 0, 1);
			Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
			for (int i = 0; i < m.height(); i++) {
				for (int j = 0; j < m.width(); j++) {
					rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
				}
			}
			// Constrain the cost
			linear(*this, costs, varValue, IRT_LQ, cost, opt.ipl());

			// Constrain the value bounds
			IntSetArgs bounds;
			for (int i = 0; i < vals.size(); i++) {
				bounds << IntSet(lowerBounds[i], upperBounds[i]);
			}
			count(*this, x, bounds, vals, opt.ipl());

			branch(*this, x, simpleBranchVar, simpleBranchVal);
			break;
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

	Script::run<CountCostsExample, DFS, InstanceOptions>(opt);

	return 0;
}
