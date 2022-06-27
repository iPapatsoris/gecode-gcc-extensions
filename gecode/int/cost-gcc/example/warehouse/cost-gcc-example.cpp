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
	int vars, fixed;
	IntSetArgs domain;
	IntArgs lowerBounds, upperBounds, vals, costs, demands;
	readInput(opt.instance(), vars, domain, vals, lowerBounds, upperBounds, costs, 
						fixed, demands);

	// Init domains
	x = IntVarArray(*this, vars);
	for (int i = 0; i < vars; i++) {
		x[i] = IntVar(*this, domain[i]);
	}

	// Local object handle for custom branching info from propagator
	BestBranch bestBranch(*this, x.size());

	total = IntVar(*this, 0, Int::Limits::max);
	// Boolean matrix of a row for each variable and a column for each value
	// Hold which variable is assigned to which value
	BoolVarArgs varValue(*this, costs.size(), 0, 1);
	Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
	for (int i = 0; i < m.height(); i++) {
		for (int j = 0; j < m.width(); j++) {
			rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
		}
	}

	// Hold open shops
	BoolVarArray open(*this, vals.size(), 0, 1);
	IntVarArray occurences(*this, vals.size());
	for (auto& o: occurences) {
		o = IntVar(*this, 0, vars);
	}
	for (int w=0; w < vals.size(); w++) {
		count(*this, x, w, IRT_EQ, occurences[w], opt.ipl());
		rel(*this, occurences[w], IRT_GQ, 1, eqv(open[w]), opt.ipl());
	}

	IntArgs fixedCosts;
	IntArgs fixedCostsArray;
	for (int f = 0; f <= fixed * vals.size(); f += fixed) {
		fixedCosts << f;
		if (f < fixed * vals.size()) {
			fixedCostsArray << fixed;
		}
	}
	minCostFlowCost = IntVar(*this, 0, Int::Limits::max);
	IntVar openCost(*this, IntSet(fixedCosts));
	linear(*this, costs, varValue, IRT_EQ, minCostFlowCost, opt.ipl());
	linear(*this, fixedCostsArray, open, IRT_EQ, openCost, opt.ipl());
	rel(*this, openCost + minCostFlowCost == total);

	Matrix<IntArgs> dd(demands, vals.size(), x.size());
	for (int w = 0; w < vals.size(); w++) {
		auto c = m.col(w);
		auto d = dd.col(w);
		linear(*this, d, c, IRT_LQ, upperBounds[w]);
	}

	auto simpleBranchVar = INT_VAR_SIZE_MIN();
	auto simpleBranchVal = INT_VAL_MIN();

	switch(opt.model()) {
		case MODEL_SINGLE:{
			countCosts(*this, x, vals, lowerBounds, upperBounds, costs, 
								minCostFlowCost, (opt.branching() == BRANCHING_CUSTOM ? 
																	&bestBranch : NULL), opt.ipl());
			if (opt.branching() == BRANCHING_CUSTOM) {
				branchBestVal(*this, x, bestBranch);
			} else {
				branch(*this, x, simpleBranchVar, simpleBranchVal);
			}
			break;
		}
		case MODEL_MULTI: {
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
	opt.branching(CountCostsExample::BRANCHING_SIMPLE,
	   	          "0", "use simple branching");
	opt.branching(CountCostsExample::BRANCHING_CUSTOM,
	    		      "1", "use custom branching");
	opt.branching(CountCostsExample::BRANCHING_CUSTOM);
	opt.ipl(IPL_DOM);
	opt.solutions(0);
	opt.parse(argc, argv);

	Script::run<CountCostsExample, BAB, InstanceOptions>(opt);

	return 0;
}
