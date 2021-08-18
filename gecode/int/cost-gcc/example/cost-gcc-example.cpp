#include "cost-gcc-example.hpp"
#include "../cost-gcc-post.hpp"
#include "LI.hpp"
#include "brancher.hpp"

CountCostsExample::CountCostsExample(const FileOptions& opt) : Script(opt) {
		int vars, cost;
		IntSetArgs domain;
		IntArgs lowerBounds, upperBounds, vals, costs;
		readInput(opt.file(), vars, domain, vals, lowerBounds, upperBounds, costs, 
							cost);
		x = IntVarArray(*this, vars);
		for (int i = 0; i < vars; i++) {
			x[i] = IntVar(*this, domain[i]);
		}

		// Local object handle, to branch using heuristic information provided by 
		// the propagator 
		LI li(*this, x.size());

		switch(opt.model()) {
			case MODEL_SINGLE:
				countCosts(*this, x, vals, lowerBounds, upperBounds, costs, cost, 
									 (opt.branch() ? &li : NULL),
									 opt.ipl());

				if (opt.branch()) {
					branchBestVal(*this, x, li);
				} else {
					branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
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

				branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
				break;
		}
	}

int main(int argc, char *argv[]) {
	FileOptions opt("Cost GCC");
	opt.model(CountCostsExample::MODEL_SINGLE,
						"single", "use single costgcc");
	opt.model(CountCostsExample::MODEL_MULTI,
						"multi", "use multiple constraints");
	opt.model(CountCostsExample::MODEL_SINGLE);
	opt.ipl(IPL_DOM);
	opt.solutions(0);
	opt.parse(argc, argv);

	Script::run<CountCostsExample, DFS, FileOptions>(opt);

	return 0;
}
