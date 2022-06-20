#include "cost-gcc-example.hpp"
#include "../../cost-gcc-post.hpp"
#include "../LI.hpp"
#include "../brancher.hpp"

 void
  CountCostsExample::constrain(const Space& _best) {
    const CountCostsExample* best =
      dynamic_cast<const CountCostsExample*>(&_best);
    if (best == nullptr)
      throw DynamicCastFailed("IntMinimizeSpace::constrain");
    rel(*this, cost(), IRT_LE, best->cost().val());
		// rel(*this, getMinCostFlowCost() <= best->getMinCostFlowCost() + getOpenCost());
		// cout << getOpenCost() << endl;
  }


CountCostsExample::CountCostsExample(const FileOptions& opt, int cost, int previousBest,
int vars, IntSetArgs domain, IntArgs lowerBounds, IntArgs upperBounds, IntArgs vals, IntArgs costs, int totalOpen, IntArgs demands) : IntMinimizeScript(opt) {
		x = IntVarArray(*this, vars);
		for (int i = 0; i < vars; i++) {
			x[i] = IntVar(*this, domain[i]);
			// cout << domain[i] << endl;
		}

		// cout << vars;
		// cout << costs;
		// cout << lowerBounds;
		// cout << upperBounds;
		// cout << vals;
		// exit(1);
		// the propagator 
		LI li(*this, x.size());
		switch(opt.model()) {
			case MODEL_SINGLE:{
				BoolVarArgs varValue(*this, costs.size(), 0, 1);
				Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
				for (int i = 0; i < m.height(); i++) {
					for (int j = 0; j < m.width(); j++) {
						rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
					}
				}

				nvalues(*this, x, IRT_LQ, totalOpen, opt.ipl());

				minCostFlowCost = IntVar(*this, 0, cost);
				linear(*this, costs, varValue, IRT_EQ, minCostFlowCost, opt.ipl());
				countCosts(*this, x, vals, lowerBounds, upperBounds, costs, minCostFlowCost, 
									 (opt.branch() ? &li : NULL),
									 opt.ipl());
				if (opt.branch()) {
					branchBestVal(*this, x, li);
				} else {
					branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
				}

				// cout << openCost << minCostFlow << total;
				break;
			}

			case MODEL_MULTI:{
				BoolVarArgs varValue(*this, costs.size(), 0, 1);
				Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
				for (int i = 0; i < m.height(); i++) {
					for (int j = 0; j < m.width(); j++) {
						rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
					}
				}

				nvalues(*this, x, IRT_LQ, totalOpen, opt.ipl());

				minCostFlowCost = IntVar(*this, 0, cost);
				linear(*this, costs, varValue, IRT_EQ, minCostFlowCost, opt.ipl());
				branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());

				// cout << openCost << minCostFlow << total;
				break;
			}
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

	int vars, fixed;
		IntSetArgs domain;
		IntArgs lowerBounds, upperBounds, vals, costs, demands;
	readInput(opt.file(), vars, domain, vals, lowerBounds, upperBounds, costs, 
							fixed);

	// cout << vars;
	// cout << domain;
	// cout << vals; 
	// cout << lowerBounds;
	// cout << upperBounds;
	// cout << costs;

	// int ij = 0;
	// for (unsigned int i = 0; i < vars; i++) {
	// 	cout << "{";
	// 	for (unsigned int j = 0; j < vals.size(); j++) {
	// 		cout << (costs[ij] == 0 ? 100000000 : costs[ij]) << ", ";
	// 		ij++;
	// 	}
	// 	cout << "}, \n";
	// }

	//  return 0;

// 	const int n_warehouses = 5;
// /// Number of stores
// const int n_stores = 10;

// /// Fixed cost for one warehouse
// const int c_fixed = 30;
// const int maxFixed = n_warehouses * c_fixed;

// /// Capacity of a single warehouse
// const int capacity[n_warehouses] = {
//   1, 4, 2, 1, 3
// };

// /// Cost for supply a store by a warehouse
// const int c_supply[n_stores][n_warehouses] = {
//   {20, 24, 11, 25, 30},
//   {28, 27, 82, 83, 74},
//   {74, 97, 71, 96, 70},
//   { 2, 55, 73, 69, 61},
//   {46, 96, 59, 83,  4},
//   {42, 22, 29, 67, 59},
//   { 1,  5, 73, 59, 56},
//   {10, 73, 13, 43, 96},
//   {93, 35, 63, 85, 46},
//   {47, 65, 55, 71, 95}
// };
	// return 1;

	// vars = n_stores;
	// for (unsigned int i = 0; i < vars; i++) {
	// 	domain << IntSet(0, n_warehouses-1);
	// }

	// for (unsigned int i = 0; i < n_warehouses; i++) {
	// 	lowerBounds << 0;
	// 	upperBounds << capacity[i];
	// 	vals << i;
	// }

	// for (unsigned int i = 0; i < n_stores; i++) {
	// 	for (unsigned int j = 0; j < n_warehouses; j++) {
	// 		costs << c_supply[i][j];
	// 	}
	// }
	const int inf = 10000000;
  int boundMinCostFlow = inf;
	int boundBAB = inf;

	//Script::run<CountCostsExample, DFS, FileOptions>(opt);
	// TODO: all this is not optimized, can be optimized to keep input data
	CountCostsExample* m = new CountCostsExample(opt, boundMinCostFlow, boundBAB, vars, domain, lowerBounds, upperBounds, vals, costs, fixed, demands);
	BAB<CountCostsExample> e(m);
	CountCostsExample *s = e.next();

	while (s != NULL) {
		boundBAB = s->cost().val();
		s->print(cout);
		cout << boundBAB << endl;
		s = e.next();
		continue;
	}
	if (boundBAB == inf) {
		cout << "no solution" << endl;
	} else {
		cout << "optimal solution " << boundBAB << endl;
	}
	delete m;

	return 0;
}
