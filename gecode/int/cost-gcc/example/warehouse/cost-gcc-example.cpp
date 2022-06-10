#include "cost-gcc-example.hpp"
#include "../../cost-gcc-post.hpp"
#include "../LI.hpp"
#include "../brancher.hpp"

CountCostsExample::CountCostsExample(const FileOptions& opt, int cost, int previousBest,
int vars, IntSetArgs domain, IntArgs lowerBounds, IntArgs upperBounds, IntArgs vals, IntArgs costs, int fixed) : IntMinimizeScript(opt) {
		x = IntVarArray(*this, vars);
		for (int i = 0; i < vars; i++) {
			x[i] = IntVar(*this, domain[i]);
			// cout << domain[i] << endl;
		}

		// cout << costs;
		// cout << lowerBounds;
		// cout << upperBounds;
		// cout << vals;
		// the propagator 
		LI li(*this, x.size());
		total = IntVar(*this, 0, previousBest-1);
		switch(opt.model()) {
			case MODEL_SINGLE:{
				BoolVarArgs varValue(*this, costs.size(), 0, 1);
				Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
				for (int i = 0; i < m.height(); i++) {
					for (int j = 0; j < m.width(); j++) {
						rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
					}
				}

				open = BoolVarArray(*this, vals.size(), 0, 1);
				for (int s=0; s<vars; s++)
					element(* this, open, x[s], 1, opt.ipl());

				IntArgs fixedCosts;
				IntArgs fixedCostsArray;
				for (unsigned int f = 0; f <= fixed * vals.size(); f += fixed) {
					fixedCosts << f;
					if (f < fixed * vals.size()) {
						fixedCostsArray << fixed;
					}
				}
				minCostFlowCost = IntVar(*this, 0, cost);
				openCost = IntVar(*this, IntSet(fixedCosts));
				linear(*this, costs, varValue, IRT_EQ, minCostFlowCost, opt.ipl());
				linear(*this, fixedCostsArray, open, IRT_EQ, openCost, opt.ipl());
				rel(*this, openCost + minCostFlowCost == total, opt.ipl());
				// rel(*this, total < previousBest, opt.ipl());
				countCosts(*this, x, vals, lowerBounds, upperBounds, costs, cost, 
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

	int vars, cost;
		IntSetArgs domain;
		IntArgs lowerBounds, upperBounds, vals, costs;
//	readInput(opt.file(), vars, domain, vals, lowerBounds, upperBounds, costs, 
			//				cost);

	const int n_warehouses = 5;
/// Number of stores
const int n_stores = 10;

/// Fixed cost for one warehouse
const int c_fixed = 30;
const int maxFixed = n_warehouses * c_fixed;

/// Capacity of a single warehouse
const int capacity[n_warehouses] = {
  1, 4, 2, 1, 3
};

/// Cost for supply a store by a warehouse
const int c_supply[n_stores][n_warehouses] = {
  {20, 24, 11, 25, 30},
  {28, 27, 82, 83, 74},
  {74, 97, 71, 96, 70},
  { 2, 55, 73, 69, 61},
  {46, 96, 59, 83,  4},
  {42, 22, 29, 67, 59},
  { 1,  5, 73, 59, 56},
  {10, 73, 13, 43, 96},
  {93, 35, 63, 85, 46},
  {47, 65, 55, 71, 95}
};
	// return 1;

	vars = n_stores;
	for (unsigned int i = 0; i < vars; i++) {
		domain << IntSet(0, n_warehouses-1);
	}

	for (unsigned int i = 0; i < n_warehouses; i++) {
		lowerBounds << 0;
		upperBounds << capacity[i];
		vals << i;
	}

	for (unsigned int i = 0; i < n_stores; i++) {
		for (unsigned int j = 0; j < n_warehouses; j++) {
			costs << c_supply[i][j];
		}
	}
	const int inf = 10000000;

	cost = inf;
				
	//Script::run<CountCostsExample, DFS, FileOptions>(opt);
	// TODO: all this is not optimized, can be optimized to keep input data
	CountCostsExample* m = new CountCostsExample(opt, cost, cost, vars, domain, lowerBounds, upperBounds, vals, costs, c_fixed);
	DFS<CountCostsExample> e(m);
	CountCostsExample *s = e.next();

	int boundBAB = inf;
	int boundMinCostFlow = inf;
	while (s != NULL) {
		s->printOpen();
		boundBAB = s->getTotal();
		boundMinCostFlow = s->getMinCostFlowCost() + maxFixed;
		s->print(cout);
		cout << boundBAB << endl;
		s = e.next();
		//continue;
		delete s;
		delete m;
		// return 1;
		m = new CountCostsExample(opt, boundMinCostFlow, boundBAB, vars, domain, lowerBounds, upperBounds, vals, costs, c_fixed);
		DFS<CountCostsExample> ee(m);
		s = ee.next();
	}
	if (boundBAB == inf) {
		cout << "no solution" << endl;
	} else {
		cout << "optimal solution " << boundBAB << endl;
	}
	delete m;

	return 0;
}
