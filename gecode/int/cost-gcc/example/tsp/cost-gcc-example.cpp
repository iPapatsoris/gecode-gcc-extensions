#include "cost-gcc-example.hpp"
#include "../../cost-gcc-post.hpp"
#include "../LI.hpp"
#include "../brancher.hpp"
#include <gecode/float.hh>

 void
  CountCostsExample::constrain(const Space& _best) {
    const CountCostsExample* best =
      dynamic_cast<const CountCostsExample*>(&_best);
    if (best == nullptr)
      throw DynamicCastFailed("IntMinimizeSpace::constrain");
    rel(*this, cost(), IRT_LE, best->cost().val());
  }


CountCostsExample::CountCostsExample(const FileOptions& opt,
int vars, IntSetArgs domain, IntArgs lowerBounds, IntArgs upperBounds, IntArgs vals, IntArgs costs, vector<int>& costsD) : IntMinimizeScript(opt) {
		x = IntVarArray(*this, vars);
		for (int i = 0; i < vars; i++) {
			x[i] = IntVar(*this, domain[i]);
			// cout << domain[i] << endl;
		}

		// Local object handle, to branch using heuristic information provided by 
		// the propagator 
		LI li(*this, x.size());
		total = IntVar(*this, 0, 100000000);
		switch(opt.model()) {
			case MODEL_SINGLE:{
				BoolVarArgs varValue(*this, costs.size(), 0, 1);
				Matrix<BoolVarArgs> m(varValue, vals.size(), x.size());
				for (int i = 0; i < m.height(); i++) {
					for (int j = 0; j < m.width(); j++) {
						rel(*this, x[i], IRT_EQ, vals[j], m(j, i), opt.ipl());
					}
				}
				
				// Constrain the cost
				linear(*this, costs, varValue, IRT_EQ, total);

				circuit(*this, x);
				countCosts(*this, x, vals, lowerBounds, upperBounds, costsD, total, 
									 (opt.branch() ? &li : NULL),
									 opt.ipl());

				if (opt.branch()) {
					branchBestVal(*this, x, li);
				} else {
					branch(*this, x, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
				}
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

	int vars;
		IntSetArgs domain;
		IntArgs lowerBounds, upperBounds, vals, costs;
		vector<int> costsD;
	readInput(opt.file(), vars, domain, vals, lowerBounds, upperBounds, costs, costsD);

	// cout << vars;
	// cout << lowerBounds;
	// cout << upperBounds;
	// cout << vals;
	//  cout << costs;
	// cout << domain;
	//  exit(1);
						
				
	//Script::run<CountCostsExample, DFS, FileOptions>(opt);
	// TODO: all this is not optimized, can be optimized to keep input data
	CountCostsExample* m = new CountCostsExample(opt, vars, domain, lowerBounds, upperBounds, vals, costs, costsD);
	BAB<CountCostsExample> e(m);
	CountCostsExample *s = e.next();

	const int inf = 10000000;
	int bound = inf;

	while (s != NULL) {
		s->print(cout);
		bound = s->cost().val();
		cout << bound << endl;
		delete s;
		s = e.next();
		continue;
	}
	if (bound == inf) {
		cout << "no solution" << endl;
	} else {
		cout << "optimal solution " << bound << endl;
	}
	delete m;

	return 0;
}
