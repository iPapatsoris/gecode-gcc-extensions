#include "sym-gcc-example.hpp"
#include "../sym-gcc-post.hpp"

#include "LI.hpp"
#include "brancher.hpp"

SymmetricGccExample::SymmetricGccExample(const InstanceOptions& opt) 
	: Script(opt), model((Model) opt.model()) {
	IntSetArgs domain;
	IntArgs vals, lowerValBounds, upperValBounds, lowerVarBounds, upperVarBounds;
	readInput(opt.instance(), varsCount, domain, vals, lowerValBounds, upperValBounds, 
						lowerVarBounds, upperVarBounds);
	valsCount = vals.size();
	
	x = SetVarArray(*this, 0);
	y = IntVarArray(*this, 0);
	z = BoolVarArray(*this, 0);

	/*
	for (auto& x: x) {
		cout << "unkown values\n";
		for (SetVarUnknownValues i(x); i(); ++i)
			std::cout << i.val() << " ";
		cout << "unkown ranges\n";
		for (SetVarUnknownRanges i(x); i(); ++i)
			std::cout << i.min() << ".." << i.max() << " ";
		cout << "lowest upper bound values\n";
		for (SetVarLubValues i(x); i(); ++i)
			std::cout << i.val() << " ";
		cout << "lowest upper bound ranges\n";
		for (SetVarLubRanges i(x); i(); ++i)
			std::cout << i.min() << ".." << i.max() << " ";
		cout << "\n\n";
	}

	dom(*this, x[4], SRT_DISJ, IntSet(4, 4));
	for (SetVarLubValues i(x[4]); i(); ++i)
			std::cout << i.val() << " ";
		cout << "\n\n";*/

	// Local object handle, to branch using heuristic information provided by 
	// the propagator 
	LI li(*this, varsCount);

	auto simpleBranchVar = SET_VAR_DEGREE_MIN();
	auto simpleBranchVal = SET_VAL_MIN_EXC();

	switch(opt.model()) {
		case MODEL_SINGLE:
			x = SetVarArray(*this, varsCount);
			for (int i = 0; i < varsCount; i++) {
				x[i] = SetVar(*this, IntSet::empty, domain[i], lowerVarBounds[i], upperVarBounds[i]);
			}
			symmetricGCC(*this, x, vals, lowerValBounds, upperValBounds, 
										lowerVarBounds, upperVarBounds, opt.branching() ? &li : NULL, 
										opt.ipl());
		//dom(*this, x[4], SRT_EQ, IntSet{3, 5});
			if (opt.branching()) {
				bestval(*this, x, li);
			} else { 
				branch(*this, x, simpleBranchVar, simpleBranchVal);
			}
			break;

		case MODEL_COUNT: {
			// Boolean flags for each variable-value assignment. Needed for cost
			y = IntVarArray(*this, varsCount * vals.size(), 0, 1);
			Matrix<IntVarArray> m(y, vals.size(), varsCount);
			// Assign 0 to the illegal domain values
			for (int i = 0; i < m.height(); i++) {
				for (int j = 0; j < m.width(); j++) {
					if (!domain[i].in(vals[j])) {
						dom(*this, m(j,i), 0, 0);
					}
				}
			}
			for (int i = 0; i < m.height(); i++) {
				count(*this, m.row(i), 1, IRT_GQ, lowerVarBounds[i], opt.ipl());
				count(*this, m.row(i), 1, IRT_LQ, upperVarBounds[i], opt.ipl());
				//count(*this, m.row(i), varBounds[i], 1);
			}
			for (int j = 0; j < m.width(); j++) {
				count(*this, m.col(j), 1, IRT_GQ, lowerValBounds[j], opt.ipl());
				count(*this, m.col(j), 1, IRT_LQ, upperValBounds[j], opt.ipl());
				//count(*this, m.col(j), bounds, ones);
			}
			branch(*this, y, INT_VAL_MAX());
			break;
		}
		case MODEL_LINEAR:
			z = BoolVarArray(*this, varsCount * vals.size(), 0, 1);
			Matrix<BoolVarArray> mB(z, vals.size(), varsCount);
			// Assign 0 to the illegal domain values
			for (int i = 0; i < mB.height(); i++) {
				for (int j = 0; j < mB.width(); j++) {
					if (!domain[i].in(vals[j])) {
						rel(*this, mB(j,i), IRT_EQ, 0);
					}
				}
			}
			for (int i = 0; i < mB.height(); i++) {
				linear(*this, mB.row(i), IRT_GQ, lowerVarBounds[i], opt.ipl());
				linear(*this, mB.row(i), IRT_LQ, upperVarBounds[i], opt.ipl());
				//count(*this, m.row(i), varBounds[i], 1);
			}
			for (int j = 0; j < mB.width(); j++) {
				linear(*this, mB.col(j), IRT_GQ, lowerValBounds[j], opt.ipl());
				linear(*this, mB.col(j), IRT_LQ, upperValBounds[j], opt.ipl());
				//count(*this, m.col(j), bounds, ones);
			}
			branch(*this, z, BOOL_VAL_MAX());
			break;
	}
}
	

int main(int argc, char *argv[]) {
	InstanceOptions opt("Cost GCC");
	opt.model(SymmetricGccExample::MODEL_SINGLE,
						"single", "use single costgcc");
	opt.model(SymmetricGccExample::MODEL_COUNT,
						"count", "use count constraints");
	opt.model(SymmetricGccExample::MODEL_LINEAR,
						"linear", "use linear constraints");
	opt.model(SymmetricGccExample::MODEL_SINGLE);
	opt.branching(SymmetricGccExample::BRANCHING_SIMPLE,
	          "0", "use simple branching");
	opt.branching(SymmetricGccExample::BRANCHING_CUSTOM,
	          "1", "use custom branching");
	opt.branching(SymmetricGccExample::BRANCHING_CUSTOM);
	opt.ipl(IPL_DOM);
	opt.solutions(0);
	opt.parse(argc, argv);

	Script::run<SymmetricGccExample, DFS, InstanceOptions>(opt);

	return 0;
}
