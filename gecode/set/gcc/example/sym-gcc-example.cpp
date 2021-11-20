#include "sym-gcc-example.hpp"
#include "../sym-gcc-post.hpp"
//#include "LI.hpp"
//#include "brancher.hpp"

SymmetricGccExample::SymmetricGccExample(const FileOptions& opt) : Script(opt) {
		int vars;
		IntSetArgs domain;
		IntArgs lowerValBounds, upperValBounds, vals, lowerVarBounds, upperVarBounds;
		readInput(opt.file(), vars, domain, vals, lowerValBounds, upperValBounds, 
						  lowerVarBounds, upperVarBounds);
		
		x = SetVarArray(*this, 0);
		y = IntVarArray(*this, 0);

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
		// LI li(*this, x.size());

		switch(opt.model()) {
			case MODEL_SINGLE:
				x = SetVarArray(*this, vars);
				for (int i = 0; i < vars; i++) {
					x[i] = SetVar(*this, IntSet::empty, domain[i], lowerVarBounds[i], upperVarBounds[i]);
				}
				symmetricGCC(*this, x, vals, lowerValBounds, upperValBounds, 
										 lowerVarBounds, upperVarBounds,
									 opt.ipl());
			//dom(*this, x[4], SRT_EQ, IntSet{3, 5});
				/*if (opt.branch()) {
					branchBestVal(*this, x, li);
				} else {*/ 
					branch(*this, x, SET_VAR_DEGREE_MIN(), SET_VAL_MIN_INC());
				//}
				break;

			case MODEL_MULTI: 
				// Constrain the value bounds
				IntSetArgs bounds;
				for (int i = 0; i < vals.size(); i++) {
					bounds << IntSet(lowerValBounds[i], upperValBounds[i]);
				}

				// Constrain the var bounds
				IntSetArgs varBounds;
				IntArgs ones;
				for (int i = 0; i < vars; i++) {
					varBounds << IntSet(lowerVarBounds[i], upperVarBounds[i]);
					ones << 1;
				}

				// Boolean flags for each variable-value assignment. Needed for cost
				y = IntVarArray(*this, vars * vals.size(), 0, 1);
				Matrix<IntVarArray> m(y, vals.size(), vars);
				// Assign 0 to the illegal domain values
				for (int i = 0; i < m.height(); i++) {
					for (int j = 0; j < m.width(); j++) {
						if (!domain[i].in(vals[j])) {
							dom(*this, m(j,i), 0, 0);
						}
					}
				}
				for (int i = 0; i < m.height(); i++) {
					count(*this, m.row(i), 1, IRT_GQ, lowerVarBounds[i]);
					count(*this, m.row(i), 1, IRT_LQ, upperVarBounds[i]);
					//count(*this, m.row(i), varBounds[i], 1);
				}
				for (int j = 0; j < m.width(); j++) {
					count(*this, m.col(j), 1, IRT_GQ, lowerValBounds[j]);
					count(*this, m.col(j), 1, IRT_LQ, upperValBounds[j]);
					//count(*this, m.col(j), bounds, ones);
				}
				branch(*this, y, INT_VAR_SIZE_MIN(), INT_VAL_MIN());
				break;
		}
	}

int main(int argc, char *argv[]) {
	FileOptions opt("Cost GCC");
	opt.model(SymmetricGccExample::MODEL_SINGLE,
						"single", "use single costgcc");
	opt.model(SymmetricGccExample::MODEL_MULTI,
						"multi", "use multiple constraints");
	opt.model(SymmetricGccExample::MODEL_SINGLE);
	opt.ipl(IPL_DOM);
	opt.solutions(0);
	opt.parse(argc, argv);

	Script::run<SymmetricGccExample, DFS, FileOptions>(opt);

	return 0;
}
