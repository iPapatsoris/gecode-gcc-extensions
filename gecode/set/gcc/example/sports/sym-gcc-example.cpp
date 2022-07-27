#include "sym-gcc-example.hpp"
#include "../../sym-gcc-post.hpp"

using namespace std;

SymmetricGccExample::SymmetricGccExample(const InstanceOptions& opt) 
	: Script(opt), model((Model) opt.model()) {
	IntSetArgs domain;
	IntArgs vals, lowerValBoundsWeek, upperValBoundsWeek, 
					lowerVarBoundsWeek, upperVarBoundsWeek, lowerValBoundsPeriod,
					upperValBoundsPeriod, lowerVarBoundsPeriod, upperVarBoundsPeriod;
	int teams = atoi(opt.instance());
	periods = teams / 2;
	weeks = teams - 1;
	int varsCount = weeks * periods;
	int valsCount = teams;

	for (int v = 0; v < valsCount; v++) {
		vals << v;
		lowerValBoundsWeek << 1;
		upperValBoundsWeek << 1;
		lowerValBoundsPeriod << 0;
		upperValBoundsPeriod << 2;
	}

	for (int w = 0; w < weeks; w++) {
		lowerVarBoundsPeriod << 2;
		upperVarBoundsPeriod << 2;
	}

	for (int p = 0; p < periods; p++) {
		lowerVarBoundsWeek << 2;
		upperVarBoundsWeek << 2;
	}

	x = SetVarArray(*this, 0);
	y = IntVarArray(*this, 0);

	auto simpleBranchVar = SET_VAR_DEGREE_MIN();
	auto simpleBranchVal = SET_VAL_MIN_EXC();

	switch(opt.model()) {
		case MODEL_SINGLE: {
			x = SetVarArray(*this, varsCount);
			for (int i = 0; i < varsCount; i++) {
				x[i] = SetVar(*this, IntSet::empty, IntSet(0, teams - 1), 2, 2);
			}
			Matrix<SetVarArray> m(x, weeks, periods);
			for (int w = 0; w < weeks; w++) {
				symmetricGCC(*this, m.col(w), vals, lowerValBoundsWeek, upperValBoundsWeek,
								 lowerVarBoundsWeek, upperVarBoundsWeek, NULL, opt.ipl()); 
			}
			for (int p = 0; p < periods; p++) {
				symmetricGCC(*this, m.row(p), vals, lowerValBoundsPeriod, 
								 upperValBoundsPeriod, lowerVarBoundsPeriod, 
								 upperVarBoundsPeriod, NULL, opt.ipl()); 
			}
			for (int i = 0; i < varsCount - 1; i++) {
				for (int j = i + 1; j < varsCount; j++)
				rel(*this, x[i], SRT_NQ, x[j]);
			}
			branch(*this, x, simpleBranchVar, simpleBranchVal);
			break;
		}
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
				// count(*this, m.row(i), 1, IRT_GQ, lowerVarBounds[i], opt.ipl());
				// count(*this, m.row(i), 1, IRT_LQ, upperVarBounds[i], opt.ipl());
				//count(*this, m.row(i), varBounds[i], 1);
			}
			for (int j = 0; j < m.width(); j++) {
				// count(*this, m.col(j), 1, IRT_GQ, lowerValBounds[j], opt.ipl());
				// count(*this, m.col(j), 1, IRT_LQ, upperValBounds[j], opt.ipl());
				//count(*this, m.col(j), bounds, ones);
			}
			branch(*this, y, INT_VAL_MAX());
			break;
		}
	}
}
	

int main(int argc, char *argv[]) {
	InstanceOptions opt("Cost GCC");
	opt.model(SymmetricGccExample::MODEL_SINGLE,
						"single", "use single symgcc");
	opt.model(SymmetricGccExample::MODEL_COUNT,
						"count", "use count constraints");
	opt.model(SymmetricGccExample::MODEL_SINGLE);
	opt.ipl(IPL_DOM);
	opt.solutions(0);
	opt.parse(argc, argv);

	Script::run<SymmetricGccExample, DFS, InstanceOptions>(opt);

	return 0;
}
