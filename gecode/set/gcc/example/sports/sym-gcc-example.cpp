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

	x = SetVarArray(*this, 0);
	y = IntVarArray(*this, 0);

	auto simpleBranchVar = SET_VAR_DEGREE_MIN();
	auto simpleBranchVal = SET_VAL_MIN_EXC();

	int valsCount = teams;
	for (int v = 0; v < valsCount; v++) {
		vals << v;
		lowerValBoundsWeek << 1;
		upperValBoundsWeek << 1;
		lowerValBoundsPeriod << 0;
		upperValBoundsPeriod << 2;
	}

	switch(opt.model()) {
		case MODEL_SINGLE: {
			int varsCount = weeks * periods;


			for (int w = 0; w < weeks; w++) {
				lowerVarBoundsPeriod << 2;
				upperVarBoundsPeriod << 2;
			}

			for (int p = 0; p < periods; p++) {
				lowerVarBoundsWeek << 2;
				upperVarBoundsWeek << 2;
			}
			x = SetVarArray(*this, varsCount);
			for (int i = 0; i < varsCount; i++) {
				x[i] = SetVar(*this, IntSet::empty, IntSet(0, teams - 1), 2, 2);
			}
			Matrix<SetVarArray> m(x, weeks, periods);
			for (int w = 0; w < weeks; w++) {
				countSet(*this, m.col(w), vals, lowerValBoundsWeek, upperValBoundsWeek,
								 lowerVarBoundsWeek, upperVarBoundsWeek, opt.ipl()); 
			}
			for (int p = 0; p < periods; p++) {
				countSet(*this, m.row(p), vals, lowerValBoundsPeriod, 
								 upperValBoundsPeriod, lowerVarBoundsPeriod, 
								 upperVarBoundsPeriod, opt.ipl()); 
			}
			for (int i = 0; i < varsCount - 1; i++) {
				for (int j = i + 1; j < varsCount; j++)
				rel(*this, x[i], SRT_NQ, x[j]);
			}
			branch(*this, x, simpleBranchVar, simpleBranchVal);
			break;
		}
		case MODEL_COUNT: {
			int varsCount = weeks * periods * 2;

			IntSetArgs boundsWeek;
			for (int i = 0; i < valsCount; i++) {
				boundsWeek << IntSet(lowerValBoundsWeek[i], upperValBoundsWeek[i]);
			}

			IntSetArgs boundsPeriod;
			for (int i = 0; i < valsCount; i++) {
				boundsPeriod << IntSet(lowerValBoundsPeriod[i], upperValBoundsPeriod[i]);
			}

			y = IntVarArray(*this, varsCount);
			for (int i = 0; i < varsCount; i++) {
				y[i] = IntVar(*this, 0, teams - 1);
			}
			Matrix<IntVarArray> m(y, weeks * 2, periods);
			for (int w = 0; w < weeks * 2; w += 2) {
				count(*this, m.col(w) + m.col(w+1), boundsWeek, vals, opt.ipl()); 
			}
			for (int p = 0; p < periods; p++) {
				count(*this, m.row(p), boundsPeriod, vals, opt.ipl()); 
			}
			for (int i = 0; i < varsCount - 2; i += 2) {
				for (int j = i + 2; j < varsCount; j += 2) {
					BoolVar tmp1(*this, 0, 1);
					BoolVar tmp2(*this, 0, 1);
					BoolVar tmp3(*this, 0, 1);
					BoolVar tmp4(*this, 0, 1);
					rel(*this, y[i], IRT_NQ, y[j], tmp1, opt.ipl());
					rel(*this, y[i + 1], IRT_NQ, y[j + 1], tmp2, opt.ipl());
					rel(*this, y[i], IRT_NQ, y[j + 1], tmp3, opt.ipl());
					rel(*this, y[i + 1], IRT_NQ, y[j], tmp4, opt.ipl());
					rel(*this, tmp1, BOT_OR, tmp2, 1, opt.ipl());
					rel(*this, tmp3, BOT_OR, tmp4, 1, opt.ipl());
				}
			}
			branch(*this, y, INT_VAR_SIZE_MIN(), INT_VAL_MAX());
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
