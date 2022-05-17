#ifndef H_COST_GCC_EXAMPLE
#define H_COST_GCC_EXAMPLE

#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include "LI.hpp"
#include "read-input.hpp"

using namespace Gecode;
using namespace std;

class FileOptions : public Options {
protected:
	Driver::StringValueOption _file;
	Driver::BoolOption _branch;
public:
	FileOptions(const char* scriptName) : 
			Options(scriptName),
			_file("file","input file name", ""),
			_branch("branch", "branch heuristic flag", true) { 
		add(_file);
		add(_branch);
	}
  string file(void) const { return _file.value(); }
	bool branch(void) const { return _branch.value(); }
};

class CountCostsExample : public Script {
protected:
	IntVarArray x;

public:
	enum {
		MODEL_SINGLE, MODEL_MULTI
	};

	CountCostsExample(const FileOptions& opt);
	CountCostsExample(CountCostsExample &s) : Script(s) {
		x.update(*this, s.x);
	}
	virtual Space *copy(void) {
		return new CountCostsExample(*this);
	}
	void print(ostream& os) const {
		os << "\tSolution: " << x << "\n";
	}
};

#endif