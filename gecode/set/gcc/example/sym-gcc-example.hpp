#ifndef H_GCC_EXAMPLE
#define H_GCC_EXAMPLE

#include <gecode/driver.hh>
#include <gecode/set.hh>
#include <gecode/minimodel.hh>
//#include "LI.hpp"
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

class SymmetricGccExample : public Script {
protected:
	SetVarArray x;

public:
	enum {
		MODEL_SINGLE, MODEL_MULTI
	};

	SymmetricGccExample(const FileOptions& opt);
	SymmetricGccExample(SymmetricGccExample &s) : Script(s) {
		x.update(*this, s.x);
	}
	virtual Space *copy(void) {
		return new SymmetricGccExample(*this);
	}
	void print(ostream& os) const {
		os << "\tSolution: " << x << "\n";
	}
};

#endif