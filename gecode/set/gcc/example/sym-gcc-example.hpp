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
	IntVarArray y;
	BoolVarArray z;
	int varsCount;
	int valsCount;


public:
	enum Model {
		MODEL_SINGLE, MODEL_COUNT, MODEL_LINEAR
	};

	Model model;

	SymmetricGccExample(const FileOptions& opt);
	SymmetricGccExample(SymmetricGccExample &s) : Script(s) {
		model = s.model;
		varsCount = s.varsCount;
		valsCount = s.valsCount;
		switch (model) {
			case MODEL_SINGLE:
				x.update(*this, s.x);
				break;
			case MODEL_COUNT:
				y.update(*this, s.y);
				break;
			case MODEL_LINEAR:
				z.update(*this, s.z);
				break;
			}
	}
	virtual Space *copy(void) {
		return new SymmetricGccExample(*this);
	}
	void print(ostream& os) const {
		return;
		switch (model) {
			case MODEL_SINGLE:
				os << "\tSolution: " << x << "\n";
				break;
			case MODEL_COUNT:
			case MODEL_LINEAR:
				os << "\tSolution: \n";
				for (int i = 0; i < varsCount; i++) {
					for (int j = 0; j < valsCount; j++) {
						if (model == MODEL_COUNT) {
							cout << y[valsCount*i + j] << " ";
						} else if (model == MODEL_LINEAR) {
							cout << z[valsCount*i + j] << " ";
						}
					}
					cout << "\n";
				}
				break;
		}
	}
};

#endif