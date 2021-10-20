#include <string>
#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

using namespace Gecode;

bool readInput(std::string fileName, int &vars, IntSetArgs &domain, 
							 IntArgs &vals, IntArgs &lowerBounds, IntArgs &upperBounds, 
							 IntArgs &lowerVarBounds, IntArgs &upperVarBounds);