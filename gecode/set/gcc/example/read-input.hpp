#include <string>
#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include <vector>
#include <unordered_set>

using namespace Gecode;

void readInput(std::string fileName, int &vars, std::vector<std::unordered_set<int> >& domain, 
							 IntArgs &vals, IntArgs &lowerValBounds, IntArgs &upperValBounds, 
							 IntArgs &lowerVarBounds, IntArgs &upperVarBounds);