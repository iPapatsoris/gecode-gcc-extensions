/*
 *  Authors:
 *    Christian Schulte <schulte@gecode.org>
 *
 *  Copyright:
 *    Christian Schulte, 2008-2019
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software, to deal in the software without restriction,
 *  including without limitation the rights to use, copy, modify, merge,
 *  publish, distribute, sublicense, and/or sell copies of the software,
 *  and to permit persons to whom the software is furnished to do so, subject
 *  to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <gecode/set.hh>
#include "BestBranch.hpp"

using namespace Gecode;

class BestVal : public Brancher {
protected:
  ViewArray<Set::SetView> x;
  mutable int start;
	BestBranch bestBranch;
  class PosVal : public Choice {
  public:
    int pos; BestValsType val;
		bool isEmptySetVal;
    PosVal(const BestVal& b, int p, BestValsType v, bool emptySet)
      : Choice(b,2), pos(p), val(v), isEmptySetVal(emptySet) {}
    virtual void archive(Archive& e) const {
      Choice::archive(e);
      //e << pos << val;
    }
  };
public:
  BestVal(Home home, ViewArray<Set::SetView>& x0, BestBranch& bestBranch)
    : Brancher(home), x(x0), start(0), bestBranch(bestBranch) {}
  static void post(Home home, ViewArray<Set::SetView>& x, BestBranch& bestBranch) {
    (void) new (home) BestVal(home,x,bestBranch);
  }
  virtual size_t dispose(Space& home) {
    (void) Brancher::dispose(home);
    return sizeof(*this);
  }
  BestVal(Space& home, BestVal& b)
    : Brancher(home,b), start(b.start) {
    x.update(home,b.x);
		bestBranch.update(home,b.bestBranch);
		//cout << "BestVal copy: " << b.bestBranch[0] << endl;
  }
  virtual Brancher* copy(Space& home) {
    return new (home) BestVal(home,*this);
  }
  virtual bool status(const Space&) const {
    for (int i=start; i<x.size(); i++)
      if (!x[i].assigned()) {
        start = i; return true;
      }
    return false;
  }
  virtual Choice* choice(Space&) {
    int p = start;
    // unsigned int s = x[p].unknownSize();
    // for (int i=start+1; i<x.size(); i++)
    //   if (!x[i].assigned() /*&& (x[i].unknownSize() < s*/) {
    //     p = i;// s = x[p].unknownSize();
    //   }

		return new PosVal(*this, p, bestBranch[p], false);

		/*for (auto& v: bestBranch[p]) {
			cout << v << endl;
			if (!x[p].contains(v) && !x[p].notContains(v)) {
				cout << "branch on var " << p << " val " << v << endl;  
				return new PosVal(*this, p, v, false);
			}
		}
		return new PosVal(*this, p, 0, true);*/
		//for (SetVarUnknownValues i(x[p]); i(); ++i) {
		//}    
  }
  virtual Choice* choice(const Space&, Archive& e) {
    int pos;
		BestValsType val;
		bool isEmptySetVal;
    //e >> pos >> val >> isEmptySetVal;
    return new PosVal(*this, pos, val, isEmptySetVal);
  }
  virtual ExecStatus commit(Space& home, 
                            const Choice& c,
                            unsigned int a) {
    const PosVal& pv = static_cast<const PosVal&>(c);
    int pos=pv.pos;

			cout << "upper bound of var " << pos << ":";
			for (SetVarLubValues i(x[pos]); i(); ++i)
				std::cout << i.val() << " ";
			cout << "cardinalities " << x[pos].cardMin() << " " << x[pos].cardMax() << endl;
			cout << "lower bound of var " << pos << ":";
			for (SetVarGlbValues i(x[pos]); i(); ++i)
	std::cout << i.val() << " " 	;
		if (!a) {
			cout << "Including to var " << pos << " set ";
			for (auto v: pv.val) {
				cout << v << " ";
			}
			cout << endl;
			for (auto v: pv.val) {
				if (me_failed(x[pos].include(home, v))) {
					cout << "Inclusion failed" << endl;
					return ES_FAILED;
				}
			}
			auto res = me_failed(x[pos].cardMax(home, x[pos].glbSize())); 
				cout << "cardinalities of var " << pos << " now " << x[pos].cardMin() << " " << x[pos].cardMax() << endl;
		return  res ? ES_FAILED : ES_OK;
		} else {
			cout << "Excluding from var " << pos << " set ";
			for (auto v: pv.val) {
				cout << v << " ";
			}
			cout << endl;
			
			if (!pv.val.size()) {
				return me_failed(x[pos].cardMin(home, 1)) ? ES_FAILED : ES_OK;
			}	
			for (auto v: pv.val) {
				if (me_failed(x[pos].exclude(home, v))) {
				cout << "trying to exclude val " << v << " from var " << pos << " failed " << endl;
				return ES_FAILED;
				}
			}
		}
		return ES_OK;

		if (a < 2) {
			if (a == 1 && !pv.val.size()) {
				return ES_NOFIX;
			} 
			cout << "Including to var " << pos << " set ";
			for (auto v: pv.val) {
				cout << v << " ";
			}
			cout << endl;
			for (auto v: pv.val) {
				if (me_failed(x[pos].include(home, v))) {
					cout << "Inclusion failed" << endl;
					return ES_FAILED;
				}
			}
			// if (me_failed(x[pos].includeI(home, pv.val))) {
			// 	return ES_FAILED; 
			// }
			if (!a) {
				cout << "Assigning var " << pos << " max cardinality to " << x[pos].glbSize() << endl; 
				auto res = me_failed(x[pos].cardMax(home, x[pos].glbSize())); 
				return  res ? ES_FAILED : ES_OK;
			} else {
				cout << "Increasing var " << pos << " min cardinality to " << x[pos].cardMin() + 1 << endl;
				auto res =	me_failed(x[pos].cardMin(home, x[pos].cardMin() + 1));
				return res ? ES_FAILED : ES_OK;
			}
		} else {
			cout << "Excluding from var " << pos << " set ";
			for (auto v: pv.val) {
				cout << v << " ";
			}
			cout << endl;
			for (auto v: pv.val) {
				if (me_failed(x[pos].exclude(home, v))) {
					cout << "exclusion failed " << endl;
					return ES_FAILED;
				}
			}
			// if (me_failed(x[pos].excludeI(home, pv.val))) {
			// 	return ES_FAILED; 
			// }
		}

		return ES_OK;
		
		// if (!pv.val.size()) {
		// 	if (a == 0)
		// 		return me_failed(x[pos].cardMax(home, 0)) ? ES_FAILED : ES_OK;
		// 	else if (a == 1)
		// 		return me_failed(x[pos].cardMin(home, 1)) ? ES_FAILED : ES_OK;
		// 	else 
		// 		return ES_FAILED;
		// } else {
		// 	int v;
		// 	for (SetVarUnknownValues i(x[pos]); i(); ++i) {
		// 		if (pv.val.find(i.val()) != pv.val.end()) {
		// 			v = i.val();
		// 			break;
		// 		}
		// 	}

		// 	if (a < 2) {
		// 		if (me_failed(x[pos].include(home, v))) {
		// 			return ES_FAILED;
		// 		}
		// 	} else {
		// 		if (me_failed(x[pos].exclude(home, v))) {
		// 			return ES_FAILED;
		// 		}
		// 	}
		// 	if (a == 0) {
		// 		return me_failed(x[pos].cardMax(home, x[pos].glbSize())) ? ES_FAILED : ES_OK;
		// 	} else if (a == 1) {
		// 		if (x[pos].cardMin() == x[pos].glbSize()) {
		// 			return me_failed(x[pos].cardMin(home, x[pos].cardMin() + 1)) ? ES_FAILED : ES_OK;
		// 		}
		// 	}
		// 	return ES_OK;
		// }
  }

  virtual void print(const Space&, const Choice& c,
                     unsigned int a,
                     std::ostream& o) const {
    const PosVal& pv = static_cast<const PosVal&>(c);
    int pos=pv.pos;
    if (a == 0)
      o << "x[" << pos << "] = ";// << pv.val;
    else
      o << "x[" << pos << "] != ";// << pv.val;
  }
};

void bestval(Home home, const SetVarArgs& x, BestBranch& bestBranch) {
  if (home.failed()) return;
  ViewArray<Set::SetView> y(home,x);
  BestVal::post(home,y,bestBranch);
}