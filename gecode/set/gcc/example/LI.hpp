#ifndef H_LI
#define H_LI

#include <gecode/kernel.hh>
#include <vector>
#include <unordered_set>

using namespace Gecode;
typedef unordered_set<int> BestValsType;

/*struct BestValsWrapper {
	bool emptySet;

}*/

// Local Object Handle for an array of integers on the heap
class LI : public LocalHandle {
protected:
  class LIO : public LocalObject {
  public:
    BestValsType* data;
    int n;
    LIO(Space& home, int n0)
     : LocalObject(home), data(heap.alloc<BestValsType>(n0)), n(n0) {
      home.notice(*this,AP_DISPOSE);
    }
    LIO(Space& home, LIO& l)
      : LocalObject(home,l), data(heap.alloc<BestValsType>(l.n)), n(l.n) {
				heap.copy(data, l.data, l.n);
			}
    virtual LocalObject* copy(Space& home) {
      return new (home) LIO(home,*this);
    }
    virtual size_t dispose(Space& home) {
			//cout << *(*data).find(0) << endl;
      home.ignore(*this,AP_DISPOSE);
      heap.free<BestValsType>(data,n);
      return sizeof(*this);
    }
  };
public:
  LI(Space& home, int n) 
    : LocalHandle(new (home) LIO(home,n)) {}
  LI(const LI& li)
    : LocalHandle(li) {}
	LI() {}
  LI& operator =(const LI& li) {
    return static_cast<LI&>(LocalHandle::operator =(li));
  }
  BestValsType operator [](int i) const {
    return static_cast<const LIO*>(object())->data[i];
  }
  BestValsType& operator [](int i) {
    return static_cast<LIO*>(object())->data[i];
  }
};

#endif