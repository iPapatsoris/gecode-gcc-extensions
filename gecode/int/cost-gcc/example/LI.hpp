#ifndef H_LI
#define H_LI

#include <gecode/kernel.hh>

using namespace Gecode;

class LI : public LocalHandle {
protected:
  class LIO : public LocalObject {
  public:
    int* data;
    int n;
    LIO(Space& home, int n0)
     : LocalObject(home), data(heap.alloc<int>(n0)), n(n0) {
      home.notice(*this,AP_DISPOSE);
    }
    LIO(Space& home, LIO& l)
      : LocalObject(home,l), data(heap.alloc<int>(l.n)), n(l.n) {
				heap.copy(data, l.data, l.n);
			}
    virtual LocalObject* copy(Space& home) {
      return new (home) LIO(home,*this);
    }
    virtual size_t dispose(Space& home) {
      home.ignore(*this,AP_DISPOSE);
      heap.free<int>(data,n);
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
  int operator [](int i) const {
    return static_cast<const LIO*>(object())->data[i];
  }
  int& operator [](int i) {
    return static_cast<LIO*>(object())->data[i];
  }
};

#endif