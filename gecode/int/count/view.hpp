/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Christian Schulte <schulte@gecode.org>
 *
 *  Copyright:
 *     Christian Schulte, 2003
 *
 *  Last modified:
 *     $Date$ by $Author$
 *     $Revision$
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
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

namespace Gecode { namespace Int { namespace Count {

  /*
   * General baseclass
   *
   */

  template<class VX, class VY, class VZ>
  forceinline
  BaseView<VX,VY,VZ>::BaseView(Home home,
                               ViewArray<VX>& x0, VY y0, VZ z0, int c0)
    : Propagator(home), x(x0), y(y0), z(z0), c(c0) {
    if (vtd(y) == VTD_INTSET)
      home.notice(*this,AP_DISPOSE);
    x.subscribe(home,*this,PC_INT_DOM);
    subscribe(home,*this,y);
    z.subscribe(home,*this,PC_INT_BND);
  }

  template<class VX, class VY, class VZ>
  forceinline
  BaseView<VX,VY,VZ>::BaseView(Space& home, bool share, BaseView<VX,VY,VZ>& p)
    : Propagator(home,share,p), c(p.c) {
    x.update(home,share,p.x);
    y.update(home,share,p.y);
    z.update(home,share,p.z);
  }

  template<class VX, class VY, class VZ>
  PropCost
  BaseView<VX,VY,VZ>::cost(const Space&, const ModEventDelta&) const {
    return PropCost::linear(PropCost::LO,x.size()+1);
  }

  template<class VX, class VY, class VZ>
  forceinline size_t
  BaseView<VX,VY,VZ>::dispose(Space& home) {
    if (vtd(y) == VTD_INTSET)
      home.ignore(*this,AP_DISPOSE);
    x.cancel(home,*this,PC_INT_DOM);
    cancel(home,*this,y);
    z.cancel(home,*this,PC_INT_BND);
    (void) Propagator::dispose(home);
    return sizeof(*this);
  }

  template<class VX, class VY, class VZ>
  forceinline void
  BaseView<VX,VY,VZ>::count(Space& home) {
    int n = x.size();
    for (int i=n; i--; )
      switch (holds(x[i],y)) {
      case RT_FALSE:
        x[i].cancel(home,*this,PC_INT_DOM); x[i]=x[--n];
        break;
      case RT_TRUE:
        x[i].cancel(home,*this,PC_INT_DOM); x[i]=x[--n];
        c--;
        break;
      case RT_MAYBE:
        break;
      default:
        GECODE_NEVER;
      }
    x.size(n);
  }

  template<class VX, class VY, class VZ>
  forceinline int
  BaseView<VX,VY,VZ>::atleast(void) const {
    return -c;
  }

  template<class VX, class VY, class VZ>
  forceinline int
  BaseView<VX,VY,VZ>::atmost(void) const {
    return x.size()-c;
  }

  template<class VX>
  forceinline bool
  shared(const IntSet&, VX) {
    return false;
  }
  template<class VX, class VY, class VZ>
  forceinline bool
  BaseView<VX,VY,VZ>::sharing(const ViewArray<VX>& x, 
                              const VY& y, const VZ& z) {
    if (shared(y,z))
      return true;
    for (int i = x.size(); i--; )
      if (shared(x[i],z))
        return true;
    return false;
  }

  /*
   * Equality
   *
   */

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  EqView<VX,VY,VZ,shr>::EqView(Home home,
                               ViewArray<VX>& x, VY y, VZ z, int c)
    : BaseView<VX,VY,VZ>(home,x,y,z,c) {}

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  EqView<VX,VY,VZ,shr>::post(Home home,
                             ViewArray<VX>& x, VY y, VZ z, int c) {
    if ((vtd(y) != VTD_VARVIEW) && z.assigned())
      return EqInt<VX,VY>::post(home,x,y,z.val()+c);
    if (sharing(x,y,z))
      (void) new (home) EqView<VX,VY,VZ,true>(home,x,y,z,c);
    else
      (void) new (home) EqView<VX,VY,VZ,false>(home,x,y,z,c);
    return ES_OK;
  }

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  EqView<VX,VY,VZ,shr>::EqView(Space& home, bool share,
                               EqView<VX,VY,VZ,shr>& p)
    : BaseView<VX,VY,VZ>(home,share,p) {}

  template<class VX, class VY, class VZ, bool shr>
  Actor*
  EqView<VX,VY,VZ,shr>::copy(Space& home, bool share) {
    return new (home) EqView<VX,VY,VZ,shr>(home,share,*this);
  }

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  EqView<VX,VY,VZ,shr>::propagate(Space& home, const ModEventDelta&) {
    count(home);

    GECODE_ME_CHECK(z.gq(home,atleast()));
    GECODE_ME_CHECK(z.lq(home,atmost()));

    if (z.assigned()) {
      if (z.val() == atleast()) {
        GECODE_ES_CHECK(post_false(home,x,y));
        return home.ES_SUBSUMED(*this);
      }
      if (z.val() == atmost()) {
        GECODE_ES_CHECK(post_true(home,x,y));
        return home.ES_SUBSUMED(*this);
      }
      if (vtd(y) != VTD_VARVIEW) {
        VY yc(y);
        GECODE_REWRITE(*this,(EqInt<VX,VY>
                              ::post(home(*this),x,yc,z.val()+c)));    
      }
    }


    if ((vtd(y) == VTD_VARVIEW) && (z.min() > 0)) {
      /* 
       * Only if the propagator is at fixpoint here, continue
       * when things are shared: the reason is that prune
       * requires that the views in x overlap with y!
       */
      if (shr && (VX::me(Propagator::modeventdelta()) != ME_INT_NONE))
        return ES_NOFIX;

      GECODE_ES_CHECK(prune(home,x,y));

      return ES_NOFIX;
    }
    
    return shr ? ES_NOFIX : ES_FIX;
  }




  /*
   * Disequality
   *
   */

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  NqView<VX,VY,VZ,shr>::NqView(Home home,
                               ViewArray<VX>& x, VY y, VZ z, int c)
    : BaseView<VX,VY,VZ>(home,x,y,z,c) {}

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  NqView<VX,VY,VZ,shr>::post(Home home,
                             ViewArray<VX>& x, VY y, VZ z, int c) {
    if (z.assigned())
      return NqInt<VX,VY>::post(home,x,y,z.val()+c);
    (void) new (home) NqView<VX,VY,VZ,shr>(home,x,y,z,c);
    return ES_OK;
  }

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  NqView<VX,VY,VZ,shr>::NqView(Space& home, bool share,
                               NqView<VX,VY,VZ,shr>& p)
    : BaseView<VX,VY,VZ>(home,share,p) {}

  template<class VX, class VY, class VZ, bool shr>
  Actor*
  NqView<VX,VY,VZ,shr>::copy(Space& home, bool share) {
    return new (home) NqView<VX,VY,VZ,shr>(home,share,*this);
  }

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  NqView<VX,VY,VZ,shr>::propagate(Space& home, const ModEventDelta&) {
    count(home);
    if (atleast() == atmost()) {
      GECODE_ME_CHECK(z.nq(home,atleast()));
      return home.ES_SUBSUMED(*this);
    }
    if (z.max() < atleast())
      return home.ES_SUBSUMED(*this);
    if (z.min() > atmost())
      return home.ES_SUBSUMED(*this);
    if (z.assigned()) {
      VY yc(y);
      GECODE_REWRITE(*this,(NqInt<VX,VY>::post(home(*this),x,yc,z.val()+c)));
    }

    return ES_FIX;
  }



  /*
   * Less or equal
   *
   */

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  LqView<VX,VY,VZ,shr>::LqView(Home home, ViewArray<VX>& x,
                               VY y, VZ z, int c)
    : BaseView<VX,VY,VZ>(home,x,y,z,c) {}

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  LqView<VX,VY,VZ,shr>::post(Home home, ViewArray<VX>& x,
                             VY y, VZ z, int c) {
    if (z.assigned())
      return LqInt<VX,VY>::post(home,x,y,z.val()+c);
    if (sharing(x,y,z))
      (void) new (home) LqView<VX,VY,VZ,true>(home,x,y,z,c);
    else
      (void) new (home) LqView<VX,VY,VZ,false>(home,x,y,z,c);
    return ES_OK;
  }

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  LqView<VX,VY,VZ,shr>::LqView(Space& home, bool share,
                               LqView<VX,VY,VZ,shr>& p)
    : BaseView<VX,VY,VZ>(home,share,p) {}

  template<class VX, class VY, class VZ, bool shr>
  Actor*
  LqView<VX,VY,VZ,shr>::copy(Space& home, bool share) {
    return new (home) LqView<VX,VY,VZ,shr>(home,share,*this);
  }

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  LqView<VX,VY,VZ,shr>::propagate(Space& home, const ModEventDelta&) {
    count(home);
    GECODE_ME_CHECK(z.gq(home,atleast()));

    if (z.max() == atleast()) {
      GECODE_ES_CHECK(post_false(home,x,y));
      return home.ES_SUBSUMED(*this);
    }

    if (x.size() == 0)
      return home.ES_SUBSUMED(*this);

    if (z.assigned()) {
      VY yc(y);
      GECODE_REWRITE(*this,(LqInt<VX,VY>::post(home(*this),x,yc,z.val()+c)));
    }

    return shr ? ES_NOFIX : ES_FIX;
  }



  /*
   * Greater or equal
   *
   */

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  GqView<VX,VY,VZ,shr>::GqView(Home home, ViewArray<VX>& x, VY y, VZ z, int c)
    : BaseView<VX,VY,VZ>(home,x,y,z,c) {}

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  GqView<VX,VY,VZ,shr>::post(Home home,
                             ViewArray<VX>& x, VY y, VZ z, int c) {
    if ((vtd(y) != VTD_VARVIEW) && z.assigned())
      return GqInt<VX,VY>::post(home,x,y,z.val()+c);
    if (sharing(x,y,z))
      (void) new (home) GqView<VX,VY,VZ,true>(home,x,y,z,c);
    else
      (void) new (home) GqView<VX,VY,VZ,false>(home,x,y,z,c);
    return ES_OK;
  }

  template<class VX, class VY, class VZ, bool shr>
  forceinline
  GqView<VX,VY,VZ,shr>::GqView(Space& home, bool share,
                               GqView<VX,VY,VZ,shr>& p)
    : BaseView<VX,VY,VZ>(home,share,p) {}

  template<class VX, class VY, class VZ, bool shr>
  Actor*
  GqView<VX,VY,VZ,shr>::copy(Space& home, bool share) {
    return new (home) GqView<VX,VY,VZ,shr>(home,share,*this);
  }

  template<class VX, class VY, class VZ, bool shr>
  ExecStatus
  GqView<VX,VY,VZ,shr>::propagate(Space& home, const ModEventDelta&) {
    count(home);

    GECODE_ME_CHECK(z.lq(home,atmost()));

    if (z.min() == atmost()) {
      GECODE_ES_CHECK(post_true(home,x,y)); 
      return home.ES_SUBSUMED(*this);
    }
    if (x.size() == 0)
      return home.ES_SUBSUMED(*this);

    if (z.assigned() && (vtd(y) != VTD_VARVIEW)) {
      VY yc(y);
      GECODE_REWRITE(*this,(GqInt<VX,VY>::post(home(*this),x,yc,z.val()+c)));
    }

    if ((vtd(y) == VTD_VARVIEW) && (z.min() > 0)) {
      /* 
       * Only if the propagator is at fixpoint here, continue
       * when things are shared: the reason is that prune
       * requires that the views in x overlap with y!
       */
      if (shr && (VX::me(Propagator::modeventdelta()) != ME_INT_NONE))
        return ES_NOFIX;

      GECODE_ES_CHECK(prune(home,x,y));

      return ES_NOFIX;
    }
    
    return shr ? ES_NOFIX : ES_FIX;
  }

}}}

// STATISTICS: int-prop

