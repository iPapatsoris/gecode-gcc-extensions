/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Christian Schulte <schulte@gecode.org>
 *
 *  Copyright:
 *     Christian Schulte, 2012
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

#include <gecode/int/branch.hh>

namespace Gecode { namespace Int { namespace Branch {

  ViewSel<IntView>*
  viewselint(Space& home, const IntVarBranch& ivb) {
    switch (ivb.select()) {
    case IntVarBranch::SEL_NONE:
      return new (home) ViewSelNone<IntView>(home,ivb);
    case IntVarBranch::SEL_RND:
      return new (home) ViewSelRnd<IntView>(home,ivb);
    default: break;
    }
    if (ivb.tbl() != NULL) {
      switch (ivb.select()) {
      case IntVarBranch::SEL_MERIT_MIN:
        return new (home) ViewSelMinTbl<MeritFunction<IntView> >(home,ivb);
      case IntVarBranch::SEL_MERIT_MAX:
        return new (home) ViewSelMaxTbl<MeritFunction<IntView> >(home,ivb);
      case IntVarBranch::SEL_MIN_MIN:
        return new (home) ViewSelMinTbl<MeritMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_MIN_MAX:
        return new (home) ViewSelMaxTbl<MeritMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_MAX_MIN:
        return new (home) ViewSelMinTbl<MeritMax<IntView> >(home,ivb);
      case IntVarBranch::SEL_MAX_MAX:
        return new (home) ViewSelMaxTbl<MeritMax<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_MIN:
        return new (home) ViewSelMinTbl<MeritSize<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_MAX:
        return new (home) ViewSelMaxTbl<MeritSize<IntView> >(home,ivb);
      case IntVarBranch::SEL_DEGREE_MIN:
        return new (home) ViewSelMinTbl<MeritDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_DEGREE_MAX:
        return new (home) ViewSelMaxTbl<MeritDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MIN:
        return new (home) ViewSelMinTbl<MeritAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MAX:
        return new (home) ViewSelMaxTbl<MeritAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MIN:
        return new (home) ViewSelMinTbl<MeritActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MAX:
        return new (home) ViewSelMaxTbl<MeritActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_DEGREE_MIN:
        return new (home) ViewSelMinTbl<MeritSizeDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_DEGREE_MAX:
        return new (home) ViewSelMaxTbl<MeritSizeDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_AFC_MIN:
        return new (home) ViewSelMinTbl<MeritSizeAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_AFC_MAX:
        return new (home) ViewSelMaxTbl<MeritSizeAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_ACTIVITY_MIN:
        return new (home) ViewSelMinTbl<MeritSizeActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_ACTIVITY_MAX:
        return new (home) ViewSelMaxTbl<MeritSizeActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MIN_MIN:
        return new (home) ViewSelMinTbl<MeritRegretMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MIN_MAX:
        return new (home) ViewSelMaxTbl<MeritRegretMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MAX_MIN:
        return new (home) ViewSelMinTbl<MeritRegretMax<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MAX_MAX:
        return new (home) ViewSelMaxTbl<MeritRegretMax<IntView> >(home,ivb);
      default:
        throw UnknownBranching("Int::branch");
      }
    } else {
      switch (ivb.select()) {
      case IntVarBranch::SEL_MERIT_MIN:
        return new (home) ViewSelMin<MeritFunction<IntView> >(home,ivb);
      case IntVarBranch::SEL_MERIT_MAX:
        return new (home) ViewSelMax<MeritFunction<IntView> >(home,ivb);
      case IntVarBranch::SEL_MIN_MIN:
        return new (home) ViewSelMin<MeritMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_MIN_MAX:
        return new (home) ViewSelMax<MeritMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_MAX_MIN:
        return new (home) ViewSelMin<MeritMax<IntView> >(home,ivb);
      case IntVarBranch::SEL_MAX_MAX:
        return new (home) ViewSelMax<MeritMax<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_MIN:
        return new (home) ViewSelMin<MeritSize<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_MAX:
        return new (home) ViewSelMax<MeritSize<IntView> >(home,ivb);
      case IntVarBranch::SEL_DEGREE_MIN:
        return new (home) ViewSelMin<MeritDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_DEGREE_MAX:
        return new (home) ViewSelMax<MeritDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MIN:
        return new (home) ViewSelMin<MeritAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MAX:
        return new (home) ViewSelMax<MeritAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MIN:
        return new (home) ViewSelMin<MeritActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MAX:
        return new (home) ViewSelMax<MeritActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_DEGREE_MIN:
        return new (home) ViewSelMin<MeritSizeDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_DEGREE_MAX:
        return new (home) ViewSelMax<MeritSizeDegree<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_AFC_MIN:
        return new (home) ViewSelMin<MeritSizeAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_AFC_MAX:
        return new (home) ViewSelMax<MeritSizeAfc<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_ACTIVITY_MIN:
        return new (home) ViewSelMin<MeritSizeActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_ACTIVITY_MAX:
        return new (home) ViewSelMax<MeritSizeActivity<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MIN_MIN:
        return new (home) ViewSelMin<MeritRegretMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MIN_MAX:
        return new (home) ViewSelMax<MeritRegretMin<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MAX_MIN:
        return new (home) ViewSelMin<MeritRegretMax<IntView> >(home,ivb);
      case IntVarBranch::SEL_REGRET_MAX_MAX:
        return new (home) ViewSelMax<MeritRegretMax<IntView> >(home,ivb);
      default:
        throw UnknownBranching("Int::branch");
      }
    }
    GECODE_NEVER;
    return NULL;
  }

  ViewSel<BoolView>*
  viewselbool(Space& home, const IntVarBranch& ivb) {
    switch (ivb.select()) {
    case IntVarBranch::SEL_NONE:
      return new (home) ViewSelNone<BoolView>(home,ivb);
    case IntVarBranch::SEL_RND:
      return new (home) ViewSelRnd<BoolView>(home,ivb);
    default: break;
    }
    if (ivb.tbl() != NULL) {
      switch (ivb.select()) {
      case IntVarBranch::SEL_MERIT_MIN:
        return new (home) ViewSelMinTbl<MeritFunction<BoolView> >(home,ivb);
      case IntVarBranch::SEL_MERIT_MAX:
        return new (home) ViewSelMaxTbl<MeritFunction<BoolView> >(home,ivb);
      case IntVarBranch::SEL_MIN_MIN:
      case IntVarBranch::SEL_MIN_MAX:
      case IntVarBranch::SEL_MAX_MIN:
      case IntVarBranch::SEL_MAX_MAX:
      case IntVarBranch::SEL_SIZE_MIN:
      case IntVarBranch::SEL_SIZE_MAX:
      case IntVarBranch::SEL_REGRET_MIN_MIN:
      case IntVarBranch::SEL_REGRET_MIN_MAX:
      case IntVarBranch::SEL_REGRET_MAX_MIN:
      case IntVarBranch::SEL_REGRET_MAX_MAX:
        return new (home) ViewSelNone<BoolView>(home,ivb);
      case IntVarBranch::SEL_DEGREE_MIN:
        return new (home) ViewSelMinTbl<MeritDegree<BoolView> >(home,ivb);
      case IntVarBranch::SEL_DEGREE_MAX:
        return new (home) ViewSelMaxTbl<MeritDegree<BoolView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MIN:
        return new (home) ViewSelMinTbl<MeritAfc<BoolView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MAX:
        return new (home) ViewSelMaxTbl<MeritAfc<BoolView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MIN:
        return new (home) ViewSelMinTbl<MeritActivity<BoolView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MAX:
        return new (home) ViewSelMaxTbl<MeritActivity<BoolView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_DEGREE_MIN:
        return new (home) ViewSelMinTbl<MeritSizeDegree<BoolView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_DEGREE_MAX:
        return new (home) ViewSelMaxTbl<MeritSizeDegree<BoolView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_AFC_MIN:
        return new (home) ViewSelMinTbl<MeritSizeAfc<BoolView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_AFC_MAX:
        return new (home) ViewSelMaxTbl<MeritSizeAfc<BoolView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_ACTIVITY_MIN:
        return new (home) ViewSelMinTbl<MeritSizeActivity<BoolView> >(home,ivb);
      case IntVarBranch::SEL_SIZE_ACTIVITY_MAX:
        return new (home) ViewSelMaxTbl<MeritSizeActivity<BoolView> >(home,ivb);
      default:
        throw UnknownBranching("Int::branch");
      }
    } else {
      switch (ivb.select()) {
      case IntVarBranch::SEL_MERIT_MIN:
        return new (home) ViewSelMin<MeritFunction<BoolView> >(home,ivb);
      case IntVarBranch::SEL_MERIT_MAX:
        return new (home) ViewSelMax<MeritFunction<BoolView> >(home,ivb);
      case IntVarBranch::SEL_MIN_MIN:
      case IntVarBranch::SEL_MIN_MAX:
      case IntVarBranch::SEL_MAX_MIN:
      case IntVarBranch::SEL_MAX_MAX:
      case IntVarBranch::SEL_SIZE_MIN:
      case IntVarBranch::SEL_SIZE_MAX:
      case IntVarBranch::SEL_REGRET_MIN_MIN:
      case IntVarBranch::SEL_REGRET_MIN_MAX:
      case IntVarBranch::SEL_REGRET_MAX_MIN:
      case IntVarBranch::SEL_REGRET_MAX_MAX:
        return new (home) ViewSelNone<BoolView>(home,ivb);
      case IntVarBranch::SEL_DEGREE_MIN:
      case IntVarBranch::SEL_SIZE_DEGREE_MAX:
        return new (home) ViewSelMin<MeritDegree<BoolView> >(home,ivb);
      case IntVarBranch::SEL_DEGREE_MAX:
      case IntVarBranch::SEL_SIZE_DEGREE_MIN:
        return new (home) ViewSelMax<MeritDegree<BoolView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MIN:
      case IntVarBranch::SEL_SIZE_AFC_MAX:
        return new (home) ViewSelMin<MeritAfc<BoolView> >(home,ivb);
      case IntVarBranch::SEL_AFC_MAX:
      case IntVarBranch::SEL_SIZE_AFC_MIN:
        return new (home) ViewSelMax<MeritAfc<BoolView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MIN:
      case IntVarBranch::SEL_SIZE_ACTIVITY_MAX:
        return new (home) ViewSelMin<MeritActivity<BoolView> >(home,ivb);
      case IntVarBranch::SEL_ACTIVITY_MAX:
      case IntVarBranch::SEL_SIZE_ACTIVITY_MIN:
        return new (home) ViewSelMax<MeritActivity<BoolView> >(home,ivb);
      default:
        throw UnknownBranching("Int::branch");
      }
    }
    GECODE_NEVER;
    return NULL;
  }

}}}


// STATISTICS: int-branch
