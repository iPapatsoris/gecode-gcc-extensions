/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Mikael Lagerkvist <lagerkvist@gecode.org>
 *
 *  Copyright:
 *     Mikael Lagerkvist, 2005
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

#include "examples/support.hh"
#include <gecode/minimodel.hh>

namespace {

  /// List of specifications
  extern const int* specs[];
  /// Number of specifications
  extern const unsigned int n_examples;

}

/**
 * \brief %Example: Nonogram
 *
 * This example solves nonograms. A nonogram is composed of a matrix of 
 * markers. For each row/column there is a specification on how many groups
 * of markers (separated by one or more unmarked spots) and their length. 
 * The objective is to find a valid assignment, which incidentally may also 
 * produce a pretty picture.
 *
 * See problem 12 at http://www.csplib.org/.
 *
 * \ingroup ExProblem
 *
 */
class Nonogram : public Example {
protected:
  /// Specification to be used
  const int* spec;
  /// Fields of board
  BoolVarArray b;

  /// Return width of board
  int width(void) const {
    return spec[0];
  }
  /// Return height of board
  int height(void) const {
    return spec[1];
  }

  /// Access position (h,w) in matrix
  BoolVar pos(int h, int w) {
    return b[h*width() + w];
  }

  /// Returns next regular expression for line starting from spos
  DFA line(int& spos) {
    REG r0(0), r1(1);
    REG border = *r0;
    REG between = +r0;
    int size = spec[spos++];
    REG r = border + r1(spec[spos],spec[spos]);
    spos++;
    for (int i=size-1; i--; spos++)
      r += between + r1(spec[spos],spec[spos]);
    return r + border;
  }


public:
  /// Construction of the model.
  Nonogram(const SizeOptions& opt)
    : spec(specs[opt.size()]), b(*this,width()*height(),0,1) {
    int spos = 2;
    Matrix<BoolVarArray> m(b, width(), height());

    // Post constraints for columns
    for (int w = 0; w < width(); ++w)
      extensional(*this, m.col(w), line(spos));

    // Post constraints for rows
    for (int h = 0; h < height(); ++h)
      extensional(*this, m.row(h), line(spos));

    // Install branchings
    branch(*this, b, INT_VAR_NONE, INT_VAL_MAX);
  }

  /// Constructor for cloning \a s
  Nonogram(bool share, Nonogram& s) : Example(share,s), spec(s.spec) {
    b.update(*this, share, s.b);
  }

  /// Copy space during cloning
  virtual Space*
  copy(bool share) {
    return new Nonogram(share,*this);
  }

  /// Print solution
  virtual void
  print(std::ostream& os) {
    for (int h = 0; h < height(); ++h) {
      os << '\t';
      for (int w = 0; w < width(); ++w)
        os << ((pos(h,w).val() == 1) ? '#' : ' ');
      os << std::endl;
    }
    os << std::endl;
  }
};


/** \brief Main-function
 *  \relates Nonogram
 */
int
main(int argc, char* argv[]) {
  SizeOptions opt("Nonogram");
  opt.size(8);
  opt.parse(argc,argv);
  if (opt.size() >= n_examples) {
    std::cerr << "Error: size must be between 0 and "
              << n_examples-1 << std::endl;
    return 1;
  }
  Example::run<Nonogram,DFS,SizeOptions>(opt);
  return 0;
}

namespace {

  /** \name Picture specifications
   *
   * A specification is given by a list of integers. The first two
   * integers (w and h) specify the number of columns and rows
   * respectively. Then w + h groups of integers follows. Each group is
   * started by the number of integers it contains (n), followed by n integers
   * specifying the sizes of the stretches of markers in that row/column.
   *
   * \relates Nonogram
   */
  //@{
  /// Specification for a heart-shaped picture.
const int heart[] =
  { 9, 9,
    // Column constraints.
    1, 3,
    2, 2, 3,
    2, 2, 2,
    2, 2, 2,
    2, 2, 2,
    2, 2, 2,
    2, 2, 2,
    2, 2, 3,
    1, 3,
    // Row constraints
    2, 2, 2,
    2, 4, 4,
    3, 1, 3, 1,
    3, 2, 1, 2,
    2, 1, 1,
    2, 2, 2,
    2, 2, 2,
    1, 3,
    1, 1
  };

  /// Specification for a bear/bunny-shaped picture.
const int bear[] =
  { 13, 8,
    // Column constraints
    1, 2,
    2, 2, 1,
    2, 3, 2,
    1, 6,
    2, 1, 4,
    1, 3,
    1, 4,
    1, 4,
    1, 4,
    1, 5,
    1, 4,
    2, 1, 3,
    1, 2,
    // Row constraints
    1, 1,
    1, 2,
    2, 4, 4,
    1, 12,
    1, 8,
    1, 9,
    2, 3, 4,
    2, 2, 2
  };

  /// Specification for a crocodile-shaped picture
const int crocodile[] =
  { 15, 9,
    // Column constraints
    1, 3,
    1, 4,
    2, 2, 2,
    2, 3, 1,
    2, 2, 3,
    2, 3, 2,
    2, 2, 3,
    2, 4, 2,
    2, 3, 2,
    1, 6,
    2, 1, 3,
    2, 1, 3,
    2, 1, 4,
    1, 5,
    1, 5,
    // Row constraints
    1, 3,
    3, 2, 3, 2,
    2, 10, 3,
    1, 15,
    5, 1, 1, 1, 1, 6,
    2, 1, 7,
    2, 1, 4,
    2, 1, 4,
    1, 4
  };

  /// Specification for an unknown picture
const int unknown[] =
  { 10, 10,
    // Column constraints
    1, 3,
    2, 2, 1,
    2, 2, 2,
    2, 2, 1,
    3, 1, 2, 1,
    2, 1, 1,
    3, 1, 4, 1,
    3, 1, 1, 2,
    2, 3, 1,
    1, 4,
    // Row constraints
    1, 3,
    2, 2, 1,
    2, 1, 1,
    2, 1, 4,
    4, 1, 1, 1, 1,
    4, 2, 1, 1, 1,
    3, 2, 1, 1,
    2, 1, 2,
    2, 2, 3,
    1, 3
  };

  /// Specification for a pinwheel-picture
const int pinwheel[] =
  { 6, 6,
    // Column constraints
    2, 1, 2,
    1, 1,
    1, 2,
    1, 2,
    1, 1,
    2, 2, 1,
    // Row constraints
    2, 2, 1,
    1, 1,
    1, 2,
    1, 2,
    1, 1,
    2, 1, 2
  };

  /// Specification for a more difficult picture.
const int difficult[] =
  { 15, 15,
    // Column constraints
    1, 3,
    1, 2,
    1, 2,
    1, 1,
    1, 2,
    1, 3,
    1, 2,
    1, 4,
    1, 3,
    1, 4,
    2, 2, 1,
    2, 1, 1,
    2, 1, 1,
    2, 1, 1,
    1, 3,
    // Row constraints
    1, 3,
    2, 1, 1,
    2, 1, 1,
    2, 1, 1,
    2, 1, 2,
    1, 5,
    1, 1,
    1, 2,
    1, 1,
    1, 1,
    2, 1, 2,
    2, 1, 2,
    2, 2, 1,
    2, 2, 2,
    1, 3
  };

  /// Specification for a non-unique picture
const int non_unique[] =
  { 11, 15,
    // Column constraints
    1, 5,
    3, 1, 2, 4,
    3, 2, 1, 3,
    4, 2, 2, 1, 1,
    4, 1, 1, 1, 1,
    2, 1, 5,
    5, 2, 1, 1, 3, 2,
    5, 2, 1, 1, 1, 1,
    3, 1, 4, 1,
    2, 1, 1,
    1, 1,
    // Row constraints
    2, 2, 2,
    2, 2, 2,
    1, 4,
    2, 1, 1,
    2, 1, 1,
    4, 1, 1, 1, 1,
    2, 1, 1,
    2, 1, 4,
    3, 1, 1, 1,
    3, 1, 1, 4,
    2, 1, 3,
    2, 1, 2,
    1, 5,
    2, 2, 2,
    2, 3, 3
  };

  /** \brief Specification for a dragonfly-picture.
   *
   * From http://www.oberlin.edu/math/faculty/bosch/pbn-page.html, where
   * it is claimed that it is hard.
   */
  const int dragonfly[] =
    { 20, 20,
      // Column constraints.
      4, 1, 1, 1, 2,
      5, 3, 1, 2, 1, 1,
      5, 1, 4, 2, 1, 1,
      4, 1, 3, 2, 4,
      4, 1, 4, 6, 1,
      3, 1, 11, 1,
      4, 5, 1, 6, 2,
      1, 14,
      2, 7, 2,
      2, 7, 2,
      3, 6, 1, 1,
      2, 9, 2,
      4, 3, 1, 1, 1,
      3, 3, 1, 3,
      3, 2, 1, 3,
      3, 2, 1, 5,
      3, 3, 2, 2,
      3, 3, 3, 2,
      3, 2, 3, 2,
      2, 2, 6,
      // Row constraints
      2, 7, 1,
      3, 1, 1, 2,
      3, 2, 1, 2,
      3, 1, 2, 2,
      3, 4, 2, 3,
      3, 3, 1, 4,
      3, 3, 1, 3,
      3, 2, 1, 4,
      2, 2, 9,
      3, 2, 1, 5,
      2, 2, 7,
      1, 14,
      2, 8, 2,
      3, 6, 2, 2,
      4, 2, 8, 1, 3,
      4, 1, 5, 5, 2,
      5, 1, 3, 2, 4, 1,
      5, 3, 1, 2, 4, 1,
      5, 1, 1, 3, 1, 3,
      4, 2, 1, 1, 2
    };

  /** \brief Specification for a picture of a flower.
   *
   * From http://www.icparc.ic.ac.uk/eclipse/examples/nono.ecl.txt, the
   * hardest instance.
   */
  const int p200[] =
    { 25, 25,
      // Column constraints
      4, 1,1,2,2,
      3, 5,5,7,
      4, 5,2,2,9,
      4, 3,2,3,9,
      5, 1,1,3,2,7,
      3, 3,1,5,
      5, 7,1,1,1,3,
      6, 1,2,1,1,2,1,
      3, 4,2,4,
      4, 1,2,2,2,
      3, 4,6,2,
      4, 1,2,2,1,
      4, 3,3,2,1,
      3, 4,1,15,
      6, 1,1,1,3,1,1,
      6, 2,1,1,2,2,3,
      4, 1,4,4,1,
      4, 1,4,3,2,
      4, 1,1,2,2,
      5, 7,2,3,1,1,
      5, 2,1,1,1,5,
      3, 1,2,5,
      4, 1,1,1,3,
      3, 4,2,1,
      1, 3,
      // Row constraints
      3, 2,2,3,
      5, 4,1,1,1,4,
      5, 4,1,2,1,1,
      7, 4,1,1,1,1,1,1,
      6, 2,1,1,2,3,5,
      6, 1,1,1,1,2,1,
      5, 3,1,5,1,2,
      6, 3,2,2,1,2,2,
      7, 2,1,4,1,1,1,1,
      6, 2,2,1,2,1,2,
      6, 1,1,1,3,2,3,
      5, 1,1,2,7,3,
      5, 1,2,2,1,5,
      5, 3,2,2,1,2,
      4, 3,2,1,2,
      3, 5,1,2,
      4, 2,2,1,2,
      4, 4,2,1,2,
      4, 6,2,3,2,
      4, 7,4,3,2,
      3, 7,4,4,
      3, 7,1,4,
      3, 6,1,4,
      3, 4,2,2,
      2, 2,1
    };

  const int *specs[] = {heart, bear, crocodile, unknown,
                        pinwheel, difficult, non_unique, dragonfly, p200};
  const unsigned n_examples = sizeof(specs)/sizeof(int*);
  //@}

}

// STATISTICS: example-any
