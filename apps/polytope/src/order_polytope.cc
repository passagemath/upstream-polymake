/* Copyright (c) 1997-2024
   Ewgenij Gawrilow, Michael Joswig, and the polymake team
   Technische Universität Berlin, Germany
   https://polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/Decoration.h"
#include "polymake/graph/graph_iterators.h"
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"

namespace polymake { namespace polytope {

namespace {

template <typename Decoration, typename SeqType>
SparseMatrix<Int> constraints(const graph::Lattice<Decoration, SeqType>& HD, bool is_extended)
{
  const Int d = HD.graph().nodes();
  // facets correspond to covering relations, including those incident with top/bottom
  const Int m = HD.graph().edges(); 
  
  const Int top = HD.top_node();
  const Int bottom = HD.bottom_node();
  SparseMatrix<Int> Ineq(m+2,d+1);
  Ineq(0, 0) = 0; Ineq(0, bottom + 1) = 1;
  Ineq(1, 0) = 1; Ineq(1, top + 1) = -1;
  Int i = 2;

  // all other edges
  for (auto e = entire(edges(HD.graph())); !e.at_end(); ++e, ++i) {
     Ineq(i,e.from_node()+1) = -1;
     Ineq(i,e.to_node()+1) = 1;
  }
  
  if (is_extended) {
    // change <= top to <= 1 
    // all other rows will have 0 in the first entry anyway
    Ineq.col(0) = Ineq.col(HD.top_node()+1);
    // >= bottom will turn into >= 0
    // skip first two rows and remove extra cols
    return Ineq.minor(sequence(2,m), ~Set<Int>{HD.top_node()+1, HD.bottom_node()+1});
  } else
    return Ineq;
}
  
template <typename Decoration, typename SeqType>
Matrix<Int> points(const graph::Lattice<Decoration, SeqType>& HD, const Array<Set<Int>>& max_anti_chains, bool is_extended)
{
  const Int d = HD.graph().nodes(); // don't count the top and bottom elements

  // the setup is chosen such that the vertices of order and chain polytopes match
  Set<Set<Int>> all_anti_chains;
  for (auto mac=entire(max_anti_chains); !mac.at_end(); ++mac) {
    all_anti_chains += all_subsets(*mac);
  }
  if (is_extended) {
    // this makes sure we still get vertices when projecting the coordinates later
    all_anti_chains -= scalar2set(HD.top_node());
    all_anti_chains -= scalar2set(HD.bottom_node());
  }

  const Int m = all_anti_chains.size();
  Matrix<Int> Pts(m,d+1);
  Pts.col(0) = ones_vector<Int>(m); // homogenizing coord

  Int i=0;
  for (auto ac=entire(all_anti_chains); !ac.at_end(); ++ac, ++i) {
    auto it=entire(*ac);

    // here the empty anti-chain corresponds to the empty filter
    if (it.at_end()) {
       continue; // nothing to do; matrix initialized as zero
    }

    // nonempty anti-chain
    graph::BFSiterator< Graph<Directed> > j(HD.graph(), *it);
    while (true) {
      while (!j.at_end()) {
        const Int node = *j;
        Pts(i,node+1) = 1;
        ++j;
      }
      ++it;
      if (!it.at_end())
        j.process(*it);
      else
        break;
    }
  }
  if (is_extended)
    return Pts.minor(All, ~Set<Int>{HD.top_node()+1, HD.bottom_node()+1});
  else
    return Pts;
}

}
      
template <typename Decoration, typename SeqType>
BigObject order_polytope(BigObject L, bool is_extended=1)
{
  const graph::Lattice<Decoration,SeqType> HD(L);

  const Int d = HD.graph().nodes();
  const Int top = HD.top_node();
  const Int bottom = HD.bottom_node();
  Set<Int> tb1, tb2;
  tb1 += 0; tb1 += d-1; tb2 += top; tb2 += bottom;
  if (tb1 != tb2)
    throw std::runtime_error("non-standard indices for top and bottom");

  const Int dim = d + 1 - 2*is_extended;
  const SparseMatrix<Int> facets = constraints(HD, is_extended);
  const Array<Set<Int>> max_anti_chains = L.give("MAXIMAL_ANTI_CHAINS");
  const Matrix<Int> vertices = points(HD, max_anti_chains, is_extended);
  const Matrix<Rational> affine_hull(0,dim);
  
  BigObject poly("Polytope<Rational>",
                 "FACETS", facets,
                 "AFFINE_HULL", affine_hull,
                 "VERTICES", vertices,
                 "CONE_DIM", dim,
                 "CONE_AMBIENT_DIM", dim);
  return poly;
}
      
UserFunctionTemplate4perl("#@category Producing a polytope from graphs"
                          "# Order polytope of a poset."
                          "# See Stanley, Discr Comput Geom 1 (1986)."
                          "# @param PartiallyOrderedSet L"
                          "# @param Bool is_extended interpret input as extended poset and ignore top and bottom node"
                          "# @return Polytope<Rational>",
                          "order_polytope<Decoration, SeqType>(Lattice<Decoration,SeqType>; $=1)");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
