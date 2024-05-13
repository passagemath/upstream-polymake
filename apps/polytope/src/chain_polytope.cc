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
#include "polymake/linalg.h"
#include "polymake/Matrix.h"
#include "polymake/SparseMatrix.h"
#include "polymake/Rational.h"
#include "polymake/graph/Lattice.h"
#include "polymake/graph/Decoration.h"
#include "polymake/Graph.h"
#include "polymake/PowerSet.h"

namespace polymake { namespace polytope {

namespace {

template <typename Decoration, typename SeqType>
SparseMatrix<Int> constraints(const graph::Lattice<Decoration, SeqType>& HD, const Array<Set<Int>>& max_chains, bool is_extended)
{
  const Int d = HD.graph().nodes();
  const Int m = max_chains.size();
  IncidenceMatrix<NonSymmetric> chains(m, d, entire(max_chains));
  auto Ineq = -1 * same_element_sparse_matrix<Int>(chains);

  // nonnegativity constraints from nontrivial poset elements
  auto mat = (zero_vector<Int>(d) | unit_matrix<Int>(d)) / (ones_vector<Int>() | Ineq);

  if (is_extended)
    return remove_zero_rows(mat.minor(All, ~Set<Int>{HD.top_node()+1, HD.bottom_node()+1}));
  else
    return mat;
}

template <typename Decoration, typename SeqType>
SparseMatrix<Int> points(const graph::Lattice<Decoration, SeqType>& HD, const Array<Set<Int>>& max_anti_chains, bool is_extended)
{
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
  const Int d = HD.graph().nodes();
  IncidenceMatrix<NonSymmetric> im(m, d, entire(all_anti_chains));
  auto mat = ones_vector<Int>(m) | same_element_sparse_matrix<Int>(im);
  if (is_extended)
    return mat.minor(All, ~Set<Int>{HD.top_node()+1, HD.bottom_node()+1});
  else
    return mat;
}

}

template <typename Decoration, typename SeqType>
BigObject chain_polytope(BigObject L, bool is_extended=1)
{
  const graph::Lattice<Decoration, SeqType> HD(L);

  const Int d = HD.graph().nodes();
  const Int top = HD.top_node();
  const Int bottom = HD.bottom_node();
  Set<Int> tb1, tb2;
  tb1 += 0; tb1 += d-1; tb2 += top; tb2 += bottom;
  if (tb1 != tb2)
    throw std::runtime_error("non-standard indices for top and bottom");

  const Int dim = d + 1 - 2*is_extended;
  const Array<Set<Int>> max_chains = L.give("MAXIMAL_CHAINS");
  const SparseMatrix<Int> facets = constraints(HD, max_chains, is_extended);
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
                          "# Chain polytope of a poset."
                          "# See Stanley, Discr Comput Geom 1 (1986)."
                          "# @param PartiallyOrderedSet L"
                          "# @param Bool is_extended interpret input as extended poset and ignore top and bottom node"
                          "# @return Polytope<Rational>",
                          "chain_polytope<Decoration,SeqType>(Lattice<Decoration,SeqType>; $=1)");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
