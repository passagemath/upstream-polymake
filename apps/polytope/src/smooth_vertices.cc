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

#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Graph.h"
#include "polymake/linalg.h"
#include "polymake/common/lattice_tools.h"

namespace polymake { namespace polytope {

Set<Int> smooth_vertices(BigObject P)
{
   const Int d = static_cast<Int>(P.give("CONE_DIM")) - 1;
   const Int ad = P.give("CONE_AMBIENT_DIM");
   const Graph<>& G = P.give("GRAPH.ADJACENCY");
   const Set<Int>& far_face = P.give("FAR_FACE");
   const Int n = G.nodes();
   auto actual_vertices = sequence(0,n) - far_face;
   const EdgeMap<Undirected,Vector<Rational>>& edge_directions = P.give("GRAPH.EDGE_DIRECTIONS");
   Set<Int> sv;

   for (const auto& v : actual_vertices) {
      ListMatrix<Vector<Integer>> M;
      for (const auto& e : G.out_edges(v)) {
         M /= common::primitive(edge_directions[e]);
      }
      if (M.rows() != d) continue; // not a simple vertex

      Integer g = 0;
      for (auto ss : all_subsets_of_k(sequence(1, ad - 1), d)) {
         g = gcd(g, abs(det(M.minor(All,ss))));
         if (g == 1)
            break;
      }
      if (g == 1)
         sv += v;
   }
   
   return sv;
}

UserFunction4perl("# @category Lattice points in cones"
                  "# Return the indices of all bounded and simple vertices where the"
                  "# edge-directions form a lattice basis."
                  "# @example"
                  "# > print smooth_vertices(cube(2, 2/3));"
                  "# | {0 1 2 3}",
                  &smooth_vertices, "smooth_vertices(Polytope)");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
