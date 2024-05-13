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

namespace polymake { namespace graph {

BigObject maximal_ranked_poset(const Array<Int>& tau) {
   const Int k = tau.size();
   const Int n = accumulate(tau, operations::add()); // sum of the entries

   Graph<Directed> g(n+2);
   NodeMap<Directed, lattice::BasicDecoration> nm(g);

   // process bottom element
   const Int bottom = 0;
   nm[bottom] = lattice::BasicDecoration(scalar2set(bottom), 0);
   for (Int v=1; v<=tau[0]; ++v) { g.add_edge(0,v); }
   // process middle part
   Int t = 1;
   for (Int i=1; i<k; ++i) {
      for (Int a=0; a<tau[i-1]; ++a) {
         const Int u = t+a;
         nm[u] = lattice::BasicDecoration(scalar2set(u), i);
         for (Int b=0; b<tau[i]; ++b) {
            const Int v = t+tau[i-1]+b;
            g.add_edge(u,v);
         }
      } 
      t += tau[i-1];
   }
   // process top element
   const Int top = n+1;
   nm[top] = lattice::BasicDecoration(scalar2set(top), k+1);
   for (Int a=0; a<tau[k-1]; ++a) {
      const Int u = t+a;
      nm[u] = lattice::BasicDecoration(scalar2set(u), k);
      g.add_edge(u, top);
   }

  return BigObject("Lattice<BasicDecoration,Sequential>",
                   "ADJACENCY", g,
                   "DECORATION", nm,
                   "BOTTOM_NODE", bottom,
                   "TOP_NODE", top);
}
      
UserFunction4perl("#@category Producing a graph"
                  "# Maximal ranked partially ordered set."
                  "# See Ahmad, Fourier, Joswig, arXiv:2309.01626"
                  "# @param Array<Int> tau"
                  "# @return PartiallyOrderedSet<BasicDecoration,Sequential>",
                  &maximal_ranked_poset, "maximal_ranked_poset(Array<Int>)");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
