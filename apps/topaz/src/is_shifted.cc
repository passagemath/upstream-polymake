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
#include "polymake/graph/Lattice.h"
#include "polymake/hash_set"

namespace polymake { namespace topaz {

using graph::Lattice;

template <typename Decoration, typename SeqType>
bool is_shifted(const BigObject HD_obj)
{
   const Lattice<Decoration, SeqType> HD(HD_obj);
   const Int d = HD.rank()-1;

   // For a simplicial complex we require that the vertices are consecutively ordered.
   // So they always form a shifted 1-collection and thus do not need to be tested.
   for (Int k = d; k > 0; --k) {
      hash_set<Set<Int>> faces(HD.nodes_of_rank(k).size());
      for (const auto n : HD.nodes_of_rank(k))
         faces.insert(HD.face(n));
      for (const auto n : HD.nodes_of_rank(k)) {
         const Set<Int>& f = HD.face(n);
         for (const auto i : f) {
            Set<Int> g = f; g -= i;
            for (Int j=0; j<i; ++j) {
               const Set<Int> shifted_f = g+j;
               if (shifted_f.size()==k && !faces.contains(shifted_f))
                  return false;
            }
         }
      }
   }

   return true;
}

FunctionTemplate4perl("is_shifted<Decoration, SeqType>(Lattice<Decoration, SeqType>)");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
