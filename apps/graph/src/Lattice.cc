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
#include "polymake/graph/maximal_chains.h"
#include "polymake/graph/max_cliques.h"

namespace polymake { namespace graph {

template <typename Decoration, typename SeqType>
Array<Set<Int>> lattice_dual_faces(BigObject lattice_obj)
{
  return Lattice<Decoration, SeqType>(lattice_obj).dual_faces();
}

template <typename Decoration, typename SeqType, typename Permutation>
BigObject lattice_permuted_faces(BigObject lattice_obj, const Permutation& perm)
{
  return static_cast<BigObject>((Lattice<Decoration, SeqType>(lattice_obj)).permuted_faces(perm));
}

template <typename Decoration, typename SeqType>
Array<Set<Int>> lattice_maximal_chains(BigObject lattice_obj)
{
  const Lattice<Decoration, SeqType> HD(lattice_obj);
  return maximal_chains(HD, false, false);
}

template <typename Decoration, typename SeqType>
Graph<Undirected> lattice_comparability_graph(BigObject lattice_obj)
{
  const Lattice<Decoration, SeqType> HD(lattice_obj);
  const Int d = HD.graph().nodes();
  const Array<Set<Int>> max_chains = lattice_obj.give("MAXIMAL_CHAINS");
  Graph<Undirected> CG(d);
  for (auto c = entire(max_chains); !c.at_end(); ++c) {
    if (c->size()>1) {
      for (auto pair = entire(all_subsets_of_k(*c,2)); !pair.at_end(); ++pair) {
        CG.edge(pair->front(), pair->back());
      }
    }
  }
  return CG;
}

FunctionTemplate4perl("lattice_dual_faces<Decoration, SeqType>(Lattice<Decoration, SeqType>)");
FunctionTemplate4perl("lattice_permuted_faces<Decoration, SeqType, Permutation>(Lattice<Decoration,SeqType>, Permutation)");
FunctionTemplate4perl("lattice_maximal_chains<Decoration, SeqType>(Lattice<Decoration, SeqType>)");
FunctionTemplate4perl("lattice_comparability_graph<Decoration, SeqType>(Lattice<Decoration, SeqType>)");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
