/* Copyright (c) 1997-2023
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

#pragma once

#include "polymake/Set.h"
#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Map.h"

namespace polymake { namespace topaz {

using IndexPair = std::pair<Int,Int>;
using IndexPairTuple = std::vector<IndexPair>;

using SimplexSet = Set<Int>; // size 2d
using BallFacet = SimplexSet;
using FacetsOfBall = Set<BallFacet>;
using CompatibleFamily = FacetsOfBall;
   
using Triangle = Set<Int>;
using BoundaryFacet = Triangle;
using FacetsOfBoundary = Set<BoundaryFacet>;

using RidgeFacet = Set<Int>;
using BoundaryOfCell = Set<BoundaryFacet>;
using CellComplex = Set<BoundaryOfCell>;

      
class AIndex {
public:
   AIndex() {}
   Int operator()(const Int i) { return i-1; }
};

class BIndex {
   Int n;

public:
   BIndex() {}
   void set_n(const Int n_) { n = n_; }
   Int operator()(const Int i) { return n + i - 1; }
};

namespace {      
   extern AIndex a_index;
   extern BIndex b_index;
}

      
class TetrahedronStorer {

   const Int n;
   Map<IndexPair, SimplexSet> tetrahedron_at;

public:
   TetrahedronStorer(const Int n_)
      : n(n_)
   {}

   SimplexSet
   tetrahedron(const Int i, const Int j) {
      if (i >= n || j >= n)
         throw std::runtime_error("illegal index in TetrahedronFacets::tetrahedron()");
      
      const IndexPair ip(i,j);
      if (tetrahedron_at.exists(ip))
         return tetrahedron_at[ip];
      SimplexSet tet{ a_index(i), a_index(i+1), b_index(j), b_index(j+1) };
                             // { a_index(i), a_index(i+1), b_index(j) } ,
                             // { a_index(i), a_index(i+1), b_index(j+1) } ,
                             // { a_index(i), b_index(j), b_index(j+1) } ,
                             // { a_index(i+1), b_index(j), b_index(j+1) } };
      // cerr << "tet(" << i << ", " << j << ") = " << tet << endl;
      tetrahedron_at[ip] = tet;
      return tet;
   }
};

      
      
namespace nsw_sphere {

SimplexSet
missing_face_of(const SimplexSet& sigma,
                const FacetsOfBoundary& delta_B);

Set<BoundaryFacet>
D_sigma_of(const SimplexSet& sigma,
           const FacetsOfBoundary& delta_B);
   
} // end namespace nsw_sphere

      
struct BallData {
   Int n;
   Int k_max;
   TetrahedronStorer ss;
   Array<FacetsOfBall> B;
   Array<FacetsOfBoundary> delta_B;
   Array<FacetsOfBall> S;
   CellComplex not_in_Bs;
   
   BallData(const Int n_,
            const Int k_max_)
      : n(n_)
      , k_max(k_max_)
      , ss(n)
      , B(k_max)
      , delta_B(k_max)
      , S(k_max)
   {}

   template<typename Output>
   friend
   Output& operator<< (GenericOutput<Output>& outs, const BallData& bd) {
      Output& os = outs.top();
      for (Int k=0; k < bd.k_max; ++k) {
         os << k << ": ";
         for (const auto& s: bd.B[k])
            os << s << " ";
         os << endl;

         os << "bd " << k << ": ";
         for (const auto& s: bd.delta_B[k])
            os << s << " ";
         os << endl;

         os << "S " << k << ": ";
         for (const auto& s: bd.S[k])
            os << s << " (mf: "
               << nsw_sphere::missing_face_of(s, bd.delta_B[k]) << ") ";
         os << endl;         
      }
      return os;
   }
};

      
namespace nsw_sphere {

bool
is_ball_data_compatible(const BallData& bd);

}

} }


// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
