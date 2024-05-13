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

namespace nsw_sphere {

   // d is the number of paths,
   // 2d-1 is the dimension of the sphere, so that
   // 2d is the number of vertices in a simplex
   
using IndexTuple = Vector<Int>;  // size d
using LexIndexCollection = Set<IndexTuple>;

using IJIndex = std::pair<Int,Int>;
using IJIndices = Array<IJIndex>; // size d
using SimplexSet = Set<Int>; // size 2d

using IndexOfLabel = Map<std::pair<Int,Int>, Int>;
using IJLabels = Array<std::pair<Int,Int>>;
      
using BoundaryLexIndexCollection = Set<IndexTuple>;
   
inline      
Int
index_of_a_i_j(const Int n,
               const Int i,
               const Int j)
{
   //   return (j-1)*n + i-1;
   return j*n + i;
}

inline
Int
index_of_new_vertex(const Int n,
                    const Int d,
                    const Int k)
{
   return index_of_a_i_j(n, n-1, d-1) + 1 + k;
}

class Simplex {      
public:
   IJIndices ij_indices;
   Int index_sum=0;
   SimplexSet vertices;

   Simplex() {}
   
   Simplex(const Int n,
           const IndexTuple& it)
      : ij_indices(it.size())
      , index_sum(0)
   {
      for (Int j=0; j < it.size(); ++j) {
         vertices += index_of_a_i_j(n, it[j]  , j);
         vertices += index_of_a_i_j(n, it[j]+1, j);
         ij_indices[j] = std::make_pair(it[j], j);
         index_sum += 1 + it[j];
      }
   }

   std::vector<Int> i_s() const {
      std::vector<Int> is;
      is.reserve(ij_indices.size());
      for (auto it = entire(ij_indices); !it.at_end(); ++it)
         is.push_back(it->first);
      return is;
   }
   
   struct Compare {
      ::pm::cmp_value operator()(const Simplex& lhs,
                                 const Simplex& rhs) const {
         return operations::cmp()(lhs.vertices, rhs.vertices);
      }
   };
   
   template<typename Output>
   friend
   Output& operator<< (GenericOutput<Output>& outs, const Simplex& s) {
      return outs.top() << "[" << s.vertices
                        << "=" << s.ij_indices
                        << " sum " << s.index_sum << "]";
   }
};

Set<Int>
removed_ridge(const Simplex& sigma,
              const Int n,
              const Int ell,
              const Int m);
   
class ShellingOrderedRidge {
public:
   IJIndices ij_indices;
   std::pair<Int,Int> ell_m;
   SimplexSet vertices;

   ShellingOrderedRidge() {}

   ShellingOrderedRidge(const Simplex& s,
                        const Int n,
                        const Int ell,
                        const Int m)
      : ij_indices(s.ij_indices)
      , ell_m(std::make_pair(ell,m))
      , vertices(removed_ridge(s, n, ell, m))
   {}

   struct Compare {
      ::pm::cmp_value operator()(const ShellingOrderedRidge& lhs,
                                 const ShellingOrderedRidge& rhs) const {
         return operations::cmp()(lhs.ell_m, rhs.ell_m);
      }
   };
   
   template<typename Output>
   friend
   Output& operator<< (GenericOutput<Output>& outs, const ShellingOrderedRidge& s) {
      return outs.top() << "[" << s.vertices
                        << "=" << s.ij_indices
                        << " l,m=" << s.ell_m << "]";
   }
};

::pm::cmp_value
def_3_4_cmp(const SimplexSet& tau,
            const SimplexSet& taup,
            const IJLabels& ij_labels,
            const Int d);

extern IJLabels const *global_ij_label_ptr;
extern Int global_d;
   
class Def34OrderedSimplexSet {
public:
   const SimplexSet tau;

   Def34OrderedSimplexSet() {}
   Def34OrderedSimplexSet(const SimplexSet& tau_)
      : tau(tau_)
   {}

   struct Compare {
      ::pm::cmp_value operator()(const Def34OrderedSimplexSet& lhs,
                                 const Def34OrderedSimplexSet& rhs) const {
         return def_3_4_cmp(lhs.tau, rhs.tau, *global_ij_label_ptr, global_d);
      }
   };
   
   template<typename Output>
   friend
   Output& operator<< (GenericOutput<Output>& outs, const Def34OrderedSimplexSet& s) {
      return outs.top() << s.tau;
   }
   
};

class ShellingOrderedSubridge38 {
   // the 38 refers to Lemma 3.8
public:
   std::pair<Int,Int> ell_r;
   SimplexSet vertices;

   ShellingOrderedSubridge38() {}

   ShellingOrderedSubridge38(const SimplexSet& rest,
                             const Int ell,
                             const Int r)
      : ell_r(std::make_pair(ell,r))
        // ell is 0-based, r is 1-based
        // so that in case 2 of Def 3.7 r can also be 0
      , vertices(rest)
   {}

   struct Compare {
      ::pm::cmp_value operator()(const ShellingOrderedSubridge38& lhs,
                                 const ShellingOrderedSubridge38& rhs) const {
         return operations::cmp()(lhs.ell_r, rhs.ell_r);
      }
   };
   
   template<typename Output>
   friend
   Output& operator<< (GenericOutput<Output>& outs, const ShellingOrderedSubridge38& s) {
      return outs.top() << "[" << s.vertices
                        << " (l=" << s.ell_r.first
                        << ", r=" << s.ell_r.second << ")]";
   }
};
   
}}}
   
namespace pm {
template<>   
struct is_ordered<polymake::topaz::nsw_sphere::Simplex> : std::true_type {};

template<>   
struct is_ordered<polymake::topaz::nsw_sphere::ShellingOrderedRidge> : std::true_type {};

template<>   
struct is_ordered<polymake::topaz::nsw_sphere::Def34OrderedSimplexSet> : std::true_type {};

template<>   
struct is_ordered<polymake::topaz::nsw_sphere::ShellingOrderedSubridge38> : std::true_type {};
   
namespace operations {

template<>
struct define_comparator<polymake::topaz::nsw_sphere::Simplex, polymake::topaz::nsw_sphere::Simplex, pm::operations::cmp, pm::is_opaque, pm::is_opaque, pm::cons<pm::is_opaque, pm::is_opaque>, pm::cmp_value> {
   typedef polymake::topaz::nsw_sphere::Simplex::Compare type;
};

template<>
struct define_comparator<polymake::topaz::nsw_sphere::ShellingOrderedRidge, polymake::topaz::nsw_sphere::ShellingOrderedRidge, pm::operations::cmp, pm::is_opaque, pm::is_opaque, pm::cons<pm::is_opaque, pm::is_opaque>, pm::cmp_value> {
   typedef polymake::topaz::nsw_sphere::ShellingOrderedRidge::Compare type;
};

template<>
struct define_comparator<polymake::topaz::nsw_sphere::Def34OrderedSimplexSet, polymake::topaz::nsw_sphere::Def34OrderedSimplexSet, pm::operations::cmp, pm::is_opaque, pm::is_opaque, pm::cons<pm::is_opaque, pm::is_opaque>, pm::cmp_value> {
   typedef polymake::topaz::nsw_sphere::Def34OrderedSimplexSet::Compare type;
};

template<>
struct define_comparator<polymake::topaz::nsw_sphere::ShellingOrderedSubridge38, polymake::topaz::nsw_sphere::ShellingOrderedSubridge38, pm::operations::cmp, pm::is_opaque, pm::is_opaque, pm::cons<pm::is_opaque, pm::is_opaque>, pm::cmp_value> {
   typedef polymake::topaz::nsw_sphere::ShellingOrderedSubridge38::Compare type;
};
   
} }
   
namespace polymake { namespace topaz {

namespace nsw_sphere {      
      
class SimplexStorer {

   const Int n;
   const Int d;
   Map<IndexTuple, Simplex> simplex_at;

public:
   SimplexStorer(const Int n_,
                 const Int d_)
      : n(n_) // length of paths
      , d(d_) // number of paths; sphere is of dim D = 2d-1
   {}

   Simplex
   simplex(const IndexTuple& it) {
      if (simplex_at.exists(it))
         return simplex_at[it];

      if (d != it.size()) {
         cerr << "received tuple " << it 
              << ", but should be of size " << d << endl;
         throw std::runtime_error("SimplexStorer::simplex: wrong tuple size");
      }
      for (const Int i: it)
         if (i >= n)
            throw std::runtime_error("illegal index in SimplexStorer::simplex()");

      Simplex s(n, it);
      simplex_at[it] = s;
      return s;
   }
};


using BallFacet = Simplex;
using BoundaryFacet = SimplexSet;
using FacetsOfBall = Set<BallFacet>;
using FacetsOfBoundary = Set<BoundaryFacet>;

SimplexSet
missing_face_of(const SimplexSet& sigma,
                const FacetsOfBoundary& delta_B);

Set<BoundaryFacet>
D_sigma_of(const SimplexSet& sigma,
           const FacetsOfBoundary& delta_B);   
   
struct dDBallData {
   Int n;
   Int d;
   Int k_max;
   SimplexStorer ss;
   Array<FacetsOfBall> B_array;
   Array<FacetsOfBoundary> delta_B_array;
   Array<FacetsOfBall> S_low_array;
   Array<FacetsOfBall> S_up_array;
   Array<FacetsOfBoundary> union_of_Ds_array;
   
   dDBallData(const Int n_,
              const Int d_,
              const Int k_max_)
      : n(n_)
      , d(d_)
      , k_max(k_max_)
      , ss(n, d)
      , B_array(k_max)
      , delta_B_array(k_max)
      , S_low_array(k_max)
      , S_up_array(k_max)
      , union_of_Ds_array(k_max)
   {}

   template<typename Output>
   friend
   Output& operator<< (GenericOutput<Output>& outs, const dDBallData& bd) {
      Output& os = outs.top();
      for (Int k=0; k < bd.k_max; ++k) {
         os << k+1 << ": ";
         for (const auto& s: bd.B_array[k])
            os << s << " ";
         os << endl;

         os << "bd " << k+1 << ": ";
         for (const auto& s: bd.delta_B_array[k])
            os << s << " ";
         os << endl;

         os << "S_low_array " << k+1 << ": " << bd.S_low_array[k]
            << endl
            << "S_up_array " << k+1 << ": " << bd.S_up_array[k]
            << endl;
      }
      return os;
   }
};

enum class TriangulationChoice {
   one, two
};      

Set<SimplexSet>
C_sigma_tilde_of_impl(const SimplexSet& F_sigma,
                      const Int o_k,
                      const SimplexSet& verts_of_D_sigma,
                      const TriangulationChoice& choice);
   
struct ModifiedDiagonals {
   TriangulationChoice choice;
   Set<SimplexSet> S_low_tilde;
   Set<SimplexSet> S_up_tilde;
   Set<SimplexSet> connecting_path;
};

void
check_lemma_2_2(const IndexOfLabel& iol,
                const dDBallData& bd,
                const IJLabels& ij_labels,
                const Int verbosity,
                bool& passed);

bool
check_constructibility_proof(const dDBallData& bd,
			     const Array<ModifiedDiagonals>& md,
                             const IJLabels& ij_labels,
                             const Int verbosity,
                             const bool output_on_error);

// bool
// check_shellability_proof(const dDBallData& bd,
//                          const Array<ModifiedDiagonals>& md,
//                          const IJLabels& ij_labels,
//                          const std::string& check_constructibility,
//                          const Int verbosity);
   
template<typename T>   
void   
add_to_boundary(const Simplex& sigma,
                Set<T>& boundary)
{
   for (auto it = entire(all_subsets_less_1(sigma.vertices)); !it.at_end(); ++it) {
      const T bf(*it);
      if (boundary.contains(bf))
         boundary -= bf;
      else
         boundary += bf;
   }
}

template<typename T>
void   
add_to_boundary(const SimplexSet& sigma,
                Set<T>& boundary)
{
   for (auto it = entire(all_subsets_less_1(sigma)); !it.at_end(); ++it) {
      const T bf(*it);
      if (boundary.contains(bf))
         boundary -= bf;
      else
         boundary += bf;
   }
}

std::string
comma_if_not_first(bool& first,
                   const std::string filler = ",");

struct Label {
   std::string text;

   Label(const Vector<Int>& v,
         const dDBallData& bd,
         std::stringstream& ss)
   {
      ss.str("");
      bool first(true);
      for (auto it = entire<indexed>(v); !it.at_end(); ++it)
         ss << comma_if_not_first(first, " ")
            << (1 + *it) << "^" << (1 + it.index());
      text = ss.str();
   }

   Label(const SimplexSet& sigma,
         const dDBallData& bd,
         std::stringstream& ss)
   {
      ss.str("");
      bool first(true);
      for (const Int i: sigma)
         ss << comma_if_not_first(first, " ")
            << (1 + (i % bd.n)) << "^" << (1 + (i / bd.n));
      text = ss.str();
   }
   
   template<typename Output>
   friend
   Output& operator<< (GenericOutput<Output>& outs, const Label& l) {
      return outs.top() << l.text;
   }
};

} // end namespace nsw_sphere
      
} }


// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
    
