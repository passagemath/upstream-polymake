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

#include "polymake/client.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/topaz/nsw_3_spheres.h"

namespace polymake { namespace topaz {

namespace {

AIndex a_index;
BIndex b_index;

} // end anonymous namespace

namespace nsw_sphere {
   
SimplexSet
missing_face_of(const SimplexSet& sigma,
                const FacetsOfBoundary& delta_B)
{
   SimplexSet intersection;
   bool started_intersecting(false);
   for (const auto bf: all_subsets_less_1(sigma)) {
      if (delta_B.contains(bf))
         continue;

      if (!started_intersecting) {
         started_intersecting = true;
         intersection = bf;
      } else
         intersection *= bf;
   }
   return intersection;
}
   
Set<BoundaryFacet>
D_sigma_of(const SimplexSet& sigma,
           const FacetsOfBoundary& delta_B)
{
   Set<BoundaryFacet> D_sigma;
   for (const auto bf: all_subsets_less_1(sigma))
      if (delta_B.contains(bf)) 
         D_sigma += bf;

   return D_sigma;
}

bool
is_ball_data_compatible(const BallData& bd)
{
   for (Int k=0; k < bd.k_max; ++k) {
      Set<Set<Int>> missing_faces;
      for (const auto& tet: bd.S[k]) {
         const Set<Int> mf(nsw_sphere::missing_face_of(tet, bd.delta_B[k]));
         for (const auto& boundary_tet: bd.delta_B[k])
            if (incl(mf, boundary_tet) <= 0) {
               cerr << "missing face found in boundary" << endl;
               return false;
            }
         missing_faces += mf;
      }
      if (missing_faces.size() != bd.S[k].size()) {
         cerr << "repeated missing face" << endl;
         return false;
      }
   }
   return true;
}

} namespace {

Set<RidgeFacet>
boundary_of_D_sigma(const Set<BoundaryFacet>& D_sigma)
{
   Set<RidgeFacet> bd;
   for (const auto& sigma: D_sigma)
      for (const auto rf: all_subsets_less_1(sigma))
         if (bd.contains(rf)) // RidgeFacet
            bd -= rf;
         else
            bd += rf;
   
   return bd;
}

Set<BoundaryFacet>   
boundary_of_cell_complex(const CellComplex& C)
{
   Set<BoundaryFacet> bd;
   for (const BoundaryOfCell& boc: C) {
      for (const BoundaryFacet& bf: boc) {
         if (bd.contains(bf))
            bd -= bf;
         else
            bd += bf;
      }
   }
   return bd;
}
   
CellComplex
S_filling(const FacetsOfBall& B,
          const FacetsOfBoundary& delta_B,
          const CompatibleFamily& S,
          const Int v)
{
   CellComplex modified_cells;
   Set<BoundaryFacet> all_facets_of_all_D_sigmas;

   for (const auto& sigma: S) {
      const Set<BoundaryFacet> D_sigma = nsw_sphere::D_sigma_of(sigma, delta_B);
      const Set<RidgeFacet> bd_D_sigma = boundary_of_D_sigma(D_sigma);
      
      // For each σ ∈ S,
      // consider the cell C_σ with boundary complex D_σ ∪ (∂ D_σ ∗ v).
      BoundaryOfCell boc;

      // we store polyhedral cells as the set of their boundary facets
      for (const auto& tau: D_sigma) {
         all_facets_of_all_D_sigmas += tau;
         boc += tau;
      }

      for (const auto& ridge: bd_D_sigma)
         boc += ridge + scalar2set(v);

      modified_cells += boc;
   }

   // For each (d − 1)-simplex τ of ∂ B that is not in any D_σ , σ ∈ S,
   // consider the simplex τ ∗ v.
   for (const auto& tau: delta_B) {
      if (all_facets_of_all_D_sigmas.contains(tau))
         continue;

      BoundaryOfCell boc { tau }; // the boundary of the cell τ ∗ v
      for (auto it = entire(all_subsets_less_1(tau)); !it.at_end(); ++it) 
         boc += *it + scalar2set(v);

      modified_cells += boc;
   }

   return modified_cells;
}

bool
check_intersection(const CellComplex& cells)
{
   if (cells.size() < 2)
      return true;
   
   for (auto it = entire(all_subsets_of_k(cells, 2)); !it.at_end(); ++it) {
      const auto intersect = it->front() * it->back();
      if (intersect.size() >= 2) {
         cerr << "bad intersection between cells with boundaries "
              << it->front() << " and "
              << it->back() << endl;
         return false;
      }
   }
   return true;
}
   

BallData
ball_data_for_4d_join_of_paths(const Int n,
                               const Int k_max,
                               const Int verbosity)
{
   BallData bd(n, k_max);

   for (Int i=1; i < n; ++i) {
      for (Int j=1; j < n; ++j) {
         const SimplexSet tet = bd.ss.tetrahedron(i,j);

         // determine k in { 1, ..., k_max } such that i+j in [4k-3,...,4k]
         const Int
            div4 = (i+j) / 4,
            mod4 = (i+j) % 4,
            k ( 0 == mod4
                ? div4
                : div4 + 1 );

         if (k > k_max) {
            // the tet doesn't go into any B,
            // so we store its boundary complex for later
            BoundaryOfCell bd_tet;
            for (auto it = entire(all_subsets_less_1(tet)); !it.at_end(); ++it)
               bd_tet += *it;

            bd.not_in_Bs += bd_tet;
            continue;
         }

         bd.B[k-1] += tet;
         for (auto it = entire(all_subsets_less_1(tet)); !it.at_end(); ++it) {
            const BoundaryFacet bf(*it);
            if (bd.delta_B[k-1].contains(bf))
               bd.delta_B[k-1] -= bf;
            else
               bd.delta_B[k-1] += bf;
         }
      }
   }

   for (Int i=1; i < n; ++i) {
      for (Int j=1; j < n; ++j) {
         const Int mod4 = (i+j) % 4;

         // make S_k if i+j = 0 or 3 mod 4
         if ( 1 == mod4 || 2 == mod4 )
            continue;
         
         const Int
            div4 = (i+j) / 4,
            k ( 0 == mod4
                ? div4
                : div4 + 1 );
         if (k > k_max)
            continue;

         const SimplexSet tet = bd.ss.tetrahedron(i,j);
         Int n_triangles_in_boundary(0);
         for (auto it = entire(all_subsets_less_1(tet)); !it.at_end(); ++it)
            if (bd.delta_B[k-1].contains(*it))
               ++n_triangles_in_boundary;

         if (n_triangles_in_boundary > 1)
            bd.S[k-1] += tet;
      }
   }

   if (verbosity)
      cerr << "ball data:\n" << bd << endl;

   if (verbosity)
      cerr << "checking compatibility of S_k... ";

   if (!nsw_sphere::is_ball_data_compatible(bd))
      throw std::runtime_error("ball data S_k should be compatible");
   if (verbosity)
      cerr << "ok" << endl;

   return bd;
}



CellComplex
modify_balls(const BallData& bd,
             Int& next_v_index,
             const Int verbosity)
{
   // We make the modified complex by starting out with all non-modified tets
   CellComplex K_prime(bd.not_in_Bs);

   if (verbosity)
      cerr << "checking intersections of not_in_Bs... ";
   if (!check_intersection(K_prime))
      throw std::runtime_error("cell complex not_in_Bs has bad intersections");
   if (verbosity)
      cerr << "ok" << endl;
   
   // and modify the B_k to B'_k
   next_v_index = b_index(bd.n) + 1;
   for (Int k=0; k < bd.k_max; ++k) {
      const auto Sf = S_filling(bd.B[k], bd.delta_B[k], bd.S[k], next_v_index++);
      if (verbosity) {
         cerr << "S_filling " << k << ":\n";
         for (const auto& cell_bd: Sf) {
            Set<Int> cell;
            for (const auto& cbd: cell_bd)
               cell += cbd;
            cerr << cell << " = " << cell_bd
                 << endl;
         }
      }
      K_prime += Sf;
   }

   if (verbosity)
      cerr << "checking intersections of K'... ";
   if (!check_intersection(K_prime))
      throw std::runtime_error("cell complex K' has bad intersections");
   if (verbosity)
      cerr << "ok" << endl;

   if (verbosity) {
      const Set<BoundaryFacet> bd_Kp(boundary_of_cell_complex(K_prime));
      cerr << "boundary of K_prime:\n" << bd_Kp << endl;
      Map<Set<Int>, Set<Set<Int>>> in_bd;
      for (Int k=0; k < bd.k_max; ++k) 
         for (const BallFacet& bf: bd.B[k]) 
            for (auto it = entire(all_subsets_less_1(bf)); !it.at_end(); ++it) 
               if (bd_Kp.contains(*it))
                  in_bd[bf] += *it;
      cerr << "bd of B in_bd:\n";
      for (const auto& kv: in_bd)
         cerr << kv.first << ": " << kv.second << endl;

      in_bd.clear();
      for (Int k=0; k < bd.k_max; ++k) 
         for (const BallFacet& bf: bd.S[k]) 
            for (auto it = entire(all_subsets_less_1(bf)); !it.at_end(); ++it) 
               if (bd_Kp.contains(*it))
                  in_bd[bf] += *it;
      cerr << "\nbd of S in_bd:\n";
      for (const auto& kv: in_bd)
         cerr << kv.first << ": " << kv.second << endl;
      cerr << endl;
   }

   return K_prime;
}

} // end anonymous namespace


BigObject
bipyramidal_3_sphere(const Int n,
                     OptionSet options)
{
   if (n<3) throw std::runtime_error("need n>=3");
   const Int verbosity = options["verbosity"];
   
   // In order to avoid as many off-by-one errors as possible,
   // we will use the same 1-based indexing as the paper,
   // and let the AIndex and BIndex classes take care of converting
   // indices to consecutive-zero-based
   
   b_index.set_n(n);
   
   const Int k_max( n % 2
                    ? n/2+1
                    : n/2 );

   if (verbosity)
      cerr << "n = " << n << ", k_max = " << k_max << endl;
   
   // Construction 1 starts with the balls B_k and the compatible systems S_k
   const BallData bd(ball_data_for_4d_join_of_paths(n, k_max, verbosity));
   Int next_v_index(0);
   const CellComplex K_prime(modify_balls(bd, next_v_index, verbosity));
   
   IncidenceMatrix<> facets(K_prime.size(), next_v_index);
   auto rit = entire(rows(facets));
   for (const BoundaryOfCell& bd_Kp: K_prime) {
      Set<Int> cell;
      for (const BoundaryFacet& tau: bd_Kp)
         cell += tau;
      *rit = cell;
      ++rit;
   }

   Array<std::string> labels(next_v_index);
   for (Int i=1; i<=n; ++i)
      labels[a_index(i)] = "A" + std::to_string(i);

   for (Int i=1; i<=n; ++i)
      labels[b_index(i)] = "B" + std::to_string(i);

   for (Int i = b_index(n) + 1; i < next_v_index; ++i)
      labels[i] = "V" + std::to_string(i - b_index(n));
   
   BigObject s("fan::PolyhedralComplex<Rational>");
     s.set_description() << "Nevo-Santos-Wilson 3-sphere with bipyramids and tetrahedra as facets from paths with " << n << " vertices" << endl;
   s.take("INPUT_POLYTOPES") << facets;
   s.take("VERTEX_LABELS") << labels;
   
   return s;
}

      
InsertEmbeddedRule("REQUIRE_APPLICATION fan\n\n");

UserFunction4perl("# @category Producing from scratch"
                  "# Create the 3-sphere with bipyramidal and tetrahedral facets from [Nevo, Santos & Wilson, Many triangulated odd-dimensional spheres, Math Ann 364 (2016), 737-762"
                  "# @param Int n an integer >= 3"
                  "# @option Int verbosity default 0"
                  "# @return fan::PolyhedralComplex<Rational>",
                  &bipyramidal_3_sphere, "bipyramidal_3_sphere($ { verbosity => 0 })");


} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
    
    
