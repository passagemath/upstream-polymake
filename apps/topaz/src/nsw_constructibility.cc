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

#include "polymake/FacetList.h"
#include "polymake/PowerSet.h"
#include "polymake/MultiDimCounter.h"
#include "polymake/ListMatrix.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/topaz/nsw_d_spheres.h"

namespace polymake { namespace topaz {

namespace nsw_sphere {

IJLabels const *global_ij_label_ptr(nullptr);
Int global_d;

  /************************
    Lemma 2.2
  *************************/
  
void
check_lemma_2_2(const IndexOfLabel& iol,
                const dDBallData& bd,
		const IJLabels& ij_labels,
                const Int verbosity,
                bool& passed)
{
   if (verbosity)
      cerr << "Checking Lemma 2.2...";
   FacetsOfBoundary boundary;
   for (Int k=0; k < bd.k_max; ++k) 
      for (const Simplex& sigma: bd.B_array[k])
	add_to_boundary(sigma, boundary);
   
   FacetsOfBoundary supposed_boundary;
   for (Int k=0; k < bd.k_max; ++k) {   
     if (verbosity > 1)
       cerr << ".";
      for (const Simplex& sigma: bd.B_array[k]) {
         for (Int ell=0; ell < bd.d; ++ell) {
            const Int i0(iol[std::make_pair(0, ell)]);
            const Int i1(iol[std::make_pair(1, ell)]);
            if (sigma.vertices.contains(i0)) 
               supposed_boundary += sigma.vertices - i1;

            const Int im2(iol[std::make_pair(bd.n-2, ell)]);
            const Int im1(iol[std::make_pair(bd.n-1, ell)]);
            if (sigma.vertices.contains(im1)) 
               supposed_boundary += sigma.vertices - im2;
         }
      }
   }

   const FacetsOfBoundary
     bms(boundary - supposed_boundary),
     smb(supposed_boundary - boundary);
   
   if (boundary.size() != supposed_boundary.size() ||
       bms.size() ||
       smb.size()) {
      passed = false;
      cerr << "\nnsw_d_spheres, constructibility: Lemma 2.2 failed with a mismatch between supposed_boundary and boundary:"
           << "\n| supposed_boundary - boundary:";
      for (const auto& s: smb) {
	cerr << "\n| " << s << ": ";
	for (const Int i: s)
           cerr << "(" << ij_labels[i] << ") ";
      }

      cerr << "\n| boundary - supposed_boundary:";
      for (const auto& s: bms) {
	cerr << "\n| " << s << ": ";
	for (const Int i: s)
           cerr << "(" << ij_labels[i] << ") ";
      }
     
     cerr << endl;
   }
   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Lemma 2.3
  *************************/
  

void
check_lemma_2_3(const dDBallData& bd,
                const Int verbosity,
                bool& passed)
{
   if (verbosity)
      cerr << "checking Lemma 2.3...";
   for (Int k=0; k < bd.k_max; ++k) {
     if (verbosity > 1)
       cerr << ".";
      for (const Simplex& sigma: bd.S_low_array[k]) {
         for (const IJIndex& ij: sigma.ij_indices) {
            const Int v(index_of_a_i_j(bd.n, ij.first+1, ij.second));
            if (!bd.delta_B_array[k].contains(sigma.vertices - scalar2set(v))) {
               cerr << "k=" << k
                    << ", sigma=" << sigma
                    << ", ij=" << ij
                    << ", v=" << v
                    << ", boundary=" << bd.delta_B_array[k]
                    << " should contain "<< (sigma.vertices - scalar2set(v))
                    << " but doesnt"
                    << endl;
            }
         }
      }

      for (const Simplex& sigma: bd.S_up_array[k]) {
         for (const IJIndex& ij: sigma.ij_indices) {
            const Int v(index_of_a_i_j(bd.n, ij.first, ij.second));
            if (!bd.delta_B_array[k].contains(sigma.vertices - scalar2set(v))) {
               passed = false;
               cerr << "k=" << k
                    << ", sigma=" << sigma
                    << ", ij=" << ij
                    << ", v=" << v
                    << ", boundary=" << bd.delta_B_array[k]
                    << " should contain "<< (sigma.vertices - scalar2set(v))
                    << " but doesnt"
                    << endl;
            }
         }
      }
   }
   if (verbosity)
      cerr << " done."
           << endl;
}

  /************************
    Lemma 3.1
  *************************/
  
Set<Int>
removed_ridge(const Simplex& sigma,
              const Int n,
              const Int ell,
              const Int m)
{
   assert(sigma.ij_indices[ell].first > 0);
   assert(sigma.ij_indices[m].first < n-2);
   assert(ell < m);
   
   Set<Int> remove;
   remove += index_of_a_i_j(n,
                            sigma.ij_indices[ell].first + 1,
                            sigma.ij_indices[ell].second);
   remove += index_of_a_i_j(n,
                            sigma.ij_indices[m].first,
                            sigma.ij_indices[m].second);

   if (!incl(remove, sigma.vertices)) {
      cerr << ", current simplex: " << sigma
           << ", ell=" << ell
           << ", m=" << m
           << ", remove=" << remove
           << ": remove not included in current simplex"
           << endl;
   }
   return sigma.vertices - remove;
}

void
check_lemma_3_1(const dDBallData& bd,
                const IJLabels& ij_labels,
                const Int verbosity,
                bool& passed)
{
   if (verbosity)
      cerr << "checking Lemma 3.1...";
   for (Int k=1; k <= bd.k_max; ++k) {
     if (verbosity > 1)
       cerr << ".";
      for (Int dp = (k-1)*(bd.d+2); dp <= k*(bd.d+2)-1; ++dp) {
         Set<Simplex> Dset;
         for (const Simplex& sigma: bd.B_array[k-1])
            if (dp == sigma.index_sum)
               Dset += sigma;
         if (Dset.size() < 2)
            continue;
         
         Set<Int> cumulative_union(Dset.front().vertices);
         auto it = entire(Dset);
         for (++it; !it.at_end(); ++it) {
            const Simplex sigma_j(*it);

            Set<Int> vertices_of_generated_complex;
            
            for (Int ell = 0; ell < bd.d-1; ++ell) {
               if (0 == sigma_j.ij_indices[ell].first) 
                  continue;

               for (Int m = ell+1; m < bd.d; ++m) {
                  if (bd.n-2 == sigma_j.ij_indices[m].first)
                     continue;
                  
                  vertices_of_generated_complex +=
                     removed_ridge(sigma_j, bd.n, ell, m);
               }
            }

            const Set<Int> verts_of_complex(sigma_j.vertices * cumulative_union);
            if (verts_of_complex != vertices_of_generated_complex) {
               passed = false;
               cerr << "k=" << k
                    << ", dp=" << dp
                    << ", Dset=" << Dset
                    << ", vertices_of_generated_complex=" << vertices_of_generated_complex
                    << ", sigma_j=" << sigma_j
                    << ", cumulative_union=" << cumulative_union
                    << ": sigma_j*cumulative_union=" << (sigma_j.vertices*cumulative_union)
                    << " differs from vertices_of_generated_complex=" << vertices_of_generated_complex
                    << endl;
            }
            cumulative_union += sigma_j.vertices;
         }
      }
   }
   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Lemma 3.2
  *************************/
  
void
check_lemma_3_2(const dDBallData& bd,
                const IJLabels& ij_labels,
                const Int verbosity,
                bool& passed)
{
   if (verbosity)
      cerr << "checking Lemma 3.2...";
   for (Int k=1; k <= bd.k_max; ++k) {
     if (verbosity > 1)
       cerr << ".";
      for (const Simplex& sigma: bd.S_low_array[k-1]) {
         Set<ShellingOrderedRidge> Delta;
         
         for (Int ell = 0; ell < bd.d-1; ++ell) {
            if (0 == sigma.ij_indices[ell].first) 
               continue;

            for (Int m = ell+1; m < bd.d; ++m) {
               if (bd.n-2 == sigma.ij_indices[m].first)
                  continue;
               Delta += ShellingOrderedRidge(sigma, bd.n, ell, m);
            }
         }
         auto it1_sentinel = Delta.end(); --it1_sentinel;
         for (auto it1 = entire(Delta); it1 != it1_sentinel; ++it1) {
            const auto& lp_mp = it1->ell_m;
            auto it2 = it1;
            for (++it2; !it2.at_end(); ++it2) {
               const auto& l_m = it2->ell_m;
               if (lp_mp.first == l_m.first) { // case 1
                  const Int mp(lp_mp.second);
                  const Int mp_index(index_of_a_i_j(bd.n, sigma.ij_indices[mp].first, sigma.ij_indices[mp].second));
                  if (it1->vertices * it2->vertices !=
                      it2->vertices - scalar2set(mp_index)) {
                     passed = false;
                     cerr << "Lemma 3.2 failed at point 1"
                          << "\n Delta=" << Delta
                          << "\n *it1=" << *it1
                          << "\n *it2=" << *it2
                          << "\n lp_mp=" << lp_mp
                          << "\n l_m=" << l_m
                          << "\n mp_index=" << mp_index
                          << "\n it1->vertices * it2->vertices=" << (it1->vertices * it2->vertices)
                          << "\n it2->vertices - scalar2set(mp_index)=" << (it1->vertices - scalar2set(mp_index))
                          << endl;
                  }
               } else if (lp_mp.first < l_m.first) {
                  const ShellingOrderedRidge F_lp_m(sigma, bd.n, lp_mp.first, l_m.second);
                  if (!Delta.contains(F_lp_m)) {
                     passed = false;
                     cerr << "Lemma 3.2 failed at point 2" << endl;
                  }

                  // It's clear by construction that (l',m) <= (l,m)

                  const SimplexSet intersect(it2->vertices * F_lp_m.vertices);
                  if (!(incl(it1->vertices * it2->vertices, intersect) <= 0)) {
                     passed = false;
                     cerr << "Lemma 3.2 failed at point 3"
                          << "\n Delta=" << Delta
                          << "\n *it1=" << *it1
                          << "\n *it2=" << *it2
                          << "\n lp_mp=" << lp_mp
                          << "\n l_m=" << l_m
                          << "\n F_lp_m=" << F_lp_m
                          << "\n it1->vertices * it2->vertices=" << (it1->vertices * it2->vertices)
                          << "\n intersect=" << intersect
                          << endl;
                  }
                  const Int lp(lp_mp.first);
                  const Int lp_index(index_of_a_i_j(bd.n, sigma.ij_indices[lp].first + 1, sigma.ij_indices[lp].second));
                  if (intersect != it2->vertices - scalar2set(lp_index)) {
                     passed = false;
                     cerr << "Lemma 3.2 failed at point 4" << endl;
                  }
               } else {
                  passed = false;
                  cerr << "Wrong ordering in Lemma 3.2" << endl;
               }
            }
         }
      }
   }
   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Definition 3.3
  *************************/
  
struct Def33Cmp {
   SimplexSet tau_r;
   Int i_r;
   Int m;
};

Def33Cmp
make_def33_cmp(const SimplexSet& tau,
               const Int r,
               const IJLabels& ij_labels)
{
   Def33Cmp def33_cmp;
   def33_cmp.i_r = 1000000000;
   def33_cmp.m = -1;
   for (const Int i: tau) {
      if (i < ij_labels.size() && // i could be an o_k
          r == ij_labels[i].second) {
         def33_cmp.tau_r += i;
         assign_min(def33_cmp.i_r, ij_labels[i].first);
         def33_cmp.m = i;
      }
   }
   return def33_cmp;
}

// enum ::pm::cmp_value { cmp_lt=-1, cmp_eq=0, cmp_gt=1, cmp_ne=cmp_gt };
::pm::cmp_value
def_3_3_gt(const SimplexSet& tau,
           const SimplexSet& taup,
           const Int r,
           const IJLabels& ij_labels)
{
   const Def33Cmp
      tau_cmp (make_def33_cmp(tau,  r, ij_labels)),
      taup_cmp(make_def33_cmp(taup, r, ij_labels));

   // trivial, not stated
   if (tau_cmp.tau_r == taup_cmp.tau_r)
      return ::pm::cmp_eq;

   // case 1
   if ( tau_cmp. tau_r.size() &&
       !taup_cmp.tau_r.size())
      return ::pm::cmp_gt;

   if (!tau_cmp. tau_r.size() &&
        taup_cmp.tau_r.size())
      return ::pm::cmp_lt;

   // case 2
   if (1 == tau_cmp .tau_r.size() &&
       1 == taup_cmp.tau_r.size())
     return operations::cmp()(tau_cmp.m, taup_cmp.m);

   // case 3
   if (2 == tau_cmp .tau_r.size() &&
       2 == taup_cmp.tau_r.size())
     return operations::cmp()(tau_cmp.i_r, taup_cmp.i_r);

   // case 4
   if (2 == tau_cmp .tau_r.size() &&
       1 == taup_cmp.tau_r.size() &&
       tau_cmp.i_r >= taup_cmp.m)
      return ::pm::cmp_gt;

   if (1 == tau_cmp .tau_r.size() &&
       2 == taup_cmp.tau_r.size() &&
       taup_cmp.i_r >= tau_cmp.m)
      return ::pm::cmp_lt;

   // case 5
   if (1 == tau_cmp .tau_r.size() &&
       2 == taup_cmp.tau_r.size() &&
       tau_cmp.m >= taup_cmp.i_r + 1)
      return ::pm::cmp_gt;

   if (2 == tau_cmp .tau_r.size() &&
       1 == taup_cmp.tau_r.size() &&
       taup_cmp.m >= tau_cmp.i_r + 1)
      return ::pm::cmp_lt;
   
   throw std::runtime_error("\nnsw_d_spheres: def_3_3_gt: inconclusive comparison");
}

  /************************
    Definition 3.4
  *************************/
  
// enum ::pm::cmp_value { cmp_lt=-1, cmp_eq=0, cmp_gt=1, cmp_ne=cmp_gt };
::pm::cmp_value
def_3_4_cmp(const SimplexSet& tau,
            const SimplexSet& taup,
            const IJLabels& ij_labels,
            const Int d)
{
   if (tau.size() != taup.size())
     cerr << "\nnsw_d_spheres: def_3_4_cmp: incomparable simplices" << endl;

   if (tau == taup)
      return ::pm::cmp_eq;

   for (Int r=0; r < d; ++r) {
      const ::pm::cmp_value cmp_result(def_3_3_gt(tau, taup, r, ij_labels));
      if (::pm::cmp_gt == cmp_result)
         return ::pm::cmp_gt;
      if (::pm::cmp_lt == cmp_result)
         return ::pm::cmp_lt;
   }
   throw std::runtime_error("\nnsw_d_spheres: def_3_4_cmp inconclusive");
}

/**********************************
    Lemma 3.5: error reporting
***********************************/

Set<Int>
cell_vertices(const dDBallData& bd,
              const SimplexSet& F,
              const Map<Vector<Int>, Int>& index_of_pt)
{
  std::vector<std::vector<Int>> pairs(bd.d);
  for (const Int i: F) {
    if (i >= bd.n * bd.d)
      continue;
    const Int coordinate(i/bd.n); // rounded down
    pairs[coordinate].push_back(i);
  }
  
  Vector<Int> limits(bd.d);
  for (Int c=0; c < bd.d; ++c)
    limits[c] = pairs[c].size();

  Set<Int> cv;
  for (MultiDimCounter<true> mdc(limits); !mdc.at_end(); ++mdc) {
    Vector<Int> coordinates(bd.d);
    for (Int c=0; c < bd.d; ++c)
      coordinates[c] = pairs[c][(*mdc)[c]] % bd.n;
    cv += index_of_pt[coordinates];
  }
  return cv;
}

std::string
comma_if_not_first(bool& first,
                   const std::string filler)
{
  if (!first)
    return filler;
  first = false;
  return "";
}
   
void
output_complex(const dDBallData& bd,
	       const SimplexSet& tau,
	       const std::vector<SimplexSet>& running_union,
	       const Set<SimplexSet>& generating_cells)
{
   BigObject obj("fan::PolyhedralComplex<Rational>");
   ListMatrix<Vector<Int>> vertices(0, bd.d+1);   
   std::stringstream ss;
   ss << "nsw_" << (2*bd.d - 1)
      << "_" << bd.n << "+";
   bool first(true);
   for (const Int i: tau)
      ss << comma_if_not_first(first, "_") << i;
   obj.set_name(ss.str());
   
   ss.str("");
   ss << "shelling complex from nsw sphere of dimension "
      << (2*bd.d - 1)
      << " on " << bd.d
      << " paths of length " << bd.n;
   obj.set_description(ss.str());
   
   Map<Vector<Int>,Int> index_of_pt;
   Int next_index(0);

   first = true;
   const Vector<Int> limits(bd.d, bd.n); // d copies of n
   for (MultiDimCounter<true> mdc(limits); !mdc.at_end(); ++mdc) {
      vertices /= *mdc;
      index_of_pt[*mdc] = next_index++;
   }
   obj.take("VERTICES") << (ones_vector<Rational>(vertices.rows()) | Matrix<Rational>(vertices));

   Array<std::string> vertex_labels(vertices.rows());
   auto vit = entire(vertex_labels);
   for (MultiDimCounter<true> mdc(limits); !mdc.at_end(); ++mdc, ++vit) 
      *vit = Label(*mdc, bd, ss).text;
   obj.take("VERTEX_LABELS") << vertex_labels;

   IncidenceMatrix<> maximal_polytopes(running_union.size(), vertices.rows());
   Set<Int> generating_indices;
   Int running_index(0);
   for (const SimplexSet& F: running_union) {
      maximal_polytopes[running_index] = cell_vertices(bd, F, index_of_pt);
      if (generating_cells.contains(F))
         generating_indices += running_index;
      ++running_index;
   }
   obj.take("MAXIMAL_POLYTOPES") << maximal_polytopes;
   obj.attach("GENERATING_INDICES") << generating_indices;
   obj.save(obj.name());

   cerr << "\nnsw_d_spheres: Output written to "
        << obj.name() << ".pcom"
	<< endl;   
}

  /************************
    Lemma 3.5
  *************************/

void
lemma_3_5_impl(const dDBallData& bd,
	       const Set<Def34OrderedSimplexSet>& ordered_set,
	       std::vector<SimplexSet>& running_union,
	       const Int verbosity,
               bool& passed,
               const bool output_on_error)
{
  std::stringstream ss;
  for (const Def34OrderedSimplexSet& tau_star_o: ordered_set) {
    // tau is tau_j in the lemma
    std::vector<Int> o;
    for (const Int v: tau_star_o.tau)
      if (v > bd.n * bd.d)
	o.push_back(v);
    if (o.size() > 1)
      cerr << "\nnsw_d_spheres: Lemma 3.5: Unexpected vertex in connecting path" << endl;
    const SimplexSet tau(o.size()
			 ? tau_star_o.tau - scalar2set(o.front())
			 : tau_star_o.tau);
    FacetList intersection_complex(bd.n);
    for (const SimplexSet& sigma: running_union)
      intersection_complex.insertMax(tau * sigma);

    Set<Int> cardinalities;
    for (const auto& f: intersection_complex)
      cardinalities += f.size();
    if (intersection_complex.size() &&
	1 != cardinalities.size()) {
       passed = false;
       Map<Vector<Int>, Int> index_of_pt;
       const Vector<Int> limits(bd.d, bd.n); // d copies of n

       if (verbosity > 1)
          cerr << "\n\n| during an attempted shelling, the simplex tau " << tau
               << " = " << Label(tau, bd, ss)
               << "\n| generated a mixed-dimensional intersection_complex:\n";

       Set<SimplexSet> generating_cells;
       for (const auto& s: intersection_complex) {
          bool printed_s(false);
          for (const auto& t: running_union)
             if (incl(s,t) <= 0) {
                if (verbosity > 1 && !printed_s) {
                   cerr << "| the intersection " << s << "\n";
                   printed_s = true;
                }
                if (verbosity > 1)
                   cerr << "| | was generated by " << t << " = " << Label(t, bd, ss) << endl;
                generating_cells += t;
             }
       }

       running_union.push_back(tau);
       cerr << "nsw_d_spheres, constructibility: Lemma 3.5 failed. ";
       if (output_on_error)
          output_complex(bd, tau, running_union, generating_cells);
       else
          cerr << endl;
    } else
       running_union.push_back(tau);
  }
  if (verbosity > 2) {
     Set<SimplexSet> generating_cells;
     generating_cells += running_union.back();
     output_complex(bd, ordered_set.back().tau, running_union, generating_cells);
  }
}

void
check_lemma_3_5(const dDBallData& bd,
                const Array<ModifiedDiagonals>& md,
                const Int verbosity,
                bool& passed,
                const bool output_on_error)
{
   if (verbosity)
      cerr << "checking Lemma 3.5...";      
   for (Int k=1; k <= bd.k_max; ++k) {
     if (verbosity > 1)
       cerr << ".";
      Set<Def34OrderedSimplexSet> ordered_connecting_path;      
      for (const SimplexSet& tau: md[k-1].connecting_path) 
         ordered_connecting_path += Def34OrderedSimplexSet(tau);

      std::vector<SimplexSet> running_union_complex;
      running_union_complex.reserve(bd.S_low_array[k-1].size() +
                                    md[k-1].connecting_path.size());
      for (const Simplex& sigma: bd.S_low_array[k-1])
         running_union_complex.push_back(sigma.vertices);
      
      lemma_3_5_impl(bd, ordered_connecting_path, running_union_complex, verbosity, passed, output_on_error);
   }

   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Lemma 3.6
  *************************/

Set<Int>
rest_case_1(const Int n,
            const Set<Int>& sigma_j_vertices,
            const IJIndex& ell_pair,
            const IJIndex& em_pair,
            bool& passed)
{
   Set<Int> rest(sigma_j_vertices);
   rest -= index_of_a_i_j(n, ell_pair.first + 1, ell_pair.second);
   rest -= index_of_a_i_j(n, em_pair.first + 1, em_pair.second);
   if (rest.size() != sigma_j_vertices.size() - 2) {
      passed = false;
      cerr << "\nnsw_d_spheres: Lemma 3.6 or Def 3.7 failed in case 1" << endl;
   }
   return rest;
}

void
lemma_3_6_case_1(Set<SimplexSet>& claimed_intersection_complex,
                 const Simplex& sigma_j,
                 const Int n,
                 bool& passed)
{
   for (auto ij_it_l = entire(sigma_j.ij_indices); !ij_it_l.at_end(); ++ij_it_l) {
      if (ij_it_l->first != 0)
         continue;
      for (auto ij_it_m = entire(sigma_j.ij_indices); !ij_it_m.at_end(); ++ij_it_m) {
         if (ij_it_m->first == 0)
            continue;
         claimed_intersection_complex += rest_case_1(n, sigma_j.vertices, *ij_it_l, *ij_it_m, passed);
      }
   }
}

Set<Int>
rest_case_2(const Int n,
            const Set<Int>& sigma_j_vertices,
            const IJIndex& ell_pair,
            bool& passed)
{
   Set<Int> rest(sigma_j_vertices);
   rest -= index_of_a_i_j(n, ell_pair.first,     ell_pair.second);
   rest -= index_of_a_i_j(n, ell_pair.first + 1, ell_pair.second);
   if (rest.size() != sigma_j_vertices.size() - 2) {
      passed = false;
      cerr << "\nnsw_d_spheres: Lemma 3.6 or Def 3.7 failed in case 2" << endl;
   }
   return rest;
}

void
lemma_3_6_case_2(Set<SimplexSet>& claimed_intersection_complex,
                 const Simplex& sigma_j,
                 const Int n,
                 bool& passed)
{
   for (auto ij_it_l = entire(sigma_j.ij_indices); !ij_it_l.at_end(); ++ij_it_l) {
      if (ij_it_l->first == 0)
         continue;
      claimed_intersection_complex += rest_case_2(n, sigma_j.vertices, *ij_it_l, passed);
   }
}

Set<Int>
rest_case_3(const Int n,
            const Set<Int>& sigma_j_vertices,
            const IJIndex& ell_pair,
            const IJIndex& em_pair,
            bool& passed)
{
   if (ell_pair.second == em_pair.second) {
      passed = false;
      cerr << "\nnsw_d_spheres: rest_case_3 called inappropriately" << endl;
   }

   Set<Int> rest(sigma_j_vertices);
   rest -= index_of_a_i_j(n, ell_pair.first,    ell_pair.second);
   rest -= index_of_a_i_j(n, em_pair. first + 1, em_pair.second);
   if (rest.size() != sigma_j_vertices.size() - 2) {
      passed = false;
      cerr << "\nnsw_d_spheres: Lemma 3.6 or Def 3.7 failed in case 3" << endl;
   }

   return rest;
}

void
lemma_3_6_case_3(Set<SimplexSet>& claimed_intersection_complex,
                 const Simplex& sigma_j,
                 const Int n,
                 bool& passed)
{
   for (auto ij_it_l = entire(sigma_j.ij_indices); !ij_it_l.at_end(); ++ij_it_l) {
      if (ij_it_l->first != n-2)
         continue;
      for (auto ij_it_m = entire(sigma_j.ij_indices); !ij_it_m.at_end(); ++ij_it_m) {
         if (ij_it_m->first == 0 ||
             ij_it_m->second == ij_it_l->second)
            continue;
         claimed_intersection_complex += rest_case_3(n, sigma_j.vertices, *ij_it_l, *ij_it_m, passed);
      }
   }
}

Set<Int>
rest_case_4(const Int n,
            const Set<Int>& sigma_j_vertices,
            const IJIndex& ell_pair,
            const IJIndex& em_pair,
            bool& passed)
{
   Set<Int> rest(sigma_j_vertices);
   rest -= index_of_a_i_j(n, ell_pair.first + 1, ell_pair.second);
   rest -= index_of_a_i_j(n, em_pair. first,     em_pair. second);
   if (rest.size() != sigma_j_vertices.size() - 2) {
      passed = false;
      cerr << "\nnsw_d_spheres: Lemma 3.6 or Def 3.7 failed in case 4" << endl;
   }

   return rest;
}

void
lemma_3_6_case_4(Set<SimplexSet>& claimed_intersection_complex,
                 const Simplex& sigma_j,
                 const Int n,
                 bool& passed)
{
   for (auto ij_it_l = entire(sigma_j.ij_indices); !ij_it_l.at_end(); ++ij_it_l) {
      if (ij_it_l->first == 0)
         continue;
      auto ij_it_m = ij_it_l;
      for (++ij_it_m; !ij_it_m.at_end(); ++ij_it_m) {
         if (ij_it_m->first == n-2 ||
             ij_it_m->second <= ij_it_l->second)
            continue;
         claimed_intersection_complex += rest_case_4(n, sigma_j.vertices, *ij_it_l, *ij_it_m, passed);
      }
   }
}

void
check_lemma_3_6(const dDBallData& bd,
                const Array<ModifiedDiagonals>& md,
                const Int verbosity,
                bool& passed)
{
   if (verbosity)
      cerr << "checking Lemma 3.6...";
   for (Int k=1; k <= bd.k_max; ++k) {
     if (verbosity > 1)
       cerr << ".";
      std::vector<SimplexSet> running_union_complex;
      running_union_complex.reserve(bd.S_low_array[k-1].size() +
                                    md[k-1].connecting_path.size() +
                                    bd.S_up_array[k-1].size());
      for (const Simplex& sigma: bd.S_low_array[k-1])
         running_union_complex.push_back(sigma.vertices);
      for (const SimplexSet& sigma: md[k-1].connecting_path)
         running_union_complex.push_back(sigma);
      
      for (const Simplex& sigma_j: bd.S_up_array[k-1]) {

         FacetList intersection_complex(bd.n);
         for (const SimplexSet& sigma: running_union_complex)
            intersection_complex.insertMax(sigma_j.vertices * sigma);
         
         Set<Int> cardinalities;
         for (const auto& f: intersection_complex)
            cardinalities += f.size();
         if (intersection_complex.size() &&
             1 != cardinalities.size()) {
            passed = false;
            cerr << "\nnsw_d_spheres: Lemma 3.6 failed at point 1 with"
		 << "\n| k=" << k
                 << "\n| S_low: " << bd.S_low_array[k-1]
                 << "\n| connecting_path: " << md[k-1].connecting_path
                 << "\n| S_up: " << bd.S_up_array[k-1]
                 << "\n| running_union_complex:\n" << running_union_complex
                 << "\n| intersection_complex:\n" << intersection_complex
                 << "\n| sigma_j: " << sigma_j
                 << endl;
         }

         Set<SimplexSet> claimed_intersection_complex;
         lemma_3_6_case_1(claimed_intersection_complex, sigma_j, bd.n, passed);
         lemma_3_6_case_2(claimed_intersection_complex, sigma_j, bd.n, passed);
         lemma_3_6_case_3(claimed_intersection_complex, sigma_j, bd.n, passed);
         lemma_3_6_case_4(claimed_intersection_complex, sigma_j, bd.n, passed);
         
         if (claimed_intersection_complex.size() !=
             intersection_complex.size()) {
            passed = false;
            cerr << "\nLemma 3.6 failed at point 2 with"
		 << "\n|  k=" << k
		 << "\n|  S_low: " << bd.S_low_array[k-1]
		 << "\n|  connecting_path: " << md[k-1].connecting_path
		 << "\n|  S_up: " << bd.S_up_array[k-1]
		 << "\n|  running_union_complex:\n" << running_union_complex
		 << "\n|  intersection_complex:\n" << intersection_complex
		 << "\n|  claimed_intersection_complex:\n" << claimed_intersection_complex
		 << "\n|  sigma_j: " << sigma_j
		 << endl;
         }
         
         for (const auto& s: intersection_complex) {
            if (!claimed_intersection_complex.contains(s)) {
               passed = false;
               cerr << "\nLemma 3.6 failed at point 3 with"
                    << "\n| k=" << k
                    << "\n|  S_low: " << bd.S_low_array[k-1]
                    << "\n|  connecting_path: " << md[k-1].connecting_path
                    << "\n|  S_up: " << bd.S_up_array[k-1]
                    << "\n|  running_union_complex:\n" << running_union_complex
                    << "\n|  intersection_complex:\n" << intersection_complex
                    << "\n|  claimed_intersection_complex:\n" << claimed_intersection_complex
                    << "\n|  sigma_j: " << sigma_j
                    << "\n|  s: " << s
                    << endl << endl;
            }
         }
         
         running_union_complex.push_back(sigma_j.vertices);
      }
   }

   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Definition 3.7
  *************************/

void
add_case_37_1(Set<ShellingOrderedSubridge38>& subridges,
              const Simplex& sigma_j,
              const Int ell,
              const Int n,
              const Int verbosity,
              bool& passed)
{
   Int r(0);
   for (auto imr_it = entire(sigma_j.ij_indices); !imr_it.at_end(); ++imr_it) 
      if (imr_it->first > 0) {
         const ShellingOrderedSubridge38 s(rest_case_1(n, sigma_j.vertices, sigma_j.ij_indices[ell], *imr_it, passed), ell, ++r); // r is 1-based
         if (verbosity > 3)
            cerr << "Def 3.7 case 1: " << s << endl;
         subridges += s; 
      }
}

void
add_case_37_2(Set<ShellingOrderedSubridge38>& subridges,
              const Simplex& sigma_j,
              const Int ell,
              const Int n,
              const Int verbosity,
              bool& passed)
{
   const ShellingOrderedSubridge38 s0(rest_case_2(n, sigma_j.vertices, sigma_j.ij_indices[ell], passed), ell, 0);
   subridges += s0;
   if (verbosity > 3)
      cerr << "Def 3.7 case 2 sigma_j = " << sigma_j
           << ", s0=" << s0 << endl;
   
   Int r(0);
   for (auto imr_it = entire(sigma_j.ij_indices); !imr_it.at_end(); ++imr_it) {
      if (imr_it->second > ell &&
          imr_it->first < n-2) {
         const ShellingOrderedSubridge38 s(rest_case_4(n, sigma_j.vertices, sigma_j.ij_indices[ell], *imr_it, passed), ell, ++r); // r is 1-based
         subridges += s;
         if (verbosity > 3)
            cerr << "Def 3.7 case 2 sigma_j = " << sigma_j
                 << ", imr = " << *imr_it
                 << ", s=" << s << endl;
      }
   }
}

void
add_case_37_3(Set<ShellingOrderedSubridge38>& subridges,
              const Simplex& sigma_j,
              const Int ell,
              const Int n,
              const Int verbosity,
              bool& passed)
{
   Int r(0);
   for (auto imr_it = entire(sigma_j.ij_indices); !imr_it.at_end(); ++imr_it) {
      if (imr_it->first > 0 &&
          imr_it->second != ell) {
         const ShellingOrderedSubridge38 s(rest_case_3(n, sigma_j.vertices, sigma_j.ij_indices[ell], *imr_it, passed), ell, ++r); // r is 1-based
         subridges += s;
         if (verbosity > 3)
            cerr << "Def 3.7 case 3.1 sigma_j = " << sigma_j
                 << ", s=" << s << endl;
      }
   }
   const Int R(r);
   if (verbosity > 3)
      cerr << "R = " << R << endl;

   const ShellingOrderedSubridge38 s0(rest_case_2(n, sigma_j.vertices, sigma_j.ij_indices[ell], passed), ell, R+1);
   subridges += s0;
   if (verbosity > 3)
      cerr << "Def 3.7 case 3.2 sigma_j = " << sigma_j
           << ", s0=" << s0 << endl;

   r = 0;
   for (auto isr_it = entire(sigma_j.ij_indices); !isr_it.at_end(); ++isr_it) {
      if (isr_it->second > ell &&
          isr_it->first < n-2) {
         const ShellingOrderedSubridge38 s(rest_case_4(n, sigma_j.vertices, sigma_j.ij_indices[ell], *isr_it, passed), ell, ++r + 1 + R); // r is 1-based
         subridges += s;
         if (verbosity > 3)
            cerr << "Def 3.7 case 3.3 sigma_j = " << sigma_j
                 << ", s=" << s << endl;
      }
   }
}

Set<ShellingOrderedSubridge38>
Def37OrderedSubridges(const Simplex& sigma_j,
                      const Int n,
                      const Int verbosity,
                      bool& passed)
{
   const Int d(sigma_j.ij_indices.size());
   Set<ShellingOrderedSubridge38> subridges;
   for (Int ell=0; ell < d; ++ell) {
      if (0 == sigma_j.ij_indices[ell].first) 
         add_case_37_1(subridges, sigma_j, ell, n, verbosity, passed);
      else if (sigma_j.ij_indices[ell].first < n-2)
         add_case_37_2(subridges, sigma_j, ell, n, verbosity, passed);
      else if (n-2 == sigma_j.ij_indices[ell].first)
         add_case_37_3(subridges, sigma_j, ell, n, verbosity, passed);
      else {
         passed = false;
         cerr << "\nnsw_d_spheres: Def 37: unexpected index tuple" << endl;
      }
   }
   return subridges;
}

  /************************
    Lemma 3.8
  *************************/

void
check_lemma_3_8(const dDBallData& bd,
                const Array<ModifiedDiagonals>& md,
                const Int verbosity,
                bool& passed)
{
   if (verbosity)
      cerr << "checking Lemma 3.8...";
   for (Int k=1; k <= bd.k_max; ++k) {
     if (verbosity > 1)
       cerr << ".";
      std::vector<SimplexSet> running_union_complex;
      running_union_complex.reserve(bd.S_low_array[k-1].size() +
                                    md[k-1].connecting_path.size() +
                                    bd.S_up_array[k-1].size());
      for (const Simplex& sigma: bd.S_low_array[k-1])
         running_union_complex.push_back(sigma.vertices);
      for (const SimplexSet& sigma: md[k-1].connecting_path)
         running_union_complex.push_back(sigma);

      for (const Simplex& sigma_j: bd.S_up_array[k-1]) {
         FacetList intersection_complex(bd.n);
         for (const SimplexSet& sigma: running_union_complex)
            intersection_complex.insertMax(sigma_j.vertices * sigma);

         Set<ShellingOrderedSubridge38> so_claimed_subridges(Def37OrderedSubridges(sigma_j, bd.n, verbosity, passed));
         Set<SimplexSet> set_claimed_subridges;
         for (const auto& s: so_claimed_subridges)
            set_claimed_subridges += s.vertices;

         if (set_claimed_subridges.size() !=
             intersection_complex.size()) {
            passed = false;
            cerr << "\nLemma 3.8 failed at point 1 with"
		 << "\n|  k=" << k
		 << "\n|  S_low: " << bd.S_low_array[k-1]
		 << "\n|  connecting_path: " << md[k-1].connecting_path
		 << "\n|  S_up: " << bd.S_up_array[k-1]
		 << "\n|  running_union_complex:\n" << running_union_complex
		 << "\n|  intersection_complex:\n" << intersection_complex
		 << "\n|  claimed_subridges:\n" << so_claimed_subridges
		 << "\n|  sigma_j: " << sigma_j
		 << endl << endl;
         }
         
         for (const auto& s: intersection_complex) {
            if (!set_claimed_subridges.contains(s)) {
               passed = false;
               cerr << "\nLemma 3.8 failed at point 2 with"
                    << "\n| k=" << k
                    << "\n| S_low: " << bd.S_low_array[k-1]
                    << "\n| connecting_path: " << md[k-1].connecting_path
                    << "\n| S_up: " << bd.S_up_array[k-1]
                    << "\n| running_union_complex:\n" << running_union_complex
                    << "\n| intersection_complex:\n" << intersection_complex
                    << "\n| claimed_subridges:\n" << so_claimed_subridges
                    << "\n| sigma_j: " << sigma_j
                    << "\n| s: " << s
                    << endl << endl;
            }
         }
         
         running_union_complex.push_back(sigma_j.vertices);

         for (auto Fpit = entire(so_claimed_subridges); !Fpit.at_end(); ++Fpit) {
            auto Fit = Fpit;
            for (++Fit; !Fit.at_end(); ++Fit) {
               const SimplexSet cap(Fit->vertices * Fpit->vertices);
               if (cap.size() == 2 * bd.d - 3)
                  continue;
               bool pass(false);
               for (auto Fppit = entire(so_claimed_subridges); !pass && !Fppit.at_end() && Fppit != Fit; ++Fppit) {
                  const SimplexSet cap2(Fit->vertices * Fppit->vertices);
                  if (cap2.size() != 2 * bd.d - 3)
                     continue;
                  if (incl(cap, cap2) <= 0)
                     pass = true;
               }
               if (!pass) {
                  passed = false;
                  cerr << "\nLemma 3.8 failed at point 3 with"
                       << "\n| k=" << k
                       << "\n| S_low: " << bd.S_low_array[k-1]
                       << "\n| connecting_path: " << md[k-1].connecting_path
                       << "\n| S_up: " << bd.S_up_array[k-1]
                       << "\n| running_union_complex:\n" << running_union_complex
                       << "\n| intersection_complex:\n" << intersection_complex
                       << "\n| claimed_subridges:\n" << so_claimed_subridges
                       << "\n| sigma_j: " << sigma_j
                       << "\n| F: " << *Fit
                       << "\n| Fp: " << *Fpit
                       << "\n| cap: " << cap
                       << endl << endl;
               }
            }
         }
      }
   }
   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Lemma 3.10
  *************************/

Set<SimplexSet>
F_sigma(const Simplex& sigma,
        const Int n)
{
   Set<SimplexSet> F_sigma;
   for (const auto& ell_pair: sigma.ij_indices) 
     if (ell_pair.first) // don't include left-lower boundaries
       F_sigma += sigma.vertices - scalar2set(index_of_a_i_j(n, ell_pair.first + 1, ell_pair.second));

   return F_sigma;
}

void
check_lemma_3_10(const dDBallData& bd,
                 const Array<ModifiedDiagonals>& md,
                 const Int verbosity,
                 bool& passed)
{
   if (verbosity)
      cerr << "checking Lemma 3.10...";
   // Tk = tilde S_k low ∪ tilde C_k ∪ tilde S_k up 

   Int total_size(0);
   for (Int k=1; k <= bd.k_max; ++k) 
      total_size +=
         md[k-1].S_low_tilde.size()
         + md[k-1].connecting_path.size()
         + md[k-1].S_up_tilde.size();
      
   std::vector<SimplexSet> running_union_complex;
   running_union_complex.reserve(total_size);
   for (const auto& s: md[0].S_low_tilde)
      running_union_complex.push_back(s);
   for (const auto& s: md[0].connecting_path)
      running_union_complex.push_back(s);
   for (const auto& s: md[0].S_up_tilde)
      running_union_complex.push_back(s);
   
   for (Int k=2; k <= bd.k_max; ++k) {
     if (verbosity > 1)
       cerr << ".";
      FacetList intersection_complex(bd.n);
      for (const auto& s: md[k-1].S_low_tilde)
         for (const SimplexSet& sigma: running_union_complex)
            intersection_complex.insertMax(s * sigma);
      for (const auto& s: md[k-1].connecting_path)
         for (const SimplexSet& sigma: running_union_complex)
            intersection_complex.insertMax(s * sigma);
      for (const auto& s: md[k-1].S_up_tilde)
         for (const SimplexSet& sigma: running_union_complex)
            intersection_complex.insertMax(s * sigma);

      Set<SimplexSet> claimed_intersection_complex;
      for (const auto& sigma: bd.S_low_array[k-1])
	claimed_intersection_complex += F_sigma(sigma, bd.n);

      if (claimed_intersection_complex.size() !=
          intersection_complex.size())  {
         passed = false;
         cerr << "\nLemma 3.10 failed at point 1 with"
              << "\n| k=" << k;
         for (Int k0 = 1; k0 <= k; ++k0)
            cerr << "\n|\n| S_low[" << k0-1 << "]: " << bd.S_low_array[k0-1]
		 << "\n| S_up [" << k0-1 << "]: " << bd.S_up_array[k0-1]
		 << "\n|\n| tilde S_low[" << k0-1 << "]:     " << md[k0-1].S_low_tilde
		 << "\n| connecting_path[" << k0-1 << "]: " << md[k0-1].connecting_path
		 << "\n| tilde S_up[" << k0-1 << "]:      " << md[k0-1].S_up_tilde;
	 cerr << "\n|\n| running_union_complex U T_k = U tilde S_low u C_k u tilde S_up, k=0.." << k-2
	      << " (omitted)"
	      << "\n| intersection_complex = running union cap T_" << k-1
	      << ":\n" << intersection_complex
	      << "\n| claimed_intersection_complex:\n" << claimed_intersection_complex
	      << "\n| intersection_complex - claimed_intersection_complex:\n" << (Set<SimplexSet>(entire(intersection_complex)) - claimed_intersection_complex)
	      << "\n| claimed_intersection_complex - intersection_complex:\n" << (claimed_intersection_complex - Set<SimplexSet>(entire(intersection_complex)))
	      << endl << endl;
      }

      for (const auto& s: intersection_complex) {
         if (!claimed_intersection_complex.contains(s)) {
            passed = false;
            cerr << "\nLemma 3.10 failed at point 2 with"
                 << "\n| k=" << k;
            for (Int k0 = 1; k0 <= k; ++k0)
	      cerr << "\n|\n| S_low[" << k0-1 << "]: " << bd.S_low_array[k0-1]
		   << "\n| S_up [" << k0-1 << "]: " << bd.S_up_array[k0-1]
		   << "\n|\n| tilde S_low[" << k0-1 << "]:     " << md[k0-1].S_low_tilde
		   << "\n| connecting_path[" << k0-1 << "]: " << md[k0-1].connecting_path
		   << "\n| tilde S_up[" << k0-1 << "]:      " << md[k0-1].S_up_tilde;
            cerr << "\n|\n| running_union_complex U T_k = U tilde S_low u C_k u tilde S_up, k=0.." << k-2
		 << ":\n" << running_union_complex
		 << "\n| intersection_complex = running union cap T_" << k-1
		 << ":\n" << intersection_complex
		 << "\n| claimed_intersection_complex:\n" << claimed_intersection_complex
		 << "\n| s: " << s
		 << endl << endl;
         }
      }
      
      for (const auto& s: md[k-1].S_low_tilde)
         running_union_complex.push_back(s);
      for (const auto& s: md[k-1].connecting_path)
         running_union_complex.push_back(s);
      for (const auto& s: md[k-1].S_up_tilde)
         running_union_complex.push_back(s);
   }

   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Lemma 3.11
  *************************/

void
check_lemma_3_11(const dDBallData& bd,
                 const Int verbosity,
                 bool& passed)
{
   if (verbosity)
      cerr << "checking Lemma 3.11...";

   for (Int k=2; k <= bd.k_max; ++k) { // S_low[0] is empty
     if (verbosity > 1)
       cerr << ".";

     const FacetsOfBall& S_low(bd.S_low_array[k-1]);

     std::vector<SimplexSet> running_union_sigmas, running_union_F_sigmas;
     running_union_sigmas  .reserve(S_low.size());
     running_union_F_sigmas.reserve(S_low.size());
     running_union_sigmas  .push_back(S_low.front().vertices);
     for (const auto& tau: F_sigma(S_low.front(), bd.n))
       running_union_F_sigmas.push_back(tau);

     auto S_low_it = entire(S_low);
     for (++S_low_it; !S_low_it.at_end(); ++S_low_it) {
       FacetList sigma_intersection(bd.n);
       for (const SimplexSet& u: running_union_sigmas)
	 sigma_intersection.insertMax(S_low_it->vertices * u);

       const Set<SimplexSet> F_sigma_j(F_sigma(*S_low_it, bd.n));

       FacetList F_sigma_intersection(bd.n);
       for (const SimplexSet& tau: F_sigma_j)
	 for (const SimplexSet& u: running_union_F_sigmas)
	   F_sigma_intersection.insertMax(tau * u);

       const Set<SimplexSet>
	 sigma_list(entire(sigma_intersection)),
	 F_sigma_list(entire(F_sigma_intersection));

       if (sigma_list != F_sigma_list) {
          passed = false;
          cerr << "\nLemma 3.11 failed at point 1 with"
               << "\n| k=" << k
               << "\n|\n| S_low[" << k-1 << "]: " << bd.S_low_array[k-1]
               << "\n|\n| running_union_sigmas: " << running_union_sigmas
               << "\n| sigma: " << *S_low_it
               << "\n| running_union_F_sigmas: " << running_union_F_sigmas
               << "\n| F_sigma: " << F_sigma_j
               << "\n| sigma_list - F_sigma_list:\n" << (sigma_list - F_sigma_list)
               << "\n| F_sigma_list - sigma_list:\n" << (F_sigma_list - sigma_list)
               << endl << endl;
       }
       
       running_union_sigmas.push_back(S_low_it->vertices);
       for (const SimplexSet& tau: F_sigma_j)
	 running_union_F_sigmas.push_back(tau);
     }
   }
  
   if (verbosity)
      cerr << " done." << endl;
}

  /************************
    Theorem 1.1(1)
  *************************/

void
check_Thm_1_1_1(const dDBallData& bd,
		const Array<ModifiedDiagonals>& md,
		const Int verbosity,
                bool& passed,
                const bool output_on_error)
{
   if (verbosity)
      cerr << "checking Theorem 1.1(1)...";
   
   // Tk = tilde S_k low ∪ tilde C_k ∪ tilde S_k up
   Set<Def34OrderedSimplexSet> boundary;
   for (Int k=1; k <= bd.k_max; ++k) {
     for (const SimplexSet& sigma: md[k-1].S_low_tilde) 
       add_to_boundary(sigma, boundary);
     for (const SimplexSet& sigma: md[k-1].S_up_tilde) 
       add_to_boundary(sigma, boundary);
     for (const SimplexSet& sigma: md[k-1].connecting_path) 
       add_to_boundary(sigma, boundary);
   }

   std::vector<SimplexSet> running_union;
   running_union.reserve(boundary.size());
   
   lemma_3_5_impl(bd, boundary, running_union, verbosity, passed, output_on_error);

   if (verbosity)
      cerr << " done." << endl;
   
}

/****************************
    Putting it all together
*****************************/

bool
check_constructibility_proof(const dDBallData& bd,
			     const Array<ModifiedDiagonals>& md,
                             const IJLabels& ij_labels,
                             const Int verbosity,
                             const bool output_on_error)
{
   bool passed(true);
   check_lemma_2_3(bd, verbosity, passed);
   check_lemma_3_1(bd, ij_labels, verbosity, passed);
   check_lemma_3_2(bd, ij_labels, verbosity, passed);
  
   global_d = bd.d;
   global_ij_label_ptr = &ij_labels;

   check_lemma_3_5(bd, md, verbosity, passed, output_on_error);
   check_lemma_3_6(bd, md, verbosity, passed);
   check_lemma_3_8(bd, md, verbosity, passed);
   check_lemma_3_10(bd, md, verbosity, passed);
   check_lemma_3_11(bd, verbosity, passed);
   check_Thm_1_1_1(bd, md, verbosity, passed, output_on_error);

   return passed;
}

} // end namespace nsw_sphere

}}


// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
    
