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
#include "polymake/MultiDimCounter.h"
#include "polymake/Rational.h"
#include "polymake/FacetList.h"
#include "polymake/topaz/nsw_d_spheres.h"

namespace polymake { namespace topaz {

namespace nsw_sphere {

   
dDBallData   
ball_data_for_dd_join_of_paths(const Int n,
                               const Int d,
                               const Int k_max,
                               const Int verbosity)
{
   if (verbosity > 3)
      cerr << "making ball data... ";
   
   dDBallData bd(n, d, k_max);
   const Vector<Int> limits(d, n-1); // d copies of n-1,
   // so that only indices from 0,...,n-2 <=> 1,...,n-1 are considered
   for (MultiDimCounter<false> mdc(limits); !mdc.at_end(); ++mdc) {
      const Int index_sum(d + accumulate(*mdc, operations::add()));
      // The "d" is to make the index_sum be 1-based like in the paper
      /*
        determine k in {1,...,k_max} such that
          (k-1)(d+2) <= index_sum <= k(d+2) - 1,
        equivalently, that
          k-1 <= index_sum/(d+2)     <= k    - 1/(d+2),
        which implies that
          k   <= index_sum/(d+2) + 1 <  k+1 
        equivalently, 
          k = floor(index_sum/(d+2) + 1)

        Sometimes, this k can exceed k_max in boundary cases,
        for example for d=2, n=7, or for d=4, n=10, in which case
        k_max = ceil(d(n-1)/(d+2)) = 6,
        but mdc = (8,8,8,8) has index sum 36, with floor(36/6 + 1) = 7
        and 36 is not less-or-equal to 6(6) - 1
      */
      const Int k(//std::min(Integer(k_max),
                           floor(Rational(index_sum + d+2, d+2)));
      if (k > k_max)
         continue;
      
      // construct and store the simplex
      // *mdc is an IndexTuple
      const Simplex s = bd.ss.simplex(*mdc);
      if (s.vertices.size() != 2*d) // it's a boundary case
         continue;
      
      bd.B_array[k-1] += s;

      // and its boundary
      add_to_boundary(s, bd.delta_B_array[k-1]);

      if (index_sum == (k-1)*(d+2))
         bd.S_low_array[k-1] += s;

      if (index_sum == k*(d+2) - 1)
         bd.S_up_array[k-1] += s;
   }
   if (verbosity > 3)
      cerr << "done" << endl;
   return bd;
}

Set<SimplexSet>
C_sigma_tilde_of_impl(const SimplexSet& F_sigma,
                      const Int o_k,
                      const SimplexSet& verts_of_D_sigma,
                      const TriangulationChoice& choice)
{

   const SimplexSet G_sigma(verts_of_D_sigma - F_sigma);
   if (!G_sigma.size())
      throw std::runtime_error("nsw_d_spheres: C_sigma_tilde_of_impl: unexpected F_sigma");
      
   const SimplexSet G_sigma_join_o_k(G_sigma + scalar2set(o_k));

   Set<SimplexSet> C_sigma_tilde;

   if (TriangulationChoice::two == choice)
      //k < choices.size() && choices[k])
      // Tσ,2 = ∂F_σ ∗ (G_σ ∗ o_k)
      for (auto it = entire(all_subsets_less_1(F_sigma)); !it.at_end(); ++it) 
         C_sigma_tilde += *it + G_sigma_join_o_k;
   else 
      // Tσ,1 = F_σ ∗ ∂(G_σ ∗ o_k)
      for (auto it = entire(all_subsets_less_1(G_sigma_join_o_k)); !it.at_end(); ++it) 
         C_sigma_tilde += *it + F_sigma;
   
   return C_sigma_tilde;
}
   
// as a side effect, calculates bd.union_of_Ds_array[k]
Set<SimplexSet>
C_sigma_tilde_of(const SimplexSet& sigma,
                       dDBallData& bd,
                 const TriangulationChoice& choice,
                 const Int k)
{
   const Set<BoundaryFacet> D_sigma(D_sigma_of(sigma, bd.delta_B_array[k]));
   bd.union_of_Ds_array[k] += D_sigma;
   
   SimplexSet verts_of_D_sigma;
   for (const auto& tau: D_sigma)
      verts_of_D_sigma += tau;
   
   const SimplexSet F_sigma(missing_face_of(verts_of_D_sigma, bd.delta_B_array[k]));   
   const Int o_k(index_of_new_vertex(bd.n, bd.d, k));

   return C_sigma_tilde_of_impl(F_sigma, o_k, verts_of_D_sigma, choice);
}

// as a side effect, calculates bd.union_of_Ds_array[k]
ModifiedDiagonals
modified_triangulation(const Int k,
                       const TriangulationChoice& choice,
                             dDBallData& bd)
{
   ModifiedDiagonals md;
   md.choice = choice;
   
   for (const Simplex& sigma: bd.S_low_array[k]) 
      md.S_low_tilde += C_sigma_tilde_of(sigma.vertices, bd, choice, k);

   for (const Simplex& sigma: bd.S_up_array[k]) 
      md.S_up_tilde += C_sigma_tilde_of(sigma.vertices, bd, choice, k);

   const Int o_k(index_of_new_vertex(bd.n, bd.d, k));
   for (const SimplexSet& tau: bd.delta_B_array[k]) 
      if (!bd.union_of_Ds_array[k].contains(tau))
         md.connecting_path += tau + scalar2set(o_k);
   
   return md;
}

std::tuple<Array<std::string>, IJLabels, IndexOfLabel>
make_labels(const Int n,
            const Int d,
            const Int k_max,
            const Int verbosity)
{
   if (verbosity > 3)
      cerr << "making labels...";
   
   Array<std::string> string_labels(convert_to<Int>(d*n + ceil(Rational(d*(n+1), d+2))));
   IndexOfLabel iol;
   IJLabels ij_labels(n*d);
   for (Int j=0; j<d; ++j)
      for (Int i=0; i<n; ++i) {
         const Int index(index_of_a_i_j(n, i, j));
         string_labels[index] =
            std::string("a") + std::to_string(i+1)
            + "(" + std::to_string(j+1) + ")";
         ij_labels[index] = std::make_pair(i, j);
         iol[std::make_pair(i,j)] = index;
      }
   
   for (Int k=0; k < k_max; ++k) 
      string_labels[index_of_new_vertex(n,d,k)] = "o" + std::to_string(k+1);

   string_labels[string_labels.size()-1] = "O";
   if (verbosity > 3)
      cerr << "done" << endl;
   return std::make_tuple(string_labels, ij_labels, iol);
}

template<typename Output>
Output&
print_labeled(GenericOutput<Output>& outs,
              const SimplexSet& s,
              const Array<std::string>& labels)
{
   outs.top() << "{ ";
   for (const auto& i: s)
      outs.top() << labels[i] << " ";
   return outs.top() << "}";
}
   
template<typename Output>
Output&
print_labeled(GenericOutput<Output>& outs,
              const dDBallData& bd,
              const Array<std::string>& labels)
{
   Output& os = outs.top();
   for (Int k=0; k < bd.k_max; ++k) {
      os << k+1 << ": ";
      for (const auto& s: bd.B_array[k])
         print_labeled(outs, s, labels);
      os << endl;

      os << "bd " << k+1 << ": ";
      for (const auto& s: bd.delta_B_array[k])
         print_labeled(outs, s, labels);
      os << endl;

      os << "S_low_array " << k+1 << ": ";
      for (const auto& s: bd.S_low_array[k])
         print_labeled(outs, s, labels);
      os << endl;

      os << "S_up_array " << k+1 << ": ";
      for (const auto& s: bd.S_up_array[k])
         print_labeled(outs, s, labels);
      os << endl;

      os << "union_of_Ds_array " << k+1 << ": ";
      for (const auto& s: bd.union_of_Ds_array[k])
         print_labeled(outs, s, labels);
      os << endl;
   }
   return os;
}


} // end namespace nsw_sphere


using namespace nsw_sphere;

BigObject
bistellar_d_sphere(const Int D,
                   const Int n,
                   OptionSet options)
{
   if (n<3) throw std::runtime_error("need n>=3");
   if (D<3 || !(D%2)) throw std::runtime_error("need D>=3 odd");
   
   // In order to avoid as many off-by-one errors as possible,
   // we will use the same 1-based indexing as the paper,
   // and let the function index_of_a_i_j take care of converting
   // indices to consecutive-zero-based

   const Int verbosity = options["verbosity"];
   const Int the_i = options["i"];
   const bool check_constructibility = options["check_constructibility"];
   const bool output_on_error = options["output_on_error"];
   
   const Int d ((D+1)/2);
   const Int k_max(ceil(Rational(d*(n-1), d+2)));
   const Int effective_k_max = std::min(Int(32), k_max);

   Vector<Int> choices(effective_k_max);
   for (Int i=0; i < effective_k_max; ++i)
      if (the_i & (Int(1) << i))
         choices[i] = 1;

   if (verbosity) {
      cerr << "dimension of sphere D=" << D
           << ", number of vertices n=" << n
           << ", number of paths d=" << d
           << ", k_max=" << k_max
           << ", serial number i=" << the_i
           << " (in reverse binary: ";

      for (const auto& i: choices)
         cerr << i;
      cerr << ")"
           << endl;
   }
   
   dDBallData bd(ball_data_for_dd_join_of_paths(n, d, k_max, verbosity));
   auto llm = make_labels(n, d, bd.k_max, verbosity);
   const Array<std::string>& labels = std::get<0>(llm);
         IJLabels& ij_labels = std::get<1>(llm);
   const IndexOfLabel& iol = std::get<2>(llm);
   if (verbosity) {
      cerr << "labels(" << labels.size() << "): ";
      for (Int i=0; i < labels.size(); ++i)
         cerr << i << ":" << labels[i] << " ";
      cerr << endl;
   }
   if (verbosity > 3)
      cerr << "ball_data:\n" << bd
           << "\ndone with ball data"
           << endl;

   bool passed_lemma_2_2(true);
   if (check_constructibility)
      nsw_sphere::check_lemma_2_2(iol, bd, ij_labels, verbosity, passed_lemma_2_2);
   
   FacetsOfBoundary boundary;
   Set<Set<Int>> facets;
   Array<ModifiedDiagonals> md(bd.k_max);
   for (Int k=0; k < bd.k_max; ++k) {
      // if k > 32, clamp the choice to one/false
      const TriangulationChoice choice(k < choices.size() && choices[k]
                                       ? TriangulationChoice::two
                                       : TriangulationChoice::one);
      md[k] = modified_triangulation(k, choice, bd);
      if (verbosity > 3) 
         cerr << "\n" << k+1 << ":\nS_low_tilde:\n";

      for (const auto& s: md[k].S_low_tilde) {
         if (verbosity > 3)
            print_labeled(cerr, s, labels);
         add_to_boundary(s, boundary);
         facets += s;
      }

      if (verbosity > 3)
         cerr << "\nS_up_tilde:\n";
         
      for (const auto& s: md[k].S_up_tilde) {
         if (verbosity > 3)
            print_labeled(cerr, s, labels);
         add_to_boundary(s, boundary);
         facets += s;
      }

      if (verbosity > 3)
         cerr << "\nconnecting_path:\n";

      for (const auto& s: md[k].connecting_path) {
         if (verbosity > 3)
            print_labeled(cerr, s, labels);
         add_to_boundary(s, boundary);
         facets += s;
      }

      if (verbosity > 3)
         cerr << endl;
   }
   
   for (const auto& b: boundary) 
      facets += b + scalar2set(labels.size()-1);
   
   BigObject s("SimplicialComplex");
   s.set_description() << "Nevo-Santos-Wilson "
                       << (2*d-1) << "-sphere on "
                       << labels.size() << " vertices, made from "
                       << d << " paths of length "
                       << n << ", with i="
                       << the_i << " corresponding to choice vector "
                       << choices
                       << ((0 == (d*(n-1)) % (d+2))
                           ? ". Cannot currently certify constructibility."
                           : "");
   s.take("FACETS") << facets;
   s.take("VERTEX_LABELS") << labels;

   if (check_constructibility &&
       check_constructibility_proof(bd, md, ij_labels, verbosity, output_on_error) &&
       passed_lemma_2_2) // put this at the end so that all checks run even if this one failed
      s.take("CONSTRUCTIBLE") << true;

   // This could be added at some point...
   // if (2 == bd.d)
   //    check_shellability_proof(bd, md, ij_labels, verbosity);
   
   return s;
}

UserFunction4perl("# @category Producing from scratch"
                  "# Create a (D = 2d-1)-sphere made from d paths of n vertices from"
                  "# [Nevo, Santos & Wilson, Many triangulated odd-dimensional spheres, Math Ann 364 (2016), 737-762."
                  "# @param Int D the dimension of the sphere, an integer >= 2"
                  "# @param Int n the number of vertices along a path, an integer >= 3"
                  "# @option Int verbosity default 0"
                  "# @option Int i the serial number of which triangulation to choose, "
                  "# where 0 <= i <= min(2^k_max - 1, 2^32 - 1), k_max = floor(d(n-1)/(d+2)), and d=(D+1)/2 the number of paths."
                  "# The value of i will be clamped to that range; default is 0"
                  "# @option Bool check_constructibility default 0:"
                  "# check that the sphere is constructible according to the lemmata in Yirong Yang, https://arxiv.org/abs/2305.03186."
                  "# The proof in that paper currently has an error whenever d+2 divides d(n-1), e.g. for (D,d,n) = (5,3,11), (7,4,10)."
                  "# @option Bool output_on_error default 1 output instances of the failed shellings in those cases" 
                  "# @return SimplicialComplex",
                  &bistellar_d_sphere, "bistellar_d_sphere($$ { verbosity=>0, i=>0, check_constructibility=>0, output_on_error=>1 })");


} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
    
    
