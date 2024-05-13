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
#include "polymake/PowerSet.h"

namespace polymake { namespace polytope {


namespace {

Int col_idx(const Int n, Int i, Int j)
// index of the set {i,j} in the lex ordering, where {0,1} is the first and {n-2,n-1} is the last
{
  if (i>j) std::swap(i,j);
  return (i*(2*n-i-1))/2 + (j-i-1);
}
  
Matrix<Rational> triangle_ineqs(const Int n)
{
  const Int nrows(3*Integer::binom(n,3)), ncols(Integer::binom(n,2));
  Matrix<Rational> Ineqs(nrows, ncols);
  Int r = 0;
  for (const auto ijk : all_subsets_of_k(sequence(0,n), 3)) {
    for (const Int k : ijk) {
      const Set<Int> ij = ijk - k;
      const Int i = ij.front();
      const Int j = ij.back();
      // homogeneous: -x_{ij} + x_{ik} + x_{jk} >= 0
      Ineqs(r, col_idx(n,i,j)) = -1;
      Ineqs(r, col_idx(n,i,k)) = Ineqs(r, col_idx(n,j,k)) = 1;
      ++r;
    }
  }
  return Ineqs;
}

Matrix<Rational> bounding_ineqs(const Int n)
{
  const Int nrows(Integer::binom(n,3)), ncols(1+Integer::binom(n,2));
  Matrix<Rational> Ineqs(nrows, ncols);
  Int r = 0;
  for (Int i=0; i<n; ++i) {
    for (Int j=i+1; j<n; ++j) {
      for (Int k=j+1; k<n; ++k) {
        // inhomogeneous: x_{ij} + x_{ik} + x_{jk} <= 2
        Ineqs(r, 0) = 2;
        Ineqs(r, 1+col_idx(n,i,j)) = -1;
        Ineqs(r, 1+col_idx(n,i,k)) = -1;
        Ineqs(r, 1+col_idx(n,j,k)) = -1;
        ++r;
      }
    }
  }
  return Ineqs;
}

}

BigObject metric_cone(const Int n)
{
  if (n<3) throw std::runtime_error("metric_cone: n>=3 required");
  const Int d(Integer::binom(n,2));
  const Matrix<Rational> Ineqs(triangle_ineqs(n));
  const Matrix<Rational> NoLinealities(0, d);
  const Vector<Rational> all_ones(ones_vector<Rational>(d));
  return BigObject("Cone<Rational>",
                   "FACETS", Ineqs,
                   "LINEALITY_SPACE", NoLinealities,
                   "CONE_AMBIENT_DIM", d,
                   "CONE_DIM", d,
                   "REL_INT_POINT", all_ones);
}


BigObject metric_polytope(const Int n)
{
  if (n<3) throw std::runtime_error("metric_polytope: n>=3 required");
  const Int d(Integer::binom(n,2));
  const Matrix<Rational> Ineqs((zero_vector<Rational>(Int(3*Integer::binom(n,3))) | triangle_ineqs(n)) / bounding_ineqs(n));
  const Matrix<Rational> NoLinealities(0, 1+d);
  const Vector<Rational> all_ones(ones_vector<Rational>(1)|Rational(1,2)*ones_vector<Rational>(d));
  return BigObject("Polytope<Rational>",
                   "FACETS", Ineqs,
                   "AFFINE_HULL", NoLinealities,
                   "CONE_AMBIENT_DIM", d+1,
                   "CONE_DIM", d+1,
                   "REL_INT_POINT", all_ones);
}

UserFunction4perl("#@category Finite metric spaces"
                  "# Computes the metric cone on for points.  The triangle inequalities define the facets."
                  "# The number of rays are known for n <= 8."
                  "# See Deza and Dutour-Sikiric (2018), doi:10.1016/j.jsc.2016.01.009"
                  "# @param Int n"
                  "# @return Cone"
                  "# @example The number of rays of the metric cone for 4 points."
                  "# > print metric_cone(4)->N_RAYS"
                  "# | 7",
                  &metric_cone, "metric_cone($)");
      
UserFunction4perl("#@category Finite metric spaces"
                  "# Computes the metric polytope on for points.  This is the metric cone bounded by one inequality per triplet of points."
                  "# The number of vertices are known for n <= 8."
                  "# See Deza and Dutour-Sikiric (2018), doi:10.1016/j.jsc.2016.01.009"
                  "# @param Int n"
                  "# @return Polytope"
                  "# @example The volume of the metric polytope for 4 points."
                  "# > print metric_polytope(4)->VOLUME"
                  "# | 2/45",
                  &metric_polytope, "metric_polytope($)");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
