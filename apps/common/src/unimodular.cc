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
#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Set.h"
#include "polymake/Rational.h"
#include "polymake/linalg.h"
#include "polymake/integer_linalg.h"


namespace polymake { namespace common {

namespace {

// The bool `common_linear_subspace` indicates whether all input cells live in
// the same linear subspace. If this is the case, all cells can be transformed
// using the same transformation, reducing the number of calls to
// `hermite_normal_form` to one.
template<typename Func>
bool inner(const Matrix<Rational>& C_in, const Array<Set<Int>>& F, bool common_linear_subspace, Func f){
   if(common_linear_subspace){
      const Matrix<Integer> C(convert_to<Integer>(C_in));
      auto hnf = hermite_normal_form(C.minor(basis_rows(C),All));
      Int rk = hnf.rank;
      const auto transform = hnf.companion.minor(All, sequence(0,hnf.rank));

      for (const auto& fi : F) {
         bool is_not_unimodular = fi.size() != rk || abs( det(C.minor(fi,All)*transform) ) != 1;
         if (f(is_not_unimodular)) {
            return false;
         }
      }
      return true;
   } else {
      const Matrix<Integer> C(convert_to<Integer>(C_in));

      for (const auto& fi : F) {
         auto hnf = hermite_normal_form(C.minor(fi, All));
         bool is_not_unimodular = fi.size() != hnf.rank;
         for(Int i=0; i<hnf.rank && !is_not_unimodular; i++){
            if(!abs_equal(hnf.hnf(i,i),1)){
               is_not_unimodular = true;
               break;
            }
         }
         if (f(is_not_unimodular)) {
            return false;
         }
      }
      return true;
   }
}

}

// coordinates are homogeneous, maximal cells not necessarily simplices
// returns false if not a triangulation
bool unimodular(const Matrix<Rational>& C_in, const Array<Set<Int>>& F, bool common_linear_subspace)
{
   auto f = [](bool b){return b;};
   return inner(C_in, F, common_linear_subspace, f);
}


Int n_unimodular(const Matrix<Rational>& C_in, const Array<Set<Int>>& F, bool common_linear_subspace)
{
   Int n_unimodular = F.size();
   auto f = [&n_unimodular](bool b){n_unimodular-=b; return false;};
   inner(C_in, F, common_linear_subspace, f);
   return n_unimodular;
}

Function4perl(&unimodular,"unimodular");
Function4perl(&n_unimodular,"n_unimodular");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
