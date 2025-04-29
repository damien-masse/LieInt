/** 
 *  \file lieGroup.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <codac>
#include "lieGroup.h"

using namespace codac2;

namespace lieInt
{

#if 0
template <unsigned int DimGroup, unsigned int DimMatrix>
IntervalMatrix LieBaseMatrix<DimGroup,DimMatrix>::representationAlgebra
		(const IntervalVector &V)
{
   assert(false);
   return IntervalMatrix::Zero(DimMatrix,DimMatrix);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
Matrix LieBaseMatrix<DimGroup,DimMatrix>::representationAlgebra
		(const Vector &V)
{
   assert(false);
   return IntervalMatrix::Zero(DimMatrix,DimMatrix);
}
#endif

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
Derived LieBaseMatrix<DimGroup,DimMatrix,Derived>::inverse()
const
{
   if (this->empty) return LieBaseMatrix<DimGroup,DimMatrix,Derived>::Empty();
   IntervalMatrix inv = inverse_enclosure(this->value);
   return LieBaseMatrix<DimGroup,DimMatrix,Derived>(inv).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
void LieBaseMatrix<DimGroup,DimMatrix,Derived>::inverse_inplace()
{
   if (this->empty) return;
   this->value = inverse_enclosure(this->value);
   this->contract();
}


template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
Derived LieBaseMatrix<DimGroup,DimMatrix,Derived>::IleftProd
		(const LieBaseMatrix<DimGroup,DimMatrix,Derived> &A) const
{
   if (this->empty || A.is_empty()) return LieBaseMatrix<DimGroup,DimMatrix,Derived>::Empty();
   IntFullPivLU pivLU(this->value);
   IntervalMatrix RA = pivLU.solve(A);
   return LieBaseMatrix<DimGroup,DimMatrix,Derived>(RA).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
Derived LieBaseMatrix<DimGroup,DimMatrix,Derived>::IrightProd
		(const LieBaseMatrix<DimGroup,DimMatrix,Derived> &A) const
{
   if (this->empty || A.is_empty()) return LieBaseMatrix<DimGroup,DimMatrix,Derived>::Empty();
   IntFullPivLU pivLU(this->value.transpose());
   IntervalMatrix RA = pivLU.solve(A.transpose());
   return LieBaseMatrix<DimGroup,DimMatrix,Derived>(RA.transpose()).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
Derived LieBaseMatrix<DimGroup,DimMatrix,Derived>::center(bool uncertLeft)
{
   if (this->empty) return LieBaseMatrix<DimGroup,DimMatrix,Derived>::Empty();
   Matrix Ct = this->value.mid();
   if (uncertLeft) {
      IntFullPivLU pivLU(Ct.transpose());
      if (pivLU.isInvertible()==BoolInterval::TRUE) {
         this->value = pivLU.solve(this->value.transpose()).transpose();
         this->contract();
      } else {
         return LieBaseMatrix<DimGroup,DimMatrix,Derived>::Identity();
      }
   } else  {
      IntFullPivLU pivLU(Ct);
      if (pivLU.isInvertible()==BoolInterval::TRUE) {
         this->value = pivLU.solve(this->value());
         this->contract();
      } else {
         return LieBaseMatrix<DimGroup,DimMatrix,Derived>::Identity();
      }
   }
   return LieBaseMatrix<DimGroup,DimMatrix,Derived>(Ct).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
std::ostream& operator<<(std::ostream& os,
                             const LieBaseMatrix<DimGroup,DimMatrix,Derived>& x) {
     if (x.empty) { os << "Lie:empty" ; return os; }
     os << "Lie:" << x.value;
     return os;
}

}
