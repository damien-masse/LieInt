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

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::inverse()
const
{
   if (this->empty) return LieBaseMatrix<DimGroup,DimMatrix>::Empty();
   IntervalMatrix inv = inverse_enclosure(this->value);
   return LieBaseMatrix<DimGroup,DimMatrix>(inv);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
void LieBaseMatrix<DimGroup,DimMatrix>::inverse_inplace()
{
   if (this->empty) return;
   this->value = inverse_enclosure(this->value);
   this->contract();
}


template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::IleftProd
		(const LieBaseMatrix &A) const
{
   if (this->empty || A.is_empty()) return LieBaseMatrix<DimGroup,DimMatrix>::Empty();
   IntFullPivLU pivLU(this->value);
   IntervalMatrix RA = pivLU.solve(A);
   return LieBaseMatrix<DimGroup,DimMatrix>(RA);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::IrightProd
		(const LieBaseMatrix &A) const
{
   if (this->empty || A.is_empty()) return LieBaseMatrix<DimGroup,DimMatrix>::Empty();
   IntFullPivLU pivLU(this->value.transpose());
   IntervalMatrix RA = pivLU.solve(A.transpose());
   return LieBaseMatrix<DimGroup,DimMatrix>(RA.transpose());
}

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::center()
{
   if (this->empty) return LieBaseMatrix<DimGroup,DimMatrix>::Empty();
   Matrix Ct = this->value.mid();
   IntFullPivLU pivLU(Ct.transpose());
   if (pivLU.isInvertible()==BoolInterval::TRUE) {
      this->value = pivLU.solve(this->value.transpose()).transpose();
      this->contract();
      return LieBaseMatrix<DimGroup,DimMatrix>(Ct);
   } else {
      return LieBaseMatrix<DimGroup,DimMatrix>();
   }
}

template <unsigned int DimGroup, unsigned int DimMatrix>
std::ostream& operator<<(std::ostream& os,
                             const LieBaseMatrix<DimGroup,DimMatrix>& x) {
     if (x.empty) { os << "Lie:empty" ; return os; }
     os << "Lie:" << x.value;
     return os;
}

}
