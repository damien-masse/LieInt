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

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix>::LieBaseMatrix() :
    value(IntervalMatrix::Identity(DimMatrix,DimMatrix))
{
}   

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix>::LieBaseMatrix(const IntervalMatrix &M) :
    value(M)
{
    assert_release(M.cols()==DimMatrix && M.rows()==DimMatrix);
    this->contractor();
}

template <unsigned int DimGroup, unsigned int DimMatrix>
void LieBaseMatrix<DimGroup,DimMatrix>::contractor() {
}

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
LieBaseMatrix<DimGroup, DimMatrix> operator*
	(const LieBaseMatrix<DimGroup,DimMatrix> &A,
         const LieBaseMatrix<DimGroup,DimMatrix> &B) {
   return LieBaseMatrix<DimGroup,DimMatrix>(A.getValue()*B.getValue());
}


template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix>& LieBaseMatrix<DimGroup,DimMatrix>::leftProd
		(const LieBaseMatrix<DimGroup,DimMatrix> &A) {
    this->value=A.getValue()*this->value;
    this->contractor();
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix>& LieBaseMatrix<DimGroup,DimMatrix>::rightProd
		(const LieBaseMatrix<DimGroup,DimMatrix> &B) {
    this->value=this->value*B.getValue();
    this->contractor();
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::inverse() 
const
{
   IntervalMatrix inv = inverse_enclosure(this->value);
   return LieBaseMatrix<DimGroup,DimMatrix>(inv);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
void LieBaseMatrix<DimGroup,DimMatrix>::inverse_inplace() 
{
   this->value = inverse_enclosure(this->value);
   this->contract();
}

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::IleftProd
		(const LieBaseMatrix &A) const
{
   IntFullPivLU pivLU(this->value);
   IntervalMatrix RA = pivLU.solve(A);
   return LieBaseMatrix<DimGroup,DimMatrix>(RA);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::IrightProd
		(const LieBaseMatrix &A) const
{
   IntFullPivLU pivLU(this->value.transpose());
   IntervalMatrix RA = pivLU.solve(A.transpose());
   return LieBaseMatrix<DimGroup,DimMatrix>(RA.transpose());
}

template <unsigned int DimGroup, unsigned int DimMatrix>
LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::center()
{
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
     os << "Lie:" << x.value;
     return os;
}

}
