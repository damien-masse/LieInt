/** 
 *  \file contract.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <codac>
#include "contract.h"

using namespace codac2;

namespace lieInt
{

bool contract_unitVector2(Eigen::Vector<Interval,2> const &CA) {
   Eigen::Vector<Interval,2> &A = const_cast<Eigen::Vector<Interval,2>& >(CA);
   SqrOp::bwd(1.0-SqrOp::fwd(A[1]),A[0]);
   SqrOp::bwd(1.0-SqrOp::fwd(A[0]),A[1]);
   return (!A[0].is_empty() && !A[1].is_empty());
}

bool contract_unitVector(IntervalVector const &CA) {
   IntervalVector &A = const_cast<IntervalVector& >(CA);
   const int d = A.size();
   IntervalVector sq(d);
   for (int i=0;i<d;i++) sq[i]=SqrOp::fwd(A[i]);
   IntervalVector partial(d-2);
   partial[d-3] = AddOp::fwd(sq[d-2],sq[d-1]);
   for (int j=d-3;j>=1;j--) {
      partial[j-1] = AddOp::fwd(partial[j],sq[j]);
   }
   AddOp::bwd(Interval(1.0),partial[0],sq[0]);
   if (sq[0].is_empty()) return false;
   for (int j=1;j<=d-3;j++) {
       AddOp::bwd(partial[j-1],partial[j],sq[j]);
       if (sq[j].is_empty()) return false;
   }
   AddOp::bwd(partial[d-3],sq[d-2],sq[d-1]);
   if (sq[d-2].is_empty() || sq[d-1].is_empty()) return false;
   for (int i=0;i<d;i++) SqrOp::bwd(sq[i],A[i]);
   return true;
}

  /** \brief constraint MMt = Id, specialised version
   * 
   *  \param A a matrix of size 2*2
   *
   *  \return nothing, modifies A
   */
bool contract_rotMatrix2(Eigen::Matrix<Interval,2,2>& A) {
   A(0,0) &= A(1,1);
   A(1,0) &= -A(0,1);
   if (!contract_unitVector2(A.topLeftCorner<1,2>())) return false;
   A(1,1) &= A(0,0);
   A(0,1) &= -A(1,0);
   return true;

}

  /** \brief constraint MMt = Id
   * 
   *  \param A a matrix of size n*n
   *
   *  \return nothing, modifies A
   */
  bool contract_rotMatrix(IntervalMatrix& A) {
     const int d = A.rows();
     assert_release(A.cols()==d);
     for (int i=0;i<d;i++) {
        if (contract_unitVector(A.row(i))) return false;
        if (contract_unitVector(A.col(i))) return false;
     }
/*     for (int i=0;i<d;i++) {
       for (int j=0;j<d;j++) {
          if (i==j) continue;
          MulOp::bwd(0.0,A.col(i).transpose(),A.col(j));
		// DOES NOT COMPILE FOR NOW 
          if (A.col(i).is_empty()) return false;
       }
     }
*/
     for (int i=0;i<d;i++) {
        if (contract_unitVector(A.row(i))) return false;
        if (contract_unitVector(A.col(i))) return false;
     }
     return true;
  } 
}
