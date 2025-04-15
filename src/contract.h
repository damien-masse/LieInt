/** 
 *  \file contract.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <codac>

using namespace codac2;

namespace lieInt
{

  template <class A>
  concept VectorInterval = 
     std::is_same<typename A::Scalar,Interval>::value &&
     A::ColsAtCompileTime==1;

  template <class A>
  concept VectorRowInterval = 
     std::is_same<typename A::Scalar,Interval>::value &&
     A::IsVectorAtCompileTime==1;

  template <class A>
  concept MatrixInterval = 
     std::is_same<typename A::Scalar,Interval>::value;


  /** \brief constraint norm_2(A) = 1.0, specialised version
   * 
   *  \param A a vector of size 2
   *
   *  \return true if non empty, false if empty
   */
template <VectorInterval Derived>
requires (Derived::RowsAtCompileTime==2)
inline bool contract_unitVector2(Eigen::MatrixBase<Derived> const &CA) {
   Eigen::MatrixBase<Derived> &A = const_cast<Eigen::MatrixBase<Derived>& >(CA);
   SqrOp::bwd(1.0-SqrOp::fwd(A[1]),A[0]);
   SqrOp::bwd(1.0-SqrOp::fwd(A[0]),A[1]);
   return (!A[0].is_empty() && !A[1].is_empty());
}


  /** \brief constraint norm_2(A) = 1.0, generic version
   * 
   *  \param A a vector of size n, n>=1
   *
   *  \return true if non empty, false if empty
   */
template <VectorRowInterval Derived>
inline bool contract_unitVector(Eigen::MatrixBase<Derived> const &CA) {
   Eigen::MatrixBase<Derived> &A = const_cast<Eigen::MatrixBase<Derived>& >(CA);
   const int d = A.size();
   assert_release(d>=1);
   if (d==1) {
      A[0] &= Interval(1.0);
      return (!A[0].is_empty());
   }  
   IntervalVector sq(d);
   for (int i=0;i<d;i++) sq[i]=SqrOp::fwd(A[i]);
   if (d>2) {
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
   } else {
     AddOp::bwd(Interval(1.0),sq[0],sq[1]);
   }
   if (sq[d-2].is_empty() || sq[d-1].is_empty()) return false;
   for (int i=0;i<d;i++) SqrOp::bwd(sq[i],A[i]);
   return true;
}
  
  /** \brief constraint MMt = Id, specialised version
   * 
   *  \param A a matrix of size 2*2
   *
   *  \return true if non empty, false if empty
   */
template <MatrixInterval Derived>
requires (Derived::RowsAtCompileTime==2 && Derived::ColsAtCompileTime==2)
inline bool contract_rotMatrix2(Eigen::MatrixBase<Derived> const &CA) {
   Eigen::MatrixBase<Derived> &A = const_cast<Eigen::MatrixBase<Derived>& >(CA);
   A(0,0) &= A(1,1);
   A(1,0) &= -A(0,1);
   if (!contract_unitVector2(A.col(0))) return false;
   A(1,1) &= A(0,0);
   A(0,1) &= -A(1,0);
   return true;
}


  /** \brief constraint AB = 0
   * 
   *  \param A and B are row or col vectors
   *
   *  \return true if non empty, false if empty
   */
template <VectorRowInterval Derived, VectorRowInterval OtherDerived>
inline bool contract_orthoVectors(Eigen::MatrixBase<Derived> const &CA,
			Eigen::MatrixBase<OtherDerived> const &CB) {
   Eigen::MatrixBase<Derived> &x1 = 
		const_cast<Eigen::MatrixBase<Derived>& >(CA);
   Eigen::MatrixBase<OtherDerived> &x2 = 
			const_cast<Eigen::MatrixBase<OtherDerived>& >(CB);
   const Eigen::Index n=x1.size();
   assert (n==x2.size());
   std::vector<Interval> sums(n), prods(n);

    // Forward propagation

    for(Eigen::Index i = 0 ; i < n ; i++)
    {
      prods[i] = x1[i]*x2[i];
      sums[i] = prods[i];
      if(i > 0) sums[i] += sums[i-1];
    } 

  // Backward propagation

     sums[n-1] &= Interval(0.0);
     if (sums[n-1].is_empty()) { x1[0].set_empty(); return false; }

    for(Eigen::Index i = n-1 ; i >= 0 ; i--)
    {
       if(i > 0) AddOp::bwd(sums[i],sums[i-1],prods[i]);
       else prods[0] &= sums[0];
       MulOp::bwd(prods[i],x1[i],x2[i]);
       if (x1[i].is_empty()) return false;
    }
    return true;
}

  /** \brief constraint MMt = Id
   * 
   *  \param A a matrix of size n*n
   *
   *  \return true if non empty, false if empty
   */
template <MatrixInterval Derived>
inline bool contract_rotMatrix(Eigen::MatrixBase<Derived> const &CA) {
   Eigen::MatrixBase<Derived> &A=const_cast<Eigen::MatrixBase<Derived> &>(CA);
   const Eigen::Index n = A.rows();
   assert (n==A.cols());
   for (int i=0;i<n;i++) {
        if (!contract_unitVector(A.row(i))) return false;
        if (!contract_unitVector(A.col(i))) return false;
   }
   for (int i=0;i<n;i++) {
     for (int j=i+1;j<n;j++) {
        if (!contract_orthoVectors(A.col(i),A.col(j))) return false;
        if (!contract_orthoVectors(A.row(i),A.row(j))) return false;
     }
   }
   for (int i=0;i<n;i++) {
        if (!contract_unitVector(A.row(i))) return false;
        if (!contract_unitVector(A.col(i))) return false;
   }
   return true;
}

}
