/** 
 *  \file lieGroup.h  
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
    /* define a basic class for matrix Lie groups */
    template <unsigned int DimGroup, unsigned int DimMatrix>
    class LieBaseMatrix {
      public:

         LieBaseMatrix(); /* full */

         explicit LieBaseMatrix(const IntervalMatrix &M);
  
         static const LieBaseMatrix Empty();
         static const LieBaseMatrix Identity();

         virtual void contractor();

         const IntervalMatrix& getValue() const;

	 static IntervalMatrix
		representationAlgebra(const IntervalVector &V);
	 static Matrix
		representationAlgebra(const Vector &V);

         template <unsigned int DG, unsigned int DM>
         friend LieBaseMatrix<DG,DM> 
		operator*(const LieBaseMatrix<DG,DM> &A,
			  const LieBaseMatrix<DG,DM> &B);
         LieBaseMatrix& leftProd(const LieBaseMatrix &A);
         LieBaseMatrix& rightProd(const LieBaseMatrix &B);

         LieBaseMatrix inverse() const;
         void inverse_inplace();
         LieBaseMatrix IleftProd(const LieBaseMatrix &A) const;
         LieBaseMatrix IrightProd(const LieBaseMatrix &A) const;
         LieBaseMatrix center();

         bool is_empty() const;
         void set_empty();

         template <unsigned int DG, unsigned int DM>
	 friend std::ostream& operator<<(std::ostream& os,
				 const LieBaseMatrix<DG,DM>& x);

      protected:
	 IntervalMatrix value;
         bool empty;
    };

    /* advanced class for matrix Lie group */
    template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
    class LieExtMatrix {
        public:
	   LieExtMatrix();
           LieExtMatrix(const Base &Left, const Base &Cent);
           LieExtMatrix(const Base &B);

           static const LieExtMatrix Empty();
           static const LieExtMatrix Identity();

           void contractor();
           
           Base getValueBase() const;
           
           template <unsigned int DG, unsigned int DM, class BS>
	   friend LieExtMatrix<DG,DM,BS>
			 operator*(const LieExtMatrix<DG,DM,BS> &A,
				    const LieExtMatrix<DG,DM,BS> &B);
           template <unsigned int DG, unsigned int DM, class BS>
	   friend LieExtMatrix<DG,DM,BS> 
			operator*(const LieExtMatrix<DG,DM,BS> &A,
				  const BS &B);
           template <unsigned int DG, unsigned int DM, class BS>
	   friend LieExtMatrix<DG,DM,BS> 
			operator*(const BS &A, const LieExtMatrix<DG,DM,BS> &B);
           LieExtMatrix& leftProd(const LieExtMatrix &A);
           LieExtMatrix& leftProd(const Base &A);
           LieExtMatrix& rightProd(const LieExtMatrix &B);
           LieExtMatrix& rightProd(const Base &B);

           bool is_empty() const;
           void set_empty();

           template <unsigned int DG, unsigned int DM, class BS>
           friend std::ostream& operator<< (std::ostream& os,
                                 const LieExtMatrix<DG,DM,BS>& x);

	   protected:
	      Base left,cent;
	      bool empty;
    };

/*** inline functions for LieBaseMatrix ***/
template <unsigned int DimGroup, unsigned int DimMatrix>
inline LieBaseMatrix<DimGroup,DimMatrix>::LieBaseMatrix() :
    value(DimMatrix,DimMatrix), empty(false)
{
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline LieBaseMatrix<DimGroup,DimMatrix>::LieBaseMatrix(const IntervalMatrix &M) :
    value(M), empty(M.is_empty())
{
    assert_release(M.cols()==DimMatrix && M.rows()==DimMatrix);
    if (!this->empty) this->contractor();
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline const LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::Empty() {
    return LieBaseMatrix<DimGroup,DimMatrix>(IntervalMatrix::Constant
                (DimMatrix,DimMatrix,Interval::empty()));
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline const LieBaseMatrix<DimGroup,DimMatrix> LieBaseMatrix<DimGroup,DimMatrix>::Identity() {
    return LieBaseMatrix<DimGroup,DimMatrix>(IntervalMatrix::Identity
                (DimMatrix,DimMatrix));
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline void LieBaseMatrix<DimGroup,DimMatrix>::contractor() {
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline const IntervalMatrix& LieBaseMatrix<DimGroup,DimMatrix>::getValue() const {
    return this->value;
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline LieBaseMatrix<DimGroup, DimMatrix> operator*
        (const LieBaseMatrix<DimGroup,DimMatrix> &A,
         const LieBaseMatrix<DimGroup,DimMatrix> &B) {
   if (A.empty || B.empty) return LieBaseMatrix<DimGroup, DimMatrix>::Empty();
   return LieBaseMatrix<DimGroup,DimMatrix>(A.getValue()*B.getValue());
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline LieBaseMatrix<DimGroup,DimMatrix>& LieBaseMatrix<DimGroup,DimMatrix>::leftProd
                (const LieBaseMatrix<DimGroup,DimMatrix> &A) {
    if (this->empty) return (*this);
    if (A.empty) {  this->set_empty(); return (*this); }
    this->value=A.getValue()*this->value;
    this->contractor();
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline LieBaseMatrix<DimGroup,DimMatrix>& LieBaseMatrix<DimGroup,DimMatrix>::rightProd
                (const LieBaseMatrix<DimGroup,DimMatrix> &A) {
    if (this->empty) return (*this);
    if (A.empty) {  this->set_empty(); return (*this); }
    this->value=this->value*A.getValue();
    this->contractor();
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline bool LieBaseMatrix<DimGroup,DimMatrix>::is_empty() const {
   return this->empty;
}

template <unsigned int DimGroup, unsigned int DimMatrix>
inline void LieBaseMatrix<DimGroup,DimMatrix>::set_empty() {
   this->empty=true;
}

/*** inline functions for LieExtMatrix ***/
template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>::LieExtMatrix() :
        left(), cent(LieExtMatrix<DimGroup,DimMatrix,Base>::Identity()),
        empty(false) {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>::LieExtMatrix(const Base& Left,
			const Base& Cent) :
        left(Left), cent(Cent), empty(left.is_empty() || cent.is_empty()) {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>::LieExtMatrix(const Base& B)
			:
        left(B), cent(), empty(left.is_empty()) {
   if (this->empty) { cent=Base::Empty(); return; }
   cent = left.center();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline const LieExtMatrix<DimGroup,DimMatrix,Base> LieExtMatrix<DimGroup,DimMatrix,Base>::Empty() {
   return LieExtMatrix(Base::Empty(), Base::Empty());
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline const LieExtMatrix<DimGroup,DimMatrix,Base> LieExtMatrix<DimGroup,DimMatrix,Base>::Identity() {
   return LieExtMatrix(Base::Identity(), Base::Identity());
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline Base LieExtMatrix<DimGroup,DimMatrix,Base>::getValueBase() const {
     if (this->is_empty()) return Base::Empty();
     Base b = left*cent;
     return b;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>
	 operator*(const LieExtMatrix<DimGroup,DimMatrix,Base> &A,
		    const LieExtMatrix<DimGroup,DimMatrix,Base> &B) {
     if (A.is_empty() || B.is_empty()) 
		return LieExtMatrix<DimGroup,DimMatrix,Base>::Empty();
     Base centR = A.cent * B.cent;
     Base leftR = A.left * A.cent.IrightProd(A.cent * B.left);
     return LieExtMatrix<DimGroup,DimMatrix,Base>(leftR,centR);
}


template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>
	 operator*(const LieExtMatrix<DimGroup,DimMatrix,Base> &A,
		    const Base &B) {
     if (A.is_empty() || B.is_empty()) 
		return LieExtMatrix<DimGroup,DimMatrix,Base>::Empty();
     LieExtMatrix<DimGroup,DimMatrix,Base> C(B);
     return C.leftProd(A);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>
	 operator*(const Base &A,
	 	const LieExtMatrix<DimGroup,DimMatrix,Base> &B) {
     if (A.is_empty() || B.is_empty()) 
		return LieExtMatrix<DimGroup,DimMatrix,Base>::Empty();
     LieExtMatrix<DimGroup,DimMatrix,Base> C(A);
     return C.rightProd(B);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>& 
	LieExtMatrix<DimGroup,DimMatrix,Base>::leftProd
                (const LieExtMatrix<DimGroup,DimMatrix,Base> &A) {
     if (this->is_empty()) return (*this);
     if (A.is_empty()) { this->set_empty(); return (*this); } 
    this->cent.leftProd(A.cent);
    this->left = A.left * A.cent.IrightProd(A.cent * this->left);
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline LieExtMatrix<DimGroup,DimMatrix,Base>& 
	LieExtMatrix<DimGroup,DimMatrix,Base>::rightProd
                (const LieExtMatrix<DimGroup,DimMatrix,Base> &B) {
    if (this->is_empty()) return (*this);
    if (B.is_empty()) { this->set_empty(); return (*this); } 
    this->left=this->left * this->cent.IrightProd(this->cent * B.left);
    this->cent.rightProd(B.cent);
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline bool LieExtMatrix<DimGroup,DimMatrix,Base>::is_empty() const {
    return this->empty;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline void LieExtMatrix<DimGroup,DimMatrix,Base>::set_empty() {
    this->empty=true;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
inline std::ostream& operator<<(std::ostream& os,
                             const LieExtMatrix<DimGroup,DimMatrix,Base>& x) {
     if (x.is_empty()) { os << "Lie:(empty)"; return os; }
     os << "Lie:(" << x.left << "*" << x.cent << ")" ;
     return os;
}

}
