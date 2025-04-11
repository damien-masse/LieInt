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
    template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
    class LieBaseMatrix {
      public:
         
         using LieIntervalMatrix = Eigen::Matrix<Interval,DimMatrix,DimMatrix>;
         using LieMatrix = Eigen::Matrix<double,DimMatrix,DimMatrix>;

         LieBaseMatrix(); /* full */

         explicit LieBaseMatrix(const LieIntervalMatrix &M);
  
         static const Derived Empty();
         static const Derived Identity();

         Derived &derived();

         virtual void contractor();

         const LieIntervalMatrix& getValue() const;

	 static LieIntervalMatrix
		representationAlgebra(const IntervalVector &V);
	 static LieMatrix
		representationAlgebraV(const Vector &V);

         template <unsigned int DG, unsigned int DM, class Der>
         friend Der operator*(const LieBaseMatrix<DG,DM,Der> &A,
		 	      const LieBaseMatrix<DG,DM,Der> &B);
         Derived& leftProd(const LieBaseMatrix &A);
         Derived& rightProd(const LieBaseMatrix &B);

         Derived inverse() const;
         void inverse_inplace();
         Derived IleftProd(const LieBaseMatrix &A) const;
         Derived IrightProd(const LieBaseMatrix &A) const;
         Derived center();

         bool is_empty() const;
         void set_empty();

         template <unsigned int DG, unsigned int DM, class Der>
	 friend std::ostream& operator<<(std::ostream& os,
				 const LieBaseMatrix<DG,DM,Der>& x);

      protected:
	 LieIntervalMatrix value;
         bool empty;
    };

    /* advanced class for matrix Lie group */
    template <unsigned int DimGroup, unsigned int DimMatrix, class Base,
	class Derived>
    class LieExtMatrix {
        public:
           using BaseLieGroup = Base;
           
	   LieExtMatrix();
           LieExtMatrix(const Base &Left, const Base &Cent);
           LieExtMatrix(const Base &B);
           LieExtMatrix(const IntervalMatrix &B);

           static const Derived Empty();
           static const Derived Identity();

           Derived &derived();

           void contractor();
           
           Base getValueBase() const;
           Base getLeft() const;
           Base getCenter() const;
           
           template <unsigned int DG, unsigned int DM, class BS, class DR>
	   friend DR operator*(const LieExtMatrix<DG,DM,BS,DR> &A,
				    const LieExtMatrix<DG,DM,BS,DR> &B);
           template <unsigned int DG, unsigned int DM, class BS, class DR>
	   friend DR operator*(const LieExtMatrix<DG,DM,BS,DR> &A,
				  const BS &B);
           template <unsigned int DG, unsigned int DM, class BS, class DR>
	   friend DR operator*(const BS &A, 
				    const LieExtMatrix<DG,DM,BS,DR> &B);
           Derived& leftProd(const LieExtMatrix &A);
           Derived& leftProd(const Base &A);
           Derived& rightProd(const LieExtMatrix &B);
           Derived& rightProd(const Base &B);

           bool is_empty() const;
           void set_empty();

           template <unsigned int DG, unsigned int DM, class BS, class DR>
           friend std::ostream& operator<< (std::ostream& os,
                                 const LieExtMatrix<DG,DM,BS,DR>& x);

	   protected:
	      Base left,cent;
	      bool empty;
    };

/*** inline functions for LieBaseMatrix ***/
template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline LieBaseMatrix<DimGroup,DimMatrix,Derived>::LieBaseMatrix() :
    value(DimMatrix,DimMatrix), empty(false)
{
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline LieBaseMatrix<DimGroup,DimMatrix,Derived>::LieBaseMatrix(const LieBaseMatrix<DimGroup,DimMatrix,Derived>::LieIntervalMatrix &M) :
    value(M), empty(M.is_empty())
{
    if (!this->empty) this->contractor();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline const Derived LieBaseMatrix<DimGroup,DimMatrix,Derived>::Empty() {
    return LieBaseMatrix<DimGroup,DimMatrix,Derived>(IntervalMatrix::Constant
                (DimMatrix,DimMatrix,Interval::empty())).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline const Derived LieBaseMatrix<DimGroup,DimMatrix,Derived>::Identity() {
    return LieBaseMatrix<DimGroup,DimMatrix,Derived>(LieBaseMatrix<DimGroup,DimMatrix,Derived>::LieIntervalMatrix::Identity()).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline void LieBaseMatrix<DimGroup,DimMatrix,Derived>::contractor() {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline const LieBaseMatrix<DimGroup,DimMatrix,Derived>::LieIntervalMatrix& LieBaseMatrix<DimGroup,DimMatrix,Derived>::getValue() const {
    return this->value;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline Derived &LieBaseMatrix<DimGroup,DimMatrix,Derived>::derived() {
   return *static_cast<Derived*>(this);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline Derived operator*
        (const LieBaseMatrix<DimGroup,DimMatrix,Derived> &A,
         const LieBaseMatrix<DimGroup,DimMatrix,Derived> &B) {
   if (A.empty || B.empty) 
	return LieBaseMatrix<DimGroup, DimMatrix,Derived>::Empty();
   return 
	LieBaseMatrix<DimGroup,DimMatrix,Derived>(A.getValue()*B.getValue()).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline Derived& LieBaseMatrix<DimGroup,DimMatrix,Derived>::leftProd
                (const LieBaseMatrix<DimGroup,DimMatrix,Derived> &A) {
    if (this->empty) return (*this).derived();
    if (A.empty) {  this->set_empty(); return (*this).derived(); }
    this->value=A.getValue()*this->value;
    this->contractor();
    return (*this).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline Derived& LieBaseMatrix<DimGroup,DimMatrix,Derived>::rightProd
                (const LieBaseMatrix<DimGroup,DimMatrix,Derived> &A) {
    if (this->empty) return (*this).derived();
    if (A.empty) {  this->set_empty(); return (*this).derived(); }
    this->value=this->value*A.getValue();
    this->contractor();
    return (*this).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Derived>
inline bool LieBaseMatrix<DimGroup,DimMatrix,Derived>::is_empty() const {
   return this->empty;
}

template <unsigned int DimGroup, unsigned int DimMatrix,class Derived>
inline void LieBaseMatrix<DimGroup,DimMatrix,Derived>::set_empty() {
   this->empty=true;
}

/*** inline functions for LieExtMatrix ***/
template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::LieExtMatrix() :
        left(), cent(Base::Identity()),
        empty(false) {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::LieExtMatrix(const Base& Left,
			const Base& Cent) :
        left(Left), cent(Cent), empty(left.is_empty() || cent.is_empty()) {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::LieExtMatrix(const Base& B)
			:
        left(B), cent(), empty(left.is_empty()) {
   if (this->empty) { cent=Base::Empty(); return; }
   cent = left.center();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::LieExtMatrix(const IntervalMatrix& B)
			:
        LieExtMatrix(Base(B)) {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline const Derived LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::Empty() {
   return LieExtMatrix(Base::Empty(), Base::Empty()).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline const Derived LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::Identity() {
   return LieExtMatrix(Base::Identity(), Base::Identity()).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Derived &LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::derived() {
   return *static_cast<Derived*>(this);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Base LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::getValueBase() const {
     if (this->is_empty()) return Base::Empty();
     Base b = left*cent;
     return b;
}
template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Base LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::getLeft() const {
     if (this->is_empty()) return Base::Empty();
     return left;
}
template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Base LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::getCenter() const {
     if (this->is_empty()) return Base::Empty();
     return cent;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Derived operator*(const LieExtMatrix<DimGroup,DimMatrix,Base,Derived> &A,
		    const LieExtMatrix<DimGroup,DimMatrix,Base,Derived> &B) {
     if (A.is_empty() || B.is_empty()) 
		return LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::Empty();
     Base centR = A.cent * B.cent;
     Base leftR = A.left * A.cent.IrightProd(A.cent * B.left);
     return LieExtMatrix<DimGroup,DimMatrix,Base,Derived>(leftR,centR).derived();
}


template <unsigned int DimGroup, unsigned int DimMatrix, class Base,class Derived>
inline Derived operator*(const LieExtMatrix<DimGroup,DimMatrix,Base,Derived> &A,
		    const Base &B) {
     if (A.is_empty() || B.is_empty()) 
		return LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::Empty();
     LieExtMatrix<DimGroup,DimMatrix,Base,Derived> C(B);
     return C.leftProd(A);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Derived operator*(const Base &A,
	 	const LieExtMatrix<DimGroup,DimMatrix,Base,Derived> &B) {
     if (A.is_empty() || B.is_empty()) 
		return LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::Empty();
     LieExtMatrix<DimGroup,DimMatrix,Base,Derived> C(A);
     return C.rightProd(B);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Derived& 
	LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::leftProd
                (const LieExtMatrix<DimGroup,DimMatrix,Base,Derived> &A) {
     if (this->is_empty()) return (*this).derived();
     if (A.is_empty()) { this->set_empty(); return (*this).derived(); } 
    this->cent.leftProd(A.cent);
    this->left = A.left * A.cent.IrightProd(A.cent * this->left);
    return (*this).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline Derived&
	LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::rightProd
                (const LieExtMatrix<DimGroup,DimMatrix,Base,Derived> &B) {
    if (this->is_empty()) return (*this).derived();
    if (B.is_empty()) { this->set_empty(); return (*this).derived(); } 
    this->left=this->left * this->cent.IrightProd(this->cent * B.left);
    this->cent.rightProd(B.cent);
    return (*this).derived();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline bool LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::is_empty() const {
    return this->empty;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base,class Derived>
inline void LieExtMatrix<DimGroup,DimMatrix,Base,Derived>::set_empty() {
    this->empty=true;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base, class Derived>
inline std::ostream& operator<<(std::ostream& os,
                             const LieExtMatrix<DimGroup,DimMatrix,Base,Derived>& x) {
     if (x.is_empty()) { os << "Lie:(empty)"; return os; }
     os << "Lie:(" << x.left << "*" << x.cent << ")" ;
     return os;
}

}
