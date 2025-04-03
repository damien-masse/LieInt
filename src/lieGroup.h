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

         LieBaseMatrix(); /* Id */

         LieBaseMatrix(const IntervalMatrix &M);

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

         template <unsigned int DG, unsigned int DM>
	 friend std::ostream& operator<<(std::ostream& os,
				 const LieBaseMatrix<DG,DM>& x);

      private:
	 IntervalMatrix value;
    };

    /* advanced class for matrix Lie group */
    template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
    class LieExtMatrix {
        public:
	   LieExtMatrix();
           LieExtMatrix(const Base &Left, const Base &Cent);
           LieExtMatrix(const Base &B);

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

           template <unsigned int DG, unsigned int DM, class BS>
           friend std::ostream& operator<< (std::ostream& os,
                                 const LieExtMatrix<DG,DM,BS>& x);

	   private:
	      Base left,cent;
    };

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>::LieExtMatrix() :
        left(), cent() {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>::LieExtMatrix(const Base& Left,
			const Base& Cent) :
        left(Left), cent(Cent) {
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>::LieExtMatrix(const Base& B)
			:
        left(B), cent() {
   cent = left.center();
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
Base LieExtMatrix<DimGroup,DimMatrix,Base>::getValueBase() const {
     Base b = left*cent;
     return b;
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>
	 operator*(const LieExtMatrix<DimGroup,DimMatrix,Base> &A,
		    const LieExtMatrix<DimGroup,DimMatrix,Base> &B) {
     Base centR = A.cent * B.cent;
     Base leftR = A.left * A.cent.IrightProd(A.cent * B.left);
     return LieExtMatrix<DimGroup,DimMatrix,Base>(leftR,centR);
}


template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>
	 operator*(const LieExtMatrix<DimGroup,DimMatrix,Base> &A,
		    const Base &B) {
     LieExtMatrix<DimGroup,DimMatrix,Base> C(B);
     return C.leftProd(A);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>
	 operator*(const Base &A,
	 	const LieExtMatrix<DimGroup,DimMatrix,Base> &B) {
     LieExtMatrix<DimGroup,DimMatrix,Base> C(A);
     return C.rightProd(B);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>& 
	LieExtMatrix<DimGroup,DimMatrix,Base>::leftProd
                (const LieExtMatrix<DimGroup,DimMatrix,Base> &A) {
    this->cent.leftProd(A.cent);
    this->left = A.left * A.cent.IrightProd(A.cent * this->left);
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
LieExtMatrix<DimGroup,DimMatrix,Base>& 
	LieExtMatrix<DimGroup,DimMatrix,Base>::rightProd
                (const LieExtMatrix<DimGroup,DimMatrix,Base> &B) {
    this->left=this->left * this->cent.IrightProd(this->cent * B.left);
    this->cent.rightProd(B.cent);
    return (*this);
}

template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
std::ostream& operator<<(std::ostream& os,
                             const LieExtMatrix<DimGroup,DimMatrix,Base>& x) {
     os << "Lie:(" << x.left << "*" << x.cent << ")" ;
     return os;
}

}
