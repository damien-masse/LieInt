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

         void contractor();

         const IntervalMatrix getValue() const;

	 static IntervalMatrix
		representationAlgebra(const IntervalVector &V);
	 static Matrix
		representationAlgebra(const Vector &V);

         template <unsigned int DG, unsigned int DM>
         friend LieBaseMatrix<DG,DM> 
		operator*(const LieBaseMatrix<DG,DM> &A,
			  const LieBaseMatrix<DG,DM> &B);
         LieBaseMatrix& leftProd(const LieBaseMatrix &A, LieBaseMatrix &B);
         LieBaseMatrix& rightProd(LieBaseMatrix &A, const LieBaseMatrix &B);

         template <unsigned int DG, unsigned int DM>
	 friend std::ostream& operator<<(std::ostream& os,
				 const LieBaseMatrix<DG,DM>& x);
    };

    /* advanced class for matrix Lie group */
    template <unsigned int DimGroup, unsigned int DimMatrix, class Base>
    class LieExtMatrix {
        public:
	   LieExtMatrix();
           LieExtMatrix(const Base &Mat);
           LieExtMatrix(const IntervalMatrix &M);

           void contractor();
           
           const Base& getValueBase() const;
           
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
           LieExtMatrix& leftProd(const LieExtMatrix &A, LieExtMatrix &B);
           LieExtMatrix& leftProd(const Base &A, LieExtMatrix &B);
           LieExtMatrix& rightProd(LieExtMatrix &A, const LieExtMatrix &B);
           LieExtMatrix& rightProd(LieExtMatrix &A, const Base &B);

           template <unsigned int DG, unsigned int DM, class BS>
           friend std::ostream& operator<< (std::ostream& os,
                                 const LieExtMatrix<DG,DM,BS>& x);
    };

}
