/** 
 *  \file lieSE2.h  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <codac>
#include "lieGroup.h"

using namespace codac2;

namespace lieInt
{
    class SE2Base : public LieBaseMatrix<3,3,SE2Base> {
      public :

        SE2Base(); /* full */

        explicit SE2Base(const LieIntervalMatrix &M);
        SE2Base(const Interval &X, const Interval &Y, const Interval &theta);

        void contractor();

        static LieIntervalMatrix
                representationAlgebra(const IntervalVector &V);
        static LieMatrix
                representationAlgebraV(const Vector &V);


        SE2Base inverse() const;
        void inverse_inplace();
        SE2Base IleftProd(const SE2Base &A) const;
        SE2Base IrightProd(const SE2Base &A) const;
        SE2Base center(bool uncertLeft=true);

        friend std::ostream& operator<<(std::ostream& os,
                                 const SE2Base& x);

        void draw(bool drawDirection=true,
                const StyleProperties& sBox = StyleProperties(),
                const StyleProperties& sPie = StyleProperties()) const;
        void draw(Figure2D &fig2D, bool drawDirection=true,
		const StyleProperties& sBox = StyleProperties(),
		const StyleProperties& sPie = StyleProperties()) const;

        
        static const std::vector<SE2Base::LieMatrix> LieAlgebraGenerators;

        const std::vector<Interval> angles() const;
        const Interval &getX() const;
        const Interval &getY() const;


    };
 
    class SE2Ext : public LieExtMatrix<3,3,SE2Base,SE2Ext> {
      public:
        using Base = SE2Base;
 
        using LieExtMatrix<3,3,SE2Base,SE2Ext>::LieExtMatrix;
    
        void draw(Figure2D &fig2D, bool drawDirection=false,
		const StyleProperties& sBox = StyleProperties(),
		const StyleProperties& sPie = StyleProperties()) const;
        void draw(bool drawDirection=false,
		const StyleProperties& sBox = StyleProperties(),
		const StyleProperties& sPie = StyleProperties()) const;

    };

/* SE2Base inline */

inline SE2Base::SE2Base() : 
	LieBaseMatrix<3,3,SE2Base>
			  ({ { {-1,1}, {-1,1}, {-oo,oo} },
                             { {-1,1}, {-1,1}, {-oo,oo} },
                             { 0     , 0     , 1        } })
{ }
inline SE2Base::SE2Base(const SE2Base::LieIntervalMatrix &M) : 
	LieBaseMatrix<3,3,SE2Base>(M)  {
   this->contractor();
}

inline SE2Base::SE2Base(const Interval &X, const Interval &Y,
				 const Interval &theta) :
	SE2Base() {
   if (X.is_empty() || Y.is_empty() || theta.is_empty()) {
     this->set_empty(); return;
   }
   this->value(0,2)=X;
   this->value(1,2)=Y;
   this->value(0,0) = this->value(1,1) = CosOp::fwd(theta);
   this->value(1,0) = SinOp::fwd(theta);
   this->value(0,1) = -this->value(1,0);
}

inline const Interval& SE2Base::getX() const {
   return this->value(0,2);
}

inline const Interval& SE2Base::getY() const {
   return this->value(1,2);
}

}
