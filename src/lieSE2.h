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
    class SE2Base : public LieBaseMatrix<3,3> {
      public :

        SE2Base(); /* full */

        explicit SE2Base(const IntervalMatrix &M);

        static const SE2Base Empty();
        static const SE2Base Identity();

        void contractor();

        static IntervalMatrix
                representationAlgebra(const IntervalVector &V);
        static Matrix
                representationAlgebra(const Vector &V);

        SE2Base inverse() const;
        void inverse_inplace();
        SE2Base IleftProd(const SE2Base &A) const;
        SE2Base IrightProd(const SE2Base &A) const;
        SE2Base center();

        friend std::ostream& operator<<(std::ostream& os,
                                 const SE2Base& x);

        void draw(bool drawDirection=true,
                const StyleProperties& sBox = StyleProperties(),
                const StyleProperties& sPie = StyleProperties());
        void draw(Figure2DInterface &fig2D, bool drawDirection=true,
		const StyleProperties& sBox = StyleProperties(),
		const StyleProperties& sPie = StyleProperties());

        
        static const std::vector<Matrix> LieAlgebraGenerators;

      private:
	explicit SE2Base(const LieBaseMatrix<3,3> &A);
        
        std::vector<Interval> angles() const;

    };
 
    class SE2Ext : public LieExtMatrix<3,3,SE2Base> {
    
        void draw(Figure2DInterface &fig2D, bool drawDirection=true,
		const StyleProperties& s = StyleProperties());

    };

/* SE2Base inline */

inline SE2Base::SE2Base() : 
	LieBaseMatrix<3,3>({ { {-1,1}, {-1,1}, {-oo,oo} },
                             { {-1,1}, {-1,1}, {-oo,oo} },
                             { 0     , 0     , 1        } })
{ }
inline SE2Base::SE2Base(const IntervalMatrix &M) : LieBaseMatrix<3,3>(M)  {
   this->contractor();
}
inline SE2Base::SE2Base(const LieBaseMatrix<3,3> &A) :
	LieBaseMatrix<3,3>(A) { }
inline const SE2Base SE2Base::Empty() {
        return SE2Base(LieBaseMatrix<3,3>::Empty());
}
inline const SE2Base SE2Base::Identity() {
      return SE2Base(LieBaseMatrix<3,3>::Empty());
}

}
