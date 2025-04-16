/** 
 *  \file lieSO3.h  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <codac>
#include "lieGroup.h"
#include "contract.h"

using namespace codac2;

namespace lieInt
{
    class SO3Base : public LieBaseMatrix<3,3,SO3Base> {
      public :

        struct CardanAngles {
          Interval yaw, pitch, roll; 
        };

        SO3Base(); /* full */

        explicit SO3Base(const LieIntervalMatrix &M);
        /* from "Euler" angles : Z_psi Y_theta X_phi 
           psi = yaw    theta = pitch    phi = roll
  (cos psi cos theta, cos psi sin theta sin phi - sin psi cos phi, sin psi sin phi + cos psi cos phi sin theta )
  (sin psi cos theta, cos psi cos phi - sin psi sin theta sin phi, cos phi sin theta sin phi - cos psi sin phi )                                    )
              (  -sin theta     , cos theta sin phi, cos theta cos phi) */
        SO3Base(const CardanAngles &ang);

        void contractor();

        static LieIntervalMatrix
                representationAlgebra(const IntervalVector &V);
        static LieMatrix
                representationAlgebraV(const Vector &V);


        SO3Base inverse() const;
        void inverse_inplace();
        SO3Base IleftProd(const SO3Base &A) const;
        SO3Base IrightProd(const SO3Base &A) const;
        SO3Base center();

        friend std::ostream& operator<<(std::ostream& os,
                                 const SO3Base& x);

        void draw3D(Figure3D &fig3D, const Vector &pos, 
                double size, 
                const StyleProperties& s1 = { Color::yellow(0.8) },
                const StyleProperties& s2 = { Color::dark_green(0.8) },
                const StyleProperties& s3 = { Color::yellow(0.3) }) const;

        static const std::vector<SO3Base::LieMatrix> LieAlgebraGenerators;

        const CardanAngles anglesCardan() const;
    };
 
    class SO3Ext : public LieExtMatrix<3,3,SO3Base,SO3Ext> {
      public:
        using Base = SO3Base;
 
        using LieExtMatrix<3,3,SO3Base,SO3Ext>::LieExtMatrix;
    
        void draw3D(Figure3D &fig3D, const Vector &pos,
                double size, 
                const StyleProperties& s1 = { Color::yellow(0.8) },
                const StyleProperties& s2 = { Color::dark_green(0.8) },
                const StyleProperties& s3 = { Color::yellow(0.3) }) const;

    };

/* SO3Base inline */

inline SO3Base::SO3Base() : 
    LieBaseMatrix<3,3,SO3Base>
	(SO3Base::LieIntervalMatrix::Constant(Interval(-1,1)))
{ }
inline SO3Base::SO3Base(const SO3Base::LieIntervalMatrix &M) : 
	LieBaseMatrix<3,3,SO3Base>(M) {
   this->contractor();
}


inline SO3Base SO3Base::inverse() const {
   if (this->empty) return SO3Base::Empty();
   return LieBaseMatrix<3,3,SO3Base>(this->value.transpose()).derived();
			/* we do not call a contractor */
}

inline void SO3Base::inverse_inplace() {
   if (this->empty) return;
   this->value = this->value.transpose();
}


inline void SO3Base::contractor() {
   if (this->empty) return;
   SO3Base::LieIntervalMatrix &M = this->value;
   if (!contract_rotMatrix(M)) this->set_empty();
}



}
