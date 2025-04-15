/** 
 *  \file lieSO3.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */


#include <codac>
#include <cmath>
#include "lieSO3.h"
#include "contract.h"

using namespace codac2;

namespace lieInt
{

const std::vector<SO3Base::LieMatrix> SO3Base::LieAlgebraGenerators
     { { {0,-1,0}, {1,0,0}, {0,0,0} },
       { {0,0,1},  {0,0,0}, {-1,0,0} },
       { {0,0,0},  {0,0,-1}, {0,1,0} } };


SO3Base::SO3Base(const SO3Base::CardanAngles &ang) :
	SO3Base() {
   if (ang.roll.is_empty() || ang.pitch.is_empty() || ang.yaw.is_empty()) {
     this->set_empty(); return;
   }
  /* (cos psi cos theta, cos psi sin theta sin phi - sin psi cos phi, sin psi sin phi + cos psi cos phi sin theta )
  (sin psi cos theta, cos psi cos phi - sin psi sin theta sin phi, cos phi sin theta sin phi - cos psi sin phi )                                    )
              (  -sin theta     , cos theta sin phi, cos theta cos phi) */
   Interval cPhi = CosOp::fwd(ang.roll);
   Interval sPhi = SinOp::fwd(ang.roll);
   Interval cPsi = CosOp::fwd(ang.yaw);
   Interval sPsi = SinOp::fwd(ang.yaw);
   Interval cTheta = CosOp::fwd(ang.pitch);
   Interval sTheta = SinOp::fwd(ang.pitch);
   this->value(0,0) = cPsi * cTheta;
   this->value(1,0) = sPsi * cTheta;
   this->value(2,0) = - sTheta;
   this->value(0,1) = cPsi * sTheta * sPhi - sPsi * cPhi;
   this->value(1,1) = cPsi * cPhi + sPsi * sPhi * sTheta;
   this->value(2,1) = cTheta * sPhi;
   this->value(0,2) = sPsi * sPhi + cPsi * cPhi * sTheta;
   this->value(1,2) = - cPsi * sPhi + sPsi * cPhi * sTheta;
   this->value(2,2) = cTheta * cPhi;
   this->contractor();
}

SO3Base::LieIntervalMatrix SO3Base::representationAlgebra(const IntervalVector &V) {
   assert(V.size()==SO3Base::LieAlgebraGenerators.size());
   SO3Base::LieIntervalMatrix R = SO3Base::LieIntervalMatrix::Zero();
   for (int i=0;i<V.size();i++) {
       R += V[i]*SO3Base::LieAlgebraGenerators[i];
   }
   return R;
}

SO3Base::LieMatrix SO3Base::representationAlgebraV(const Vector &V) {
   assert(V.size()==SO3Base::LieAlgebraGenerators.size());
   SO3Base::LieMatrix R = SO3Base::LieMatrix::Zero();
   for (int i=0;i<V.size();i++) {
       R += V[i]*SO3Base::LieAlgebraGenerators[i];
   }
   return R;
}


SO3Base SO3Base::IleftProd(const SO3Base &A) const {
     if (this->empty || A.is_empty()) return SO3Base::Empty();
     SO3Base R = this->inverse();
     R.rightProd(A);
     return R;
}

SO3Base SO3Base::IrightProd(const SO3Base &A) const {
     if (this->empty || A.is_empty()) return SO3Base::Empty();
     SO3Base R = this->inverse();
     R.leftProd(A);
     return R;
}

SO3Base SO3Base::center() {
     if (this->empty) return SO3Base::Empty();
     Matrix Ct = this->value.mid();
     if ((Ct(2,1)==0.0 && Ct(2,2)==0.0) || (Ct(1,0)==0.0 && Ct(0,0)==0.0)) {
          return SO3Base::Identity();
     }
     double phi = std::atan2(Ct(2,1),Ct(2,2));
     double psi = std::atan2(Ct(1,0),Ct(0,0));
     double theta = -std::asin(Ct(2,0)/
		std::sqrt(Ct(0,0)*Ct(0,0)+Ct(1,0)*Ct(1,0)+Ct(2,0)*Ct(2,0)));
     SO3Base::CardanAngles crd { Interval(psi), Interval(theta), Interval(phi) };
     SO3Base cent(crd);
     SO3Base::LieIntervalMatrix &M = this->value;
     M = M * cent.getValue().transpose();
     return cent;
}

const SO3Base::CardanAngles SO3Base::anglesCardan() const {
     if (this->empty) 
	return { Interval::empty(), Interval::empty(), Interval::empty() };
     const SO3Base::LieIntervalMatrix &Ct = this->value;
     SO3Base::CardanAngles ret;
     ret.roll = Atan2Op::fwd(Ct(2,1),Ct(2,2));
     ret.yaw = Atan2Op::fwd(Ct(1,0),Ct(0,0));
     ret.pitch = -AsinOp::fwd(Ct(2,0));
     return ret;
}


std::ostream& operator<<(std::ostream& os, const SO3Base& x) {
     if (x.is_empty()) { os << "SO3:empty"; return os; }
     const SO3Base::LieIntervalMatrix &M = x.getValue();
     os << "SO3:" << M ;
     return os;
}

void SO3Base::draw3D(Figure3D &fig3D, const Vector& pos,
		double size, const StyleProperties& s) const {
     if (this->empty) return;
     SO3Base save(*this);
     SO3Base cent = save.center(); /* center modify save ... */
     fig3D.draw_hull(pos,cent.getValue().mid(),size,s);
     fig3D.draw_box(pos+size*this->value.col(0),s);
}

void SO3Ext::draw3D(Figure3D &fig3D, const Vector& pos,
                double size, const StyleProperties& s) const {
     if (this->empty) return;
     fig3D.draw_uncertain_hull(pos,
		this->cent.getValue(),this->left.getValue(),size,s);
}


}
