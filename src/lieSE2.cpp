/** 
 *  \file lieSE2.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */


#include <codac>
#include "lieSE2.h"

using namespace codac2;

namespace lieInt
{

const std::vector<Matrix> SE2Base::LieAlgebraGenerators
     { { {0,-1,0}, {1,0,0}, {0,0,0} },
       { {0,0,1},  {0,0,0}, {0,0,0} },
       { {0,0,0},  {0,0,1}, {0,0,0} } };

void SE2Base::contractor() {
   if (this->empty) return;
   IntervalMatrix &M = this->value;
   /* contract last row */
   M(2,0) &= 0.0;
   M(2,1) &= 0.0;
   M(2,2) &= 1.0;
   /* contract rotation */
   M(0,0) &= M(1,1);
   M(1,0) &= -M(0,1);
   SqrOp::bwd(1.0-SqrOp::fwd(M(1,0)),M(0,0));
   SqrOp::bwd(1.0-SqrOp::fwd(M(0,0)),M(1,0));
   M(1,1) = M(0,0);
   M(0,1) = -M(1,0);
   /* emptiness */
   if (M.is_empty()) M.set_empty(); 
}

IntervalMatrix SE2Base::representationAlgebra(const IntervalVector &V) {
   assert(V.size()==SE2Base::LieAlgebraGenerators.size());
   IntervalMatrix R = IntervalMatrix::Zero(3,3);
   for (int i=0;i<V.size();i++) {
       R += V[i]*SE2Base::LieAlgebraGenerators[i];
   }
   return R;
}

Matrix SE2Base::representationAlgebra(const Vector &V) {
   assert(V.size()==SE2Base::LieAlgebraGenerators.size());
   Matrix R = Matrix::Zero(3,3);
   for (int i=0;i<V.size();i++) {
       R += V[i]*SE2Base::LieAlgebraGenerators[i];
   }
   return R;
}

SE2Base SE2Base::inverse() const {
   if (this->empty) return SE2Base::Empty();
   IntervalMatrix M = this->value;
   M(0,1) = -M(0,1);
   M(1,0) = -M(1,0);
   M.block<2,1>(0,2) = -M.topLeftCorner<2,2>()*M.block<2,1>(0,2);
   return SE2Base(M);
}

void SE2Base::inverse_inplace() {
   if (this->empty) return;
   IntervalMatrix &M = this->value;
   M(0,1) = -M(0,1);
   M(1,0) = -M(1,0);
   M.block<2,1>(0,2) = -M.topLeftCorner<2,2>()*M.block<2,1>(0,2);
}

SE2Base SE2Base::IleftProd(const SE2Base &A) const {
     if (this->empty || A.is_empty()) return SE2Base::Empty();
     SE2Base R = this->inverse();
     R.rightProd(A);
     return R;
}

SE2Base SE2Base::IrightProd(const SE2Base &A) const {
     if (this->empty || A.is_empty()) return SE2Base::Empty();
     SE2Base R = this->inverse();
     R.leftProd(A);
     return R;
}

std::vector<Interval> SE2Base::angles() const {
     if (this->empty) return std::vector<Interval>();
     const IntervalMatrix &M = this->value;
     Interval AngC = AcosOp::fwd(M(0,0));
     Interval AngS = AsinOp::fwd(M(1,0));
     /* [0,pi/2] : AngC & AngS */
     /* [0,-pi/2] : AngS & (-AngC) */
     /* [pi/2,pi] : (pi-AngS) & AngC */
     /* [-pi/2,-pi/2] : (-pi-AngS) & (-AngC) */
     /* we start from -pi and turn to pi */
     std::vector<Interval> ret;
     Interval A1 = AngS & (-AngC);
     Interval A2 = (-Interval::pi()-AngS) & (-AngC);
     if (!A2.is_empty()) { 
        if (A1.intersects(A2)) { A1 |= A2; }
        else { ret.push_back(A2); }
     }
     A2 = AngS & AngC;
     if (!A1.is_empty()) { 
        if (A1.intersects(A2)) { A2 |= A1; }
        else { ret.push_back(A1); }
     }
     A1 = (Interval::pi()-AngS) & AngC;
     if (!A2.is_empty()) { 
        if (A1.intersects(A2)) { A1 |= A2; }
        else { ret.push_back(A2); }
     }
     if (!A1.is_empty()) {
        if ((ret.size()>0) && ((Interval::two_pi() + ret[0]).intersects(A1))) {
             ret[0] = (Interval::two_pi() + ret[0]) | A1;
        } else ret.push_back(A1);
     }
     return ret;
}


std::ostream& operator<<(std::ostream& os, const SE2Base& x) {
     if (x.is_empty()) { os << "SE2:empty"; return os; }
     const IntervalMatrix &M = x.getValue();
     os << "SE2:(x" << M(0,2) << ";y" << M(1,2) << ";th";
     std::vector<Interval> LAng = x.angles();
     for (auto ang : LAng) os << ang;
     os << ")";
     return os;
}

void SE2Base::draw(Figure2DInterface &fig2D, bool drawDirection,
                const StyleProperties& sBox, const StyleProperties& sPie) {
     IntervalVector B(2);
     B[0]=this->value(0,2);
     B[1]=this->value(1,2);
     fig2D.draw_box(B,sBox);
     if (!drawDirection) return;
     Vector Cent(B.mid());
     double r = B.diam().minCoeff()*0.4;
     std::vector<Interval> LAng = this->angles();
     for (auto ang : LAng) {
          Vector P1 { r*cos(ang.lb()) , r*sin(ang.lb()) };
          Vector P2 { r*cos(ang.ub()) , r*sin(ang.ub()) };
          fig2D.draw_polyline({Cent+P1, Cent, Cent+P2},0,sPie);
          fig2D.draw_pie(B.mid(),0.4*B.diam().minCoeff(),ang,sPie);
     }
}

void SE2Base::draw(bool drawDirection, 
		const StyleProperties& sBox, const StyleProperties& sPie) {
     this->draw(*DefaultView::selected_fig(),drawDirection,sBox,sPie);
}

}
