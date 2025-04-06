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

const std::vector<SE2Base::LieMatrix> SE2Base::LieAlgebraGenerators
     { { {0,-1,0}, {1,0,0}, {0,0,0} },
       { {0,0,1},  {0,0,0}, {0,0,0} },
       { {0,0,0},  {0,0,1}, {0,0,0} } };

void SE2Base::contractor() {
   if (this->empty) return;
   SE2Base::LieIntervalMatrix &M = this->value;
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

SE2Base::LieIntervalMatrix SE2Base::representationAlgebra(const IntervalVector &V) {
   assert(V.size()==SE2Base::LieAlgebraGenerators.size());
   SE2Base::LieIntervalMatrix R = SE2Base::LieIntervalMatrix::Zero();
   for (int i=0;i<V.size();i++) {
       R += V[i]*SE2Base::LieAlgebraGenerators[i];
   }
   return R;
}

SE2Base::LieMatrix SE2Base::representationAlgebra(const Vector &V) {
   assert(V.size()==SE2Base::LieAlgebraGenerators.size());
   SE2Base::LieMatrix R = SE2Base::LieMatrix::Zero();
   for (int i=0;i<V.size();i++) {
       R += V[i]*SE2Base::LieAlgebraGenerators[i];
   }
   return R;
}

SE2Base SE2Base::inverse() const {
   if (this->empty) return SE2Base::Empty();
   SE2Base::LieIntervalMatrix M = this->value;
   M(0,1) = -M(0,1);
   M(1,0) = -M(1,0);
   M.block<2,1>(0,2) = -M.topLeftCorner<2,2>()*M.block<2,1>(0,2);
   return SE2Base(M);
}

void SE2Base::inverse_inplace() {
   if (this->empty) return;
   SE2Base::LieIntervalMatrix &M = this->value;
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

SE2Base SE2Base::center() {
     if (this->empty) return SE2Base::Empty();
     IntervalMatrix Ct = this->value.mid();
     if (Ct(0,0)==0.0 && Ct(1,0)==0.0) {
          return SE2Base::Identity();
     }
     Interval a = Atan2Op::fwd(Ct(1,0),Ct(0,0));
     Ct(0,0) = Ct(1,1) = CosOp::fwd(a);
     Ct(1,0) = SinOp::fwd(a);
     Ct(0,1) = -Ct(1,0);
     /* M = M'.Ct   =>    M' = M Ct^{-1} */
     SE2Base::LieIntervalMatrix &M = this->value;
     Interval p1 = -Ct(0,2)*Ct(0,0)-Ct(1,2)*Ct(1,0);
     Interval p2 = Ct(0,2)*Ct(1,0)-Ct(1,2)*Ct(0,0);
     M(0,2) += M(0,0)*p1 + M(0,1)*p2;
     M(1,2) += M(1,0)*p1 + M(1,1)*p2;
     M(0,0) = Ct(0,0)*M(1,1)-Ct(1,0)*M(0,1);
     M(1,0) = -Ct(1,0)*M(1,1)-Ct(0,0)*M(0,1);
     SqrOp::bwd(1.0-SqrOp::fwd(M(1,0)),M(0,0));
     SqrOp::bwd(1.0-SqrOp::fwd(M(0,0)),M(1,0));
     M(0,1) = -M(1,0);
     M(1,1) = M(0,0);
     return SE2Base(Ct);
}

const std::vector<Interval> SE2Base::angles() const {
     if (this->empty) return std::vector<Interval>();
     const SE2Base::LieIntervalMatrix &M = this->value;
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
     const SE2Base::LieIntervalMatrix &M = x.getValue();
     os << "SE2:(x" << M(0,2) << ";y" << M(1,2) << ";th";
     const std::vector<Interval> LAng = x.angles();
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
     const std::vector<Interval> LAng = this->angles();
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

inline static Vector rot2(double theta) {
     return { cos(theta) , sin(theta) };
}

void SE2Ext::draw(Figure2DInterface &fig2D, bool drawDirection,
          const StyleProperties& sBox, const StyleProperties& SPie) {
     if (this->empty) return;
     IntervalVector B(2);
     B[0]=this->left.getValue()(0,2);
     B[1]=this->left.getValue()(1,2);
     IntervalVector Trans(2);
     Trans[0]=this->cent.getValue()(0,2);
     Trans[1]=this->cent.getValue()(1,2);
     Interval angleCent = Atan2Op::fwd(Trans[1],Trans[0]);
     Interval rhoCent = SqrtOp::fwd(SqrOp::fwd(Trans[0])+SqrOp::fwd(Trans[1]));
     const std::vector<Interval> LAng = this->left.angles();
     for (auto ang : LAng) {
          Interval sumAngle= ang+angleCent;
          IntervalVector boxUp = B+rhoCent*rot2(sumAngle.ub());
          if (cos(sumAngle.ub())>0.0) {
              fig2D.draw_polyline({
				{boxUp[0].lb(), boxUp[1].ub()},
				{boxUp[0].ub(), boxUp[1].ub()} },0,sBox);
          } 
          else if (cos(sumAngle.ub())<0.0) {
              fig2D.draw_polyline({
				{boxUp[0].lb(), boxUp[1].lb()},
				{boxUp[0].ub(), boxUp[1].lb()} },0,sBox);
          }
          if (sin(sumAngle.ub())>0.0) {
              fig2D.draw_polyline({
				{boxUp[0].lb(), boxUp[1].lb()},
				{boxUp[0].lb(), boxUp[1].ub()} },0,sBox);
          } 
          else if (sin(sumAngle.ub())<0.0) {
              fig2D.draw_polyline({
				{boxUp[0].ub(), boxUp[1].lb()},
				{boxUp[0].ub(), boxUp[1].ub()} },0,sBox);
          }
          IntervalVector boxDown = B+rhoCent*rot2(sumAngle.lb());
          if (cos(sumAngle.lb())>0.0) {
              fig2D.draw_polyline({
				{boxDown[0].lb(), boxDown[1].lb()},
				{boxDown[0].ub(), boxDown[1].lb()} },0,sBox);
          } 
          else if (cos(sumAngle.lb())<0.0) {
              fig2D.draw_polyline({
				{boxDown[0].lb(), boxDown[1].ub()},
				{boxDown[0].ub(), boxDown[1].ub()} },0,sBox);
          }
          if (sin(sumAngle.lb())>0.0) {
              fig2D.draw_polyline({
				{boxDown[0].ub(), boxDown[1].lb()},
				{boxDown[0].ub(), boxDown[1].ub()} },0,sBox);
          } 
          else if (sin(sumAngle.lb())<0.0) {
              fig2D.draw_polyline({
				{boxDown[0].lb(), boxDown[1].lb()},
				{boxDown[0].lb(), boxDown[1].ub()} },0,sBox);
          }
          fig2D.draw_pie(B.ub(),rhoCent,sumAngle,sBox);
          fig2D.draw_pie(B.lb(),rhoCent,sumAngle,sBox);
          fig2D.draw_pie({B[0].ub(),B[1].lb()},rhoCent,sumAngle,sBox);
          fig2D.draw_pie({B[0].lb(),B[1].ub()},rhoCent,sumAngle,sBox);
          if (CosOp::fwd(sumAngle).contains(1.0)) 
             fig2D.draw_polyline({
				{B[0].ub()+rhoCent.ub(),B[1].lb()},
				{B[0].ub()+rhoCent.ub(),B[1].ub()}
				 },0,sBox);
          if (CosOp::fwd(sumAngle).contains(-1.0)) 
             fig2D.draw_polyline({
				{B[0].lb()-rhoCent.ub(),B[1].lb()},
				{B[0].lb()-rhoCent.ub(),B[1].ub()}
				 },0,sBox);
          if (SinOp::fwd(sumAngle).contains(1.0)) 
             fig2D.draw_polyline({
				{B[0].lb(),B[1].ub()+rhoCent.ub()},
				{B[0].ub(),B[1].ub()+rhoCent.ub()}
				 },0,sBox);
          if (SinOp::fwd(sumAngle).contains(-1.0)) 
             fig2D.draw_polyline({
				{B[0].lb(),B[1].lb()-rhoCent.ub()},
				{B[0].ub(),B[1].lb()-rhoCent.ub()}
				 },0,sBox);
     }
}

void SE2Ext::draw(bool drawDirection, 
		const StyleProperties& sBox, const StyleProperties& sPie) {
     this->draw(*DefaultView::selected_fig(),drawDirection,sBox,sPie);
}

}
