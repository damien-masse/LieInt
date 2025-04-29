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
#include "contract.h"

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
   if (M(2,0).is_empty() || M(2,1).is_empty() || M(2,2).is_empty()) {
      M.set_empty(); return;
   }
   /* contract rotation */
   if (!contract_rotMatrix2(M.topLeftCorner<2,2>())) { M.set_empty(); }
}

SE2Base::LieIntervalMatrix SE2Base::representationAlgebra(const IntervalVector &V) {
   assert(V.size()==SE2Base::LieAlgebraGenerators.size());
   SE2Base::LieIntervalMatrix R = SE2Base::LieIntervalMatrix::Zero();
   for (int i=0;i<V.size();i++) {
       R += V[i]*SE2Base::LieAlgebraGenerators[i];
   }
   return R;
}

SE2Base::LieMatrix SE2Base::representationAlgebraV(const Vector &V) {
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

SE2Base SE2Base::center(bool uncertLeft) {
     if (this->empty) return SE2Base::Empty();
     IntervalMatrix Ct = this->value.mid();
     if (Ct(0,0)==0.0 && Ct(1,0)==0.0) {
          return SE2Base::Identity();
     }
     Interval a = Atan2Op::fwd(Ct(1,0),Ct(0,0));
     Ct(0,0) = Ct(1,1) = CosOp::fwd(a);
     Ct(1,0) = SinOp::fwd(a);
     Ct(0,1) = -Ct(1,0);
     /* inverse de position : position de Ct^{-1} */
     Interval p1 = -Ct(0,2)*Ct(0,0)-Ct(1,2)*Ct(1,0); 
     Interval p2 = Ct(0,2)*Ct(1,0)-Ct(1,2)*Ct(0,0);
     SE2Base::LieIntervalMatrix &M = this->value;
     if (uncertLeft) {
     /* M = M'.Ct   =>    M' = M Ct^{-1} */
       M(0,2) += M(0,0)*p1 + M(0,1)*p2;
       M(1,2) += M(1,0)*p1 + M(1,1)*p2;
    } else {
       /* M = Ct * M'  => M' = Ct^{-1} M */
       M(0,2) = M(0,2)*Ct(0,0) + M(1,2)*Ct(1,0) + p1;
       M(1,2) = M(0,2)*Ct(0,1) + M(1,2)*Ct(1,1) + p2;
    }
    M(0,0) = Ct(0,0)*M(1,1)-Ct(1,0)*M(0,1);
    M(1,0) = -Ct(1,0)*M(1,1)-Ct(0,0)*M(0,1);
    if (!contract_unitVector2(M.topLeftCorner<2,1>())) {
      M.set_empty(); return SE2Base::Empty();
    }
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

void SE2Base::draw(Figure2D &fig2D, bool drawDirection,
                const StyleProperties& sBox, const StyleProperties& sPie) const {
     IntervalVector B(2);
     B[0]=this->getX();
     B[1]=this->getY();
     fig2D.draw_box(B,sBox);
     if (!drawDirection) return;
     Vector Cent(B.mid());
     double r = B.diam().minCoeff()*0.4;
     const std::vector<Interval> LAng = this->angles();
     for (auto ang : LAng) {
          Vector P1 { r*cos(ang.lb()) , r*sin(ang.lb()) };
          Vector P2 { r*cos(ang.ub()) , r*sin(ang.ub()) };
          fig2D.draw_polyline({Cent+P1, Cent, Cent+P2},0,sPie);
          fig2D.draw_pie(Cent,r,ang,sPie);
     }
}

void SE2Base::draw(bool drawDirection, 
		const StyleProperties& sBox, const StyleProperties& sPie) const {
     this->draw(*DefaultView::selected_fig(),drawDirection,sBox,sPie);
}

inline static Eigen::Matrix<double,2,1> rot2(double theta) {
     return { cos(theta) , sin(theta) };
}

void SE2Ext::draw(Figure2D &fig2D, bool drawDirection,
          const StyleProperties& sBox, const StyleProperties& sPie) const {
     if (this->empty) return;
     auto TransL = this->left.getValue().topRows<2>().mid();
     auto &B = this->cent.getValue().topRightCorner<2,1>();
     auto &TransR = this->right.getValue().topRightCorner<2,1>();
     Interval rhoRight = SqrtOp::fwd(SqrOp::fwd(TransR[0])+SqrOp::fwd(TransR[1]));
     Interval angleRight = Interval::zero();
     if (rhoRight.ub()>0.0) angleRight = Atan2Op::fwd(TransR[1],TransR[0]);
     Interval angleAfter = Atan2Op::fwd(this->right.getValue()(1,0),
		this->right.getValue()(0,0));
     Interval angleBefore = Atan2Op::fwd(this->left.getValue()(1,0),
		this->left.getValue()(0,0));
     const std::vector<Interval> CAng = this->cent.angles();
     for (auto ang : CAng) {
          Interval sumAngle= ang+angleRight;
          IntervalVector boxUp = B+rhoRight*rot2(sumAngle.ub());
          if (cos(sumAngle.ub())>0.0) {
              Vector v1 = 
		TransL * Vector({boxUp[0].lb(), boxUp[1].ub(), 1.0});
              Vector v2 = 
		TransL * Vector({boxUp[0].ub(), boxUp[1].ub(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          } 
          else if (cos(sumAngle.ub())<0.0) {
              Vector v1 = 
		TransL * Vector({boxUp[0].lb(), boxUp[1].lb(), 1.0});
              Vector v2 = 
		TransL * Vector({boxUp[0].ub(), boxUp[1].lb(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          }
          if (sin(sumAngle.ub())>0.0) {
              Vector v1 = 
		TransL * Vector({boxUp[0].lb(), boxUp[1].lb(), 1.0});
              Vector v2 = 
		TransL * Vector({boxUp[0].lb(), boxUp[1].ub(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          } 
          else if (sin(sumAngle.ub())<0.0) {
              Vector v1 = 
		TransL * Vector({boxUp[0].ub(), boxUp[1].lb(), 1.0});
              Vector v2 = 
		TransL * Vector({boxUp[0].ub(), boxUp[1].ub(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          }
	  if (drawDirection) {
              Vector Cnt = TransL * 
			Vector({boxUp[0].mid(), boxUp[1].mid(), 1.0});
              double r = B.diam().minCoeff()*0.4;
              Interval ang2 = angleBefore + ang.ub()+angleAfter;
              Vector P1 { r*cos(ang2.lb()) , r*sin(ang2.lb()) };
              Vector P2 { r*cos(ang2.ub()) , r*sin(ang2.ub()) };
              fig2D.draw_pie(Cnt,Interval(0,r),ang2,sPie);
          }
          IntervalVector boxDown = B+rhoRight*rot2(sumAngle.lb());
          if (cos(sumAngle.lb())>0.0) {
              Vector v1 = 
		TransL * Vector({boxDown[0].lb(), boxDown[1].lb(), 1.0});
              Vector v2 = 
		TransL * Vector({boxDown[0].ub(), boxDown[1].lb(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          } 
          else if (cos(sumAngle.lb())<0.0) {
              Vector v1 = 
		TransL * Vector({boxDown[0].lb(), boxDown[1].ub(), 1.0});
              Vector v2 = 
		TransL * Vector({boxDown[0].ub(), boxDown[1].ub(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          }
          if (sin(sumAngle.lb())>0.0) {
              Vector v1 = 
		TransL * Vector({boxDown[0].ub(), boxDown[1].lb(), 1.0});
              Vector v2 = 
		TransL * Vector({boxDown[0].ub(), boxDown[1].ub(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          } 
          else if (sin(sumAngle.lb())<0.0) {
              Vector v1 = 
		TransL * Vector({boxDown[0].lb(), boxDown[1].lb(), 1.0});
              Vector v2 = 
		TransL * Vector({boxDown[0].lb(), boxDown[1].ub(), 1.0});
              fig2D.draw_polyline( { v1,v2 }, sBox);
          }
	  if (drawDirection) {
              Vector Cnt = TransL * 
			Vector({boxDown[0].mid(), boxDown[1].mid(), 1.0});
              double r = B.diam().minCoeff()*0.4;
              Interval ang2 = angleBefore + ang.lb()+angleAfter;
              Vector P1 { r*cos(ang2.lb()) , r*sin(ang2.lb()) };
              Vector P2 { r*cos(ang2.ub()) , r*sin(ang2.ub()) };
              fig2D.draw_pie(Cnt,Interval(0,r),ang2,sPie);
          }
          fig2D.draw_pie(TransL * Vector({B[0].ub(),B[1].ub(),1.0}),
			rhoRight,angleBefore+sumAngle,sBox);
          fig2D.draw_pie(TransL * Vector({B[0].lb(),B[1].lb(),1.0}),
			rhoRight,angleBefore+sumAngle,sBox);
          fig2D.draw_pie(TransL * Vector({B[0].ub(),B[1].lb(),1.0}),
			rhoRight,angleBefore+sumAngle,sBox);
          fig2D.draw_pie(TransL * Vector({B[0].lb(),B[1].ub(),1.0}),
			rhoRight,angleBefore+sumAngle,sBox);
          if (CosOp::fwd(sumAngle).contains(1.0)) 
             fig2D.draw_line( TransL* Vector({B[0].ub()+rhoRight.ub(),B[1].lb(),1.0}),
			TransL* Vector({B[0].ub()+rhoRight.ub(),B[1].ub(),1.0}) ,sBox);
          if (CosOp::fwd(sumAngle).contains(-1.0)) 
             fig2D.draw_line( TransL*Vector({B[0].lb()-rhoRight.ub(),B[1].lb(),1.0}),
			TransL*Vector({B[0].lb()-rhoRight.ub(),B[1].ub(),1.0}) ,sBox);
          if (SinOp::fwd(sumAngle).contains(1.0)) 
             fig2D.draw_line( TransL*Vector({B[0].lb(),B[1].ub()+rhoRight.ub(),1.0}),
		 	TransL*Vector({B[0].ub(),B[1].ub()+rhoRight.ub(),1.0}) ,sBox);
          if (SinOp::fwd(sumAngle).contains(-1.0)) 
             fig2D.draw_line( TransL*Vector({B[0].lb(),B[1].lb()-rhoRight.ub(),1.0}),
			 TransL*Vector({B[0].ub(),B[1].lb()-rhoRight.ub(),1.0}) ,sBox);
     }
}

void SE2Ext::draw(bool drawDirection, 
		const StyleProperties& sBox, const StyleProperties& sPie) const {
     this->draw(*DefaultView::selected_fig(),drawDirection,sBox,sPie);
}

}
