/** 
 *  \file funcSE2.h functions on SE2
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <codac>
#include "lieGroup.h"
#include "funcLie.h"
#include "lieSE2.h"
#include "funcSE2.h"

using namespace codac2;


namespace lieInt {
  
/* params for SE2 functions */
ScalarVar SE2_theta;
ScalarVar SE2_x;
ScalarVar SE2_y;

FunctionSE2Params::FunctionSE2Params(const VectorExpr &expr) :
   func({SE2_theta,SE2_x,SE2_y},expr)
{
}

IntervalVector FunctionSE2Params::eval(const SE2Base &R) const {
   const std::vector<Interval> angles = R.angles();
   const Interval& X = R.getX();
   const Interval& Y = R.getY();
   IntervalVector ret = IntervalVector::Constant(3,Interval::empty());
   for (auto ang : angles) {
       ret |= this->func(ang,X,Y);
   }
   return ret;
}

IntervalVector FunctionSE2Params::eval(const SE2Ext &R) const {
   const std::vector<Interval> angles = R.getLeft().angles();
   const Interval& X = R.getLeft().getX();
   const Interval& Y = R.getLeft().getY();
   Interval angAfter = Atan2Op::fwd(R.getCenter().getValue()(1,0),
                R.getCenter().getValue()(0,0));
   IntervalVector ret = IntervalVector::Constant(3,Interval::empty());
   const Interval &X2 = R.getCenter().getX();
   const Interval &Y2 = R.getCenter().getY();
   for (auto ang : angles) {
       ret |= this->func(ang+angAfter,X+X2*CosOp::fwd(ang)-Y2*SinOp::fwd(ang),
			Y+X2*SinOp::fwd(ang)+Y2*CosOp::fwd(ang));
   }
   return ret;
}

std::vector<IntervalVector> FunctionSE2Params::boundFO(const SE2Ext &R) const {
   IntervalVector enclosure=this->eval(R);
   

}

}
