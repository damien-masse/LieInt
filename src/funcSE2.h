/** 
 *  \file funcSE2.h functions on SE2
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once
#include <codac>
#include "lieGroup.h"
#include "funcLie.h"
#include "lieSE2.h"

using namespace codac2;


namespace lieInt {
  
  /* params for SE2 functions */
  extern ScalarVar SE2_theta;
  extern ScalarVar SE2_x;
  extern ScalarVar SE2_y;

  /* space-dependant function using parametric analytic expressions */
  class FunctionSE2Params : FunctionLie<SE2Base,SE2Ext> {
      public:
          FunctionSE2Params(const VectorExpr &expr);
          /* evaluation on one interval */
          IntervalVector eval(const SE2Base &R) const;
          IntervalVector eval(const SE2Ext &R) const;
          /* first-order bounding */
          std::vector<IntervalVector> boundFO(const SE2Ext &R) const;
      
      private:
          AnalyticFunction<VectorType> func;
  } 

}
