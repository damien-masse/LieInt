/** 
 *  \file funcLie.h functions from Lie group to Lie algebra
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


namespace lieInt {

  /* variable time for functions */
  extern ScalarVar lie_time;

  /* time-dependant only function : R -> LieAlg */
  class FunctionTAlg {
    public:
      virtual ~FunctionTAlg() = default;
      /* evaluation on one interval */
      virtual IntervalVector eval(const Interval &T) const = 0;
      /* polynomial taylor approximation, degree d */
      virtual std::vector<IntervalVector> eval(const Interval &T, int d) const = 0;
  };

  /* space-dependant only function : SE2 -> LieAlg(T) */
  template <class TBase, class TExt>
  class FunctionLie {
    public:
      virtual ~FunctionLie() = default;
      /* evaluation on one space interval */
      virtual IntervalVector eval(const TBase &R) const = 0;
      virtual IntervalVector eval(const TExt &R) const = 0;
      /* first-order bounding : [R0] and [R1] such that
         forall r \in R.Left R.Cent, f(LeftCent) \in [R0] + [R1](Cent-Id) */
      virtual std::vector<IntervalVector> boundFO(const TExt &R) const = 0;
  };

  /* full-dependant only function : SE2 -> LieAlg(T) */
  template <class TBase, class TExt>
  class FunctionTLie {
    public:
      virtual ~FunctionTLie() = default;
      /* evaluation on one space interval */
      virtual IntervalVector eval(const Interval &T, const TBase &R) const = 0;
      virtual IntervalVector eval(const Interval &T, const TExt &R) const = 0;
      /* first-order bounding : [R0] and [R1] such that
         forall r \in R.Left R.Cent, f(LeftCent) \in [R0] + [R1](Left-Id) */
      virtual std::vector<IntervalVector> boundFO(const Interval &T,
						const TExt &R) const = 0;
  };

  /* time-dependant function from analytic function expression */
  class FuncTAnalytic : public FunctionTAlg {
      public:
         FuncTAnalytic(const VectorExpr &expr);
         /* evaluation on one interval */
         IntervalVector eval(const Interval &T) const;
         /* polynomial taylor approximation, degree d */
         std::vector<IntervalVector> eval(const Interval &T, int d) const;
      
      private:
         AnalyticFunction<VectorType> func;
  };

  inline FuncTAnalytic::FuncTAnalytic(const VectorExpr &expr) :
	func({lie_time},expr) 
  { }

  inline IntervalVector FuncTAnalytic::eval(const Interval &T) const
  { 
     return func.eval(T);
  }

  inline std::vector<IntervalVector> FuncTAnalytic::eval(const Interval &T, int d) const
  { 
     if (d==0) { std::vector<IntervalVector> ret;
		 ret.push_back(func.eval(T)); return ret; }
     assert(d==1);
     std::vector<IntervalVector> ret;
     ret.push_back (func.eval(EvalMode::NATURAL,T.lb()));
     ret.push_back (func.diff(T));
     return ret;
  }

  
}
