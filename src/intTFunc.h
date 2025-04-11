/** 
 *  \file intTFunc.h integration of time-dependant functions
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once
#include <codac>
#include "lieGroup.h"
#include "intExp.h"
#include "lieSE2.h"

using namespace codac2;


namespace lieInt {

  /* tube simulation */
  template <class V>
  struct TimeSlice {
     Interval time;
     V value;
  }; 

  template <class V>
  struct TimeTube {
     Interval timeRange;
     std::vector<struct TimeSlice<V>> tube;
  };

  template<class T> 
  inline TimeTube<IntervalVector>
		 buildRegularTube(const AnalyticFunction<T> &f, 	
			const Interval &range, unsigned int nbslices) {
     TimeTube<IntervalVector> ret = { range , 
		std::vector<TimeSlice<IntervalVector>>(nbslices) };
     double time = range.lb();
     for (unsigned int i=0;i<nbslices;i++) {
        double nexttime = (i<nbslices-1 ?
		range.lb()+(range.diam()/nbslices)*(i+1) : range.ub());
        ret.tube[i].time=Interval(time, nexttime);
        ret.tube[i].value=f.eval(EvalMode::NATURAL,ret.tube[i].time);
        time=nexttime;
     }
     return ret;
  }

  template<class T> 
  inline TimeTube<std::vector<IntervalVector>>
		 buildPolynomialTube(const AnalyticFunction<T> &f, int degree,
			const Interval &range, unsigned int nbslices) {
     assert(degree==1);
     TimeTube<std::vector<IntervalVector>> ret = { range , 
		std::vector<TimeSlice<std::vector<IntervalVector>>>(nbslices) };
     double time = range.lb();
     for (unsigned int i=0;i<nbslices;i++) {
        double nexttime = (i<nbslices-1 ?
		range.lb()+(range.diam()/nbslices)*(i+1) : range.ub());
        ret.tube[i].time=Interval(time, nexttime);
        ret.tube[i].value.push_back
			(f.eval(EvalMode::NATURAL,ret.tube[i].time.lb()));
        ret.tube[i].value.push_back
			(f.diff(ret.tube[i].time));
        time=nexttime;
     }
     return ret;
  }

  template <class LEG>
//  requires std::is_base_of_v<LieExtMatrix,LEG>
  class odeIntegrator {
     public:
       static TimeTube<LEG> integrate_tube(const TimeTube<IntervalVector> &TIn,
                const LEG &startvalue, bool gates=false);
       static TimeTube<LEG> integrate_Ptube(const TimeTube<std::vector<IntervalVector>> &TIn, const LEG &startvalue, bool gates=false);

       static int maxScaling;
       static int maxTaylor;
       static double tauScaling;
       static double tauTaylor;
  };



  template <class LEG>
  inline TimeTube<LEG> odeIntegrator<LEG>::integrate_tube(const TimeTube<IntervalVector> &TIn, const LEG &startvalue, bool gates) {
      unsigned int nbslices = TIn.tube.size();
      TimeTube<LEG> ret = { TIn.timeRange , std::vector<TimeSlice<LEG>>(nbslices) };
      auto currentVal = LEG::Identity();
      for (unsigned int i=0;i<nbslices;i++) {
         Interval timeval=TIn.tube[i].time - TIn.tube[i].time.lb();
         Interval timediam=TIn.tube[i].time.diam();
         {
           IntervalMatrix repAlgIn = 
		LEG::Base::representationAlgebra(timeval*TIn.tube[i].value);
           IntervalMatrix valIn =
		encloseIntExp(repAlgIn, maxTaylor, maxScaling,
			tauScaling, tauTaylor);
           LEG bIn(valIn);
           ret.tube[i].time=TIn.tube[i].time;
           if (!gates) ret.tube[i].value=(startvalue*currentVal)*bIn;
         }
         {
           IntervalMatrix repAlgOut = 
		LEG::Base::representationAlgebra(TIn.tube[i].value*timediam);
           IntervalMatrix valOut =
		encloseIntExp(repAlgOut, maxTaylor, maxScaling,
			tauScaling, tauTaylor);
           LEG bIn(valOut);
           currentVal.rightProd(bIn);
           if (gates) ret.tube[i].value=startvalue*currentVal;
         }
      }
      return ret;
   } 

  template <class LEG>
  inline TimeTube<LEG> odeIntegrator<LEG>::integrate_Ptube(const TimeTube<std::vector<IntervalVector>> &TIn, const LEG &startvalue, bool gates) {
      unsigned int nbslices = TIn.tube.size();
      TimeTube<LEG> ret = { TIn.timeRange , std::vector<TimeSlice<LEG>>(nbslices) };
      auto currentVal = LEG::Identity();
      int degPol=TIn.tube[0].value.size();
      for (unsigned int i=0;i<nbslices;i++) {
         Interval timeval=TIn.tube[i].time - TIn.tube[i].time.lb();
         Interval timediam=TIn.tube[i].time.diam();
         {
	   std::vector<IntervalMatrix> repPoly(degPol);
           Interval tval=timeval;
	   for (int z=0;z<degPol;z++) {
 	      repPoly[z] = 
		LEG::Base::representationAlgebra(tval*TIn.tube[i].value[z]);
              tval = tval*timeval;
           }
           IntervalMatrix valIn =
		encloseIntExp(repPoly, maxTaylor, tauTaylor);
           LEG bIn(valIn);
           ret.tube[i].time=TIn.tube[i].time;
           if (!gates) ret.tube[i].value=(startvalue*currentVal)*bIn;
         }
         {
	   std::vector<IntervalMatrix> repPoly(degPol);
           Interval tval=timediam;
	   for (int z=0;z<degPol;z++) {
 	      repPoly[z] = 
		LEG::Base::representationAlgebra(tval*TIn.tube[i].value[z]);
              tval = tval*timediam;
           }
           IntervalMatrix valOut =
		encloseIntExp(repPoly, maxTaylor, tauTaylor);
           LEG bIn(valOut);
           currentVal.rightProd(bIn);
           if (gates) ret.tube[i].value=startvalue*currentVal;
         }
      }
      return ret;
   } 


   template <class T>
   int odeIntegrator<T>::maxTaylor=20;
   template <class T>
   int odeIntegrator<T>::maxScaling=4;
   template <class T>
   double odeIntegrator<T>::tauScaling=1e-6;
   template <class T>
   double odeIntegrator<T>::tauTaylor=1e-9;


   using SE2Integrator = odeIntegrator<SE2Ext>;

}
