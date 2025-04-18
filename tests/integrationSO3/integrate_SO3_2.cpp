/** 
 *  \file intlie.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */


#include <codac>
#include <cmath>
#include "intLie.h"

using namespace lieInt;
using namespace codac2;

int main() { 

   ScalarVar t;

   {

   AnalyticFunction f({t}, vec(cos(t)+Interval(-0.005,0.005),ScalarExpr(0.5)+Interval(-0.005,0.005),sin(t)+Interval(-0.005,0.005)));
   TimeTube<std::vector<IntervalVector>> PT = buildPolynomialTube(f,1,Interval(0.0,10.0), 2000);

//   SO3Integrator::maxTaylor=10;
//   SO3Integrator::maxScaling=0;
   TimeTube<SO3Ext> TRes = 
	SO3Integrator::integrate_Ptube(PT, 
		SO3Ext::Identity(),true);

   for (int i=0;i<2000;i++) {
      if (i%20==19)
         std::cout << TRes.tube[i].time << " : " << 
		TRes.tube[i].value.getValueBase() << "\n";
   }
   }

   {

   AnalyticFunction f({t}, vec(cos(t),ScalarExpr(0.5),sin(t)));
   TimeTube<std::vector<IntervalVector>> PT = buildPolynomialTube(f,1,Interval(0.0,10.0), 2000);

//   SO3Integrator::maxTaylor=10;
//   SO3Integrator::maxScaling=0;
   TimeTube<SO3Ext> TRes = 
	SO3Integrator::integrate_Ptube(PT, 
                SO3Ext( {{{0.995999,1}, {-0.02,0.02}, {-0.02,0.02}},
                         {{-0.02,0.02}, {0.995999,1}, {-0.02,0.02}},
                         {{-0.02,0.02}, {-0.02,0.02}, {0.995999,1}} }),true);

   for (int i=0;i<2000;i++) {
      if (i%20==19)
         std::cout << TRes.tube[i].time << " : " << 
		TRes.tube[i].value.getValueBase() << "\n";
   }
   }


   return 0;
}
