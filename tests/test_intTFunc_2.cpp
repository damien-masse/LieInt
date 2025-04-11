/** 
 *  \file intlie.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */


#include <codac>
#include "../intLie.h"

using namespace lieInt;
using namespace codac2;

int main() { 

   DefaultView::set_window_properties({600,600},{300,300});

   ScalarVar t;

   AnalyticFunction f({t}, vec(M_PI*pow(sin(t),2.0),ScalarExpr(1.0),ScalarExpr(0.0)));
   TimeTube<IntervalVector> T = buildRegularTube(f,Interval(0.0,10.0),3001);

   TimeTube<SE2Ext> TRes = 
	SE2Integrator::integrate_tube(T, 
		SE2Ext::Identity(), true);

   for (int i=0;i<3001;i++) {
     if (i%100==0) std::cout << TRes.tube[i].time << " : " << TRes.tube[i].value << "\n";
     if (i%100==0) TRes.tube[i].value.draw(false);
   }

   TimeTube<std::vector<IntervalVector>> PT = buildPolynomialTube(f,1,Interval(0.0,10.0), 301);


   TimeTube<SE2Ext> TRes2 = 
	SE2Integrator::integrate_Ptube(PT, 
		SE2Ext::Identity(), true);

   for (int i=0;i<301;i++) {
     if (i%10==0) {
         std::cout << TRes2.tube[i].time << " : " << TRes2.tube[i].value << "\n";
     }
      if (i%10==0) TRes2.tube[i].value.draw(false,StyleProperties(Color::red()));
   }


   return 0;
}
