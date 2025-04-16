/** 
 *  \file intlie.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */


#include <codac>
#include "intLie.h"

using namespace lieInt;
using namespace codac2;

int main() { 

   DefaultView::set_window_properties({600,600},{300,300});

   ScalarVar t;

   AnalyticFunction f({t}, vec(M_PI*pow(sin(t),2.0),ScalarExpr(1.0),ScalarExpr(0.0)));
   TimeTube<IntervalVector> T = buildRegularTube(f,Interval(0.0,10.0),3000);

   TimeTube<SE2Ext> TRes = 
	SE2Integrator::integrate_tube(T, 
		SE2Ext::Identity());

   for (int i=0;i<3000;i++) {
     std::cout << TRes.tube[i].time << " : " << TRes.tube[i].value << "\n";
     if (i%10==0) TRes.tube[i].value.draw(false);
   }


   return 0;
}
