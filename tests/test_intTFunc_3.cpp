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
#include "../intLie.h"

using namespace lieInt;
using namespace codac2;

int main() { 

   Figure3D f3D("test_intTFunc_3");

   f3D.draw_axes(2.0);
   ScalarVar t;

   AnalyticFunction f({t}, vec(ScalarExpr(0.0),ScalarExpr(Interval(0.9,1.0)),ScalarExpr(0.0)));
   TimeTube<IntervalVector> T = buildRegularTube(f,Interval(0.0,10.0),3000);

   TimeTube<SO3Ext> TRes = 
	SO3Integrator::integrate_tube(T, 
		SO3Ext::Identity());

   for (int i=0;i<3000;i++) {
//     std::cout << TRes.tube[i].time << " : " << TRes.tube[i].value << "\n";
        
     if (i%30==0) {
        Vector place { std::cos(i/100.0)*10.0, std::sin(i/100.0)*10.0 ,
				i/100.0 };
        TRes.tube[i].value.draw3D(f3D,place, 3.0, { Color::red(0.9) });
     }
   }


   return 0;
}
