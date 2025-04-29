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

   ScalarVar &t = lie_time;

   SE2Ext currentpos = SE2Ext::Identity();

   for (int i=0;i<3;i++) {

   {
     FuncTAnalytic f(vec(ScalarExpr(Interval(0.95,1.05)),ScalarExpr(Interval(0.95,1.05)),ScalarExpr(Interval(-0.05,0.05))));
     TimeTube<IntervalVector> T = buildRegularTube(f,Interval(0.0,1.0),201);

     TimeTube<SE2Ext> TRes = 
  	SE2Integrator::integrate_tube(T, 
		currentpos, true);

     for (int i=0;i<201;i++) {
//     if (i%10==0) std::cout << TRes.tube[i].time << " : " << TRes.tube[i].value << "\n";
       if (i%10==0) TRes.tube[i].value.draw(false);
     }
     std::cout << TRes.tube[200].value.getCenter().getValue() << "\n";
     currentpos=TRes.tube[200].value.getLeft()*TRes.tube[200].value.getRight();
     std::cout << currentpos << "\n";
   }

   {
     FuncTAnalytic f(vec(ScalarExpr(Interval(-1.05,-0.95)),ScalarExpr(Interval(0.95,1.05)),ScalarExpr(Interval(-0.05,0.05))));
     TimeTube<IntervalVector> T = buildRegularTube(f,Interval(0.0,1.0),201);

     TimeTube<SE2Ext> TRes = 
  	SE2Integrator::integrate_tube(T, 
		currentpos, true);

     for (int i=0;i<201;i++) {
//     if (i%10==0) std::cout << TRes.tube[i].time << " : " << TRes.tube[i].value << "\n";
       if (i%10==0) TRes.tube[i].value.draw(false);
     }
     std::cout << TRes.tube[200].value.getCenter().getValue() << "\n";
     currentpos=TRes.tube[200].value.getLeft()*TRes.tube[200].value.getRight();
     std::cout << currentpos << "\n";
   }

   }

   return 0;
}
