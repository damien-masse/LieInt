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
  {
     IntervalMatrix M  {
      { {-1.0,0.0}, {0.0,0.4}, {1,4} },
      { {-0.4,-0.0}, {-1.9,0.0}, {2,3} },
      { 0.0, 0.0, 1.0 } };
   
     SE2Base E(M);

     E.draw();

     SE2Base Ct = E.center();
     SE2Base N = E*Ct;
     N.draw();
     
     SE2Ext Ext(SE2Base::Identity(),E,Ct);
     Ext.draw(true, StyleProperties(Color::red()));
  }

  {
     SE2Base C(Interval(1.2),Interval(-2.1),Interval(-1.7,-1.6));

     SE2Base M(Interval(-0.3,0.3), Interval(-0.4,0.2), Interval(-0.9,1.5));
 
     SE2Ext Ext(C,M,SE2Base::Identity());

     Ext.draw(true,StyleProperties(Color::blue()));
  }

  {
     SE2Base C = SE2Base::Identity();

     SE2Base M(Interval(-1.0,1.0), Interval(-1.0,1.0), Interval(0.0,0.0));
 
     SE2Ext Ext(SE2Base::Identity(),M,C);


//     (M*C).draw();
     Ext.draw(true,StyleProperties(Color::red()));
  }

}
