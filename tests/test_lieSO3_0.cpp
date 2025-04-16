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
  {
     Figure3D f3D("test_lieSO3_0");
     f3D.draw_axes(2.0);
     for (double i=-1.0;i<=1.0;i+=0.25) 
     for (double j=0.0;j<=1.0;j+=0.125) 
     for (double k=-1.0;k<=1.0;k+=0.25) {
        double thetac = -M_PI/2.0+j*M_PI;
        double phic = i*M_PI;
        double psic = k*M_PI;
        SO3Base::CardanAngles ang {
         psic+Interval(-0.1,0.1),
         thetac+Interval(-0.05,0.05),
         phic+Interval(-0.1,0.1) };
         SO3Base E(ang);
         if (j!=0.5) 
		E.draw3D(f3D,{ i*10.0, k*10.0, j*20.0 }, 1.5);
         else { 
                std::cout << E << "\n";
		E.draw3D(f3D,{ i*10.0, k*10.0, j*20.0 }, 1.5,
			{ Color::red(0.9) });
         }
     } 
   

  }
}
