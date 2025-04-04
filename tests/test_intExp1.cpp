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
   IntervalMatrix M  {
      { 0.0, M_PI, {0.0,M_PI} },
      { -M_PI, 0.0, 0.0 },
      { 0.0, 0.0, 0.0 } };
   
   for (int s=0;s<=5;s++) {
     for (int n=2;n<=20;n++) {
        IntervalMatrix U = encloseIntExp(M,n,s);
        double res = std::max(U(0,2).diam(),U(1,2).diam());
        std::cout << s << ";" << n << ";" << res << "\n";
     }
   }
}
