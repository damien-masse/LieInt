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
  {
     IntervalMatrix M  {
      { {-1.0,0.0}, {0.0,0.4}, {0.0,M_PI} },
      { {-0.4,-0.0}, {-1.9,0.0}, {0.0,1.0} },
      { 0.0, 0.0, 1.0 } };
   
     SE2Base E(M);
     std::cout << E << "\n";

     E.draw();
  }

  {
     IntervalMatrix M  {
      { {-1.0,0.9}, {-0.9,0.9}, {0.0,M_PI} },
      { {-0.8,0.9}, {-1.0,0.9}, {-2.0,-1.0} },
      { 0.0, 0.0, 1.0 } };
   
     SE2Base E(M);
     std::cout << E << "\n";

     E.draw(true,StyleProperties::inside(),StyleProperties(Color::red()));
  }

}
