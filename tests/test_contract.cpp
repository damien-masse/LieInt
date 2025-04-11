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
     Eigen::Matrix<Interval,2,2> M  {
     { {0.3,0.6}, {-0.9,-0.75} },
      { {0.75,0.9}, {0.2,0.7} } };

     std::cout << contract_rotMatrix2(M) << "\n";

     std::cout << M << "\n";


     contract_unitVector2(M.col(1));
 
     std::cout << M << "\n";
}
