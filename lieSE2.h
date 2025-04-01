/** 
 *  \file lieSE2.h  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <codac>

using namespace codac2;

namespace lieInt
{
    class SE2Base : public LieBaseMatrix<3,3> {

    };

    class SE2Ext : public LieExtMatrix<3,3,SE2Base> {
    
    };

}
