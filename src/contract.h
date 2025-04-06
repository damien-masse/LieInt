/** 
 *  \file contract.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <codac>

using namespace codac2;

namespace lieInt
{

  /** \brief constraint norm_2(A) = 1.0, specialised version
   * 
   *  \param A a vector of size 2
   *
   *  \return true if non empty, false if empty
   */
  bool contract_unitVector2(Eigen::Vector<Interval,2>& A);

  /** \brief constraint norm_2(A) = 1.0, generic version
   * 
   *  \param A a vector of size n
   *
   *  \return true if non empty, false if empty
   */
  bool contract_unitVector(IntervalVector& A);

  /** \brief constraint MMt = Id, specialised version
   * 
   *  \param A a matrix of size 2*2
   *
   *  \return true if non empty, false if empty
   */
  bool contract_rotMatrix2(Eigen::Matrix<Interval,2,2>& A);

  /** \brief constraint MMt = Id
   * 
   *  \param A a matrix of size n*n
   *
   *  \return true if non empty, false if empty
   */
  bool contract_rotMatrix(IntervalMatrix& A);

}
