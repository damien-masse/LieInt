/** 
 *  \file intExp.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <codac>

namespace lieInt
{

  /** \brief computation of IntExp(A) with a truncated Taylor expansion
   * 
   *  \param A a square matrix
   *  \param n (max) order of the truncated expansion
   *  \param tau precision of the expansion (stops when ||A||^k/(k!) < tau).
   *         (default 0.0)
   *
   *  \return the enclosure of IntExp(A)
   */
  IntervalMatrix encloseIntExp(const IntervalMatrix& A, int n, double tau=0.0);

  /** \brief computation of IntExp(A) with a truncated Taylor expansion and
   *  scaling and squaring
   * 
   *  \param A a square matrix
   *  \param n (max) order of the truncated expansion
   *  \param S (max) number of scaling and squaring
   *  \param tauSS precision of scaling and squaring 
   *     (no step when ||A|| < tauSS) (default 0.0)
   *  \param tauExp precision of the expansion (stops when ||A||^k/(k!) < tau).
   *         (default 0.0)
   *
   *  \return the enclosure of IntExp(A)
   */
  IntervalMatrix encloseIntExp(const IntervalMatrix& A, int n, int S,
		   	       double tauSS=0.0, double tauExp=0.0);

}
