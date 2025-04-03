/** 
 *  \file intExp.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Damien Massé
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <codac>
#include <cmath>

using namespace codac2;

namespace lieInt
{

/* norm computation */
inline static Interval infinite_norm(const IntervalMatrix &A) {
   return A.cwiseAbs().rowwise().sum().maxCoeff();
}

IntervalMatrix encloseIntExp(const IntervalMatrix& A, int n, double tau) {
   int dim=A.rows();
   IntervalMatrix expA = IntervalMatrix::Zero(dim,dim);
   double nrmA = infinite_norm(A).ub();
   assert(nrmA<n);
   IntervalMatrix Tn=A; /* successive terms */
   int nbsteps=1;
   double computedtau = tau+1;
   while (nbsteps<=n && (nrmA<(nbsteps+1) || computedtau>tau)) {
     expA += Tn;
     nbsteps++;
     Tn *= (Interval(1.0)/nbsteps)*A;
     if (nrmA<(nbsteps+1))
        computedtau = infinite_norm(Tn).ub()/(1.0-nrmA/(nbsteps+1));
   }
   Interval eps=infinite_norm(Tn)/(1.0-infinite_norm(A)/(nbsteps+1));
   expA.inflate(eps.ub());
   return expA+IntervalMatrix::Identity(dim,dim);
}

IntervalMatrix encloseIntExp(const IntervalMatrix& A, int n, int S,
                               double tauSS=0.0, double tauExp=0.0) {
   int nbSS=0;
   double nrmA = infinite_norm(A).ub();
   IntervalMatrix dA=A;
   while (nbSS<S && nrmA>tauSS) {
      dA *= 0.5;
      nbSS++;
      nrmA *= 0.5;
   }
   IntervalMatrix expA = encloseIntExp(dA,n,tauExp);
   for (int i=0;i<nbSS;i++) {
      expA = expA*expA;
   }
   return expA;
}

}
