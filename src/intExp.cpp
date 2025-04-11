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


/* Approche de la borne :
   pour T_n = sum(i=1^d) Bi T{n-i} / n
   on fixe Q >= 2 (Bi/n)^{1/i}
   Si on cherche Q < \tau, on obtient n >= Bi/((\tau/2)^i)
   donc : 
      => calcul de n en fonction de tau,
         détermination et calcul de Q.
         retour ...
*/
static int calculNQ(const std::vector<Interval>& norms, double tau, Interval &Q) {
   double tau2=tau/2.0;
   double ptau2=tau2;
   int d = norms.size();
   int n=d;
   for (int i=0;i<d;i++) {
       int n2 = (norms[i].ub())/ptau2;
       if (n2>=n) n=n2+1;
       ptau2*=tau2;
   }
   Q=0.0;
   
   for (int i=0;i<d;i++) {
      Q = max(Q,2.0*root(norms[i]/n,i+1));
   }  
   return n;
}

IntervalMatrix encloseIntExp(const std::vector<IntervalMatrix>& P,
                                         int n, double tau=0.0) {
   const int sizeP = P.size();
   assert(sizeP>=0);
   const int dim = P[0].rows();
   /* determination de Q */
   std::vector<Interval> norms(sizeP);
   for (int i=0;i<sizeP;i++) norms[i]=infinite_norm(P[i]);
   Interval Q;
   int nbstepsmin = calculNQ(norms,0.9,Q);
   Q = Q.ub();

   std::vector<IntervalMatrix> terms = std::vector<IntervalMatrix>(sizeP,
		IntervalMatrix::Zero(dim,dim));
   terms[0]=IntervalMatrix::Identity(dim,dim);
   IntervalMatrix Somme = IntervalMatrix::Zero(dim,dim);
   int i=1;
   Interval maxnorm(1.0);
   std::cout << "nbstepsmin : " << nbstepsmin << "\n";
   while ((i<=nbstepsmin) || (maxnorm.ub()>=tau && i<=n)) {
        maxnorm = Interval(0.0);
        Interval K(1.0);
        K /= (double)i;
        int place = i%sizeP;
        if (i>=sizeP) {
           terms[place] = terms[place]*P[sizeP-1];
        }
        int start=std::min(sizeP-2,i-1);
        Interval PQ=Q;
        for (int j=0;j<=start;j++) {
            int p2 = (i-j-1)%sizeP;
            p2 = p2%sizeP;
            terms[place] += terms[p2]*P[j];
            maxnorm=max(maxnorm,infinite_norm(terms[p2])/PQ);
            PQ=Q*PQ;
        }
        terms[place] *= K;
        maxnorm=max(maxnorm,infinite_norm(terms[place]));
        Somme = Somme+terms[place];
	i=i+1;
   }
   std::cout << "i final :" << i << "\n";
   std::cout << "Q :" << Q << "\n";
   std::cout << "maxnorm :" << maxnorm << "\n";
   double eps = ((maxnorm*Q)/(1-Q)).ub();
   Somme += IntervalMatrix::Constant(dim,dim,Interval(-eps,eps));
//   std::cout << "Somme : " << Somme << "\n";
   return Somme+IntervalMatrix::Identity(dim,dim);
}

}


