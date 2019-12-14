#include <RcppArmadillo.h>

#include <cmath>
#include "funs.h"
#include "tree.h"
#include "info.h"
#include "rng.h"
#include <map>
#ifdef MPIBART
#include "mpi.h"
#endif

using Rcpp::Rcout;
using namespace arma;
using namespace Rcpp;

// PAVA projection for monotonicity.  Fx must be sorted in order of increasing x.
// alpha_ti gives means for each time point; need to add these back in to project to
// monotonicity correctly.

void pava(bool increasing, vec& fx){

   //======================================================
   // Setup.
   //======================================================
   int n = fx.size();
   double fxhat[n];

   // Initialize vector.
   fxhat[0] = fx[0];

   //======================================================
   // Monotone increasing projection.
   //======================================================

   if(increasing){
      for(int i = 1; i < n; i++) {
         fxhat[i] = fx[i];
         // Now sweep backwards from i and find the adjacent pool of violators,
         // i.e. points where fx[i] < fx[i+1] for an adjacent string of indices.
         // Track the running numerator and denominator of the mean fx value in the pool
         int n_sum = 1;
         double f_sum = fxhat[i];
         int j = i;

         while( (j>0) && (fxhat[j] < fxhat[j-1]) ) {
            f_sum += fxhat[j-1];
            n_sum++;
            fxhat[j-1] = f_sum/n_sum;
            j--;
         }

         // Now sweep forward and replace up through fx[i]
         // with the mean of this pool
         for(int k=j; k <= i; k++) {
            fxhat[k] = f_sum/n_sum;
         }
      }
   }

   //======================================================
   // Monotone decreasing projection.
   // Note: This is just the monotone increasing algorithm,
   // working from n to 1, instead of 1 to n.
   //======================================================

   if(!increasing){
      for(int i = (n-2); i >= 0; i--) {
         fxhat[i] = fx[i];
         // Now sweep backwards from i and find the adjacent pool of violators,
         // i.e. points where fx[i+1] < fx[i] for an adjacent string of indices.
         // Track the running numerator and denominator of the mean fx value in the pool
         int n_sum = 1;
         double f_sum = fxhat[i];
         int j = i;

        // Rcpp::Rcout << "For obs i = " << i << ":" << std::endl;
         while( (j<n) && (fxhat[j] < fxhat[j+1]) ) {
            f_sum += fxhat[j+1];
            n_sum++;
            fxhat[j+1] = f_sum/n_sum;
            //Rcpp::Rcout << j << std::endl;
            j++;
         }

         // Now sweep forward and replace up through fx[i]
         // with the mean of this pool
         for(int k=j; k>=i; k--){
            fxhat[k] = f_sum/n_sum;
         }
      }
   }

   //======================================================
   // Update fx.
   //======================================================
   for(int i=0; i<n; i++){
      fx[i] = fxhat[i];
   }
}

//--------------------------------------------------
//fit for multiple data points, not by reference, pava-projecting the entire vector before retrieving fit.
// We pass alpha_t (tlen-length vec of overall means at each time point, because projection
// needs to be done with means included, else just projecting deviations from the mean.)
void fit_pava(tree& t, xinfo& xi, dinfo& di, vec& fv, vec& alpha_t)
{
   double *xx;
   tree::tree_cp bn;
   fv.resize(di.n);

   arma::uvec id; //idx of current obs t value.

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);

      // Find index of mu (time point) corresponding to each obs.  Use this mu.
      fv[i] = bn->getm(di.t[i], di.tref);
   }
}
