#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>
#include "tree.h"
#include "info.h"
#include "rng.h"
#include "GIGrvg.h"

using namespace arma;

//--------------------------------------------------
// PAVA projection for monotonicity.
void pava(bool increasing, vec& fx);

//--------------------------------------------------
// New fit_i_pava function.  This is all we need for pava updates, since the projection
// is done post-hoc.
// Notes:
// 1. We input the alpha_t vector since projection needs to be done on the
// re-centered vector, else just projecting differences from overall mean.
// 2. We need the eta scaling factor, since need to project after multiplying.
//--------------------------------------------------
template<class T>
double fit_i_pava(T i, std::vector<tree>& t, xinfo& xi, dinfo& di, vec& alpha_t, double& muscale, bool& incr)
{
   double *xx;
   double fv = 0.0;
   tree::tree_cp bn;
   xx = di.x + i*di.p;

   // Define a new mu vector to track the entire t trajectory for the corresponding xi.
   // We will track the trajectory, summing across trees in tree vector t,
   // then we project.  Begin with alpha_t, and add in tree fits.
   vec mu = alpha_t;

   for (size_t j=0; j<t.size(); ++j) {
      bn = t[j].bn(xx,xi);  // Get bottom node.
      mu += bn -> getm();   // Get mu vec for bottom node, and add to running sum.
   }

   // Multiply mu by scaling factor, and perform pava projection.
   mu = mu * muscale;
   pava(incr, mu);

   // Find vec index and return fitted value (fv).
   uvec idx = find(di.t[i]==di.tref);
   fv = as_scalar(mu(idx));

   return fv;
}

