/* This is called by the wrapper function tsbart().  Documentation retained at bottom of function.*/
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "rtn.hpp"
#include "pava.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List tsbartFit(arma::vec y,
               arma::vec tgt,
               arma::vec tpred,
               arma::vec x,
               arma::vec xpred,
               List xinfo_list,
               int nburn, int nsim, int ntree=200,
               double lambda=-999, double sigq=.9, double sighat=-999,
               double nu=3, double ecross=1,
               double base_tree=.95, double power_tree=2.0,
               double con_sd = 1,
               bool use_fscale = false,
               CharacterVector treef_name_="tsbtrees.txt",
               bool save_trees = false,
               bool silent_mode = false,
               bool monotone = false,
               bool incr = true,
               arma::vec round_ind = 0, double round_c2 = 0)
{

   // Begin output stream.
   std::string treef_name = as<std::string>(treef_name_);
   std::ofstream treef(treef_name.c_str());

   // If save_trees = false, set badbit so no output generated.
   if(save_trees==false){
      treef.setstate(std::ios_base::badbit);
   }

   // Random number generator, used in all draws.
   RNGScope scope;
   RNG gen;

   if(silent_mode==false){ Rcout << "\n*****Into bart main\n"; }

   //*****************************************************************************
   // Sample sizes.
   //*****************************************************************************

   // All data, is and oos.
   size_t n = y.size();
   size_t np = tpred.size();


   //*****************************************************************************
   // Read, format info about times
   //*****************************************************************************

   vec tref = unique(join_cols(tgt,tpred));  // Vector of unique time points, sorted in asc. order.
   int tlen = tref.n_elem;                   // Number of unique time points (T).

   if(silent_mode==false){
      Rcout << "unique times: " << endl << tref << endl;
      Rcout << "Number of total time points T: " << tlen << endl;
   }

   //*****************************************************************************
   // Read, format y
   //*****************************************************************************

   //------------------------------------------------------
   // Calculate alpha_t (ybar_t) vector, for conditional-on-time mean-centering data.
   //------------------------------------------------------

   vec alpha_t = zeros(tlen);    // Holds ybar at each t.
   vec n_vec = zeros(tlen);      // Holds sample sizes at each t.
   uvec idx;                     // Holds indices where ti = t, for calculating ybar_t.

   // Iterate through time points.
   for(int i = 0; i < tlen; i++){

      // All observations.
      idx = find(tgt == tref(i));
      n_vec(i) = idx.size();

      if(n_vec(i) > 0){
         alpha_t(i) = sum(y(idx)) / n_vec(i);
      }
   }

   //------------------------------------------------------
   // Mean-center data, conditional on the mean for each time point.
   //------------------------------------------------------

   // Create n-length vector with corresponding ybar_i for each y_i.
   // Will let us subtract off time-means and add back later.
   vec alpha_ti(y.size());
   idx.reset();

   for(size_t i = 0; i < y.size(); i++){
      idx = find(tref == tgt(i));
      alpha_ti(i) = as_scalar(alpha_t(idx));
   }

   // Subtract off means (in-place subtraction).
   y = y - alpha_ti;


   //--------------------------------------------------------------------------------------
   // Sufficient statistics for y - scalar and vector versions.

   //---------------------------------------------
   // Scalar versions.
   //---------------------------------------------

   sinfo allys;        //sufficient stats for all of y, use to initialize the bart trees.

   allys.n = n;
   allys.sy = sum(y);
   allys.sy2 = sum(y % y);

   //---------------------------------------------
   // Vector updates (T-length time vector).
   //---------------------------------------------

   allys.n_vec = n_vec;
   allys.sy_vec = zeros(tlen);

   // Stores vectors of sample means.
   vec ybar_vec(tlen);

   // Iterate through time points.
   for(int i = 0; i < tlen; i++){

      idx = find(tgt == tref(i));

      if(n_vec(i) > 0){
         allys.sy_vec(i) = sum(y(idx));
      }
   }

   ybar_vec = allys.sy_vec / allys.n_vec;

   //*****************************************************************************
   //* Read, format X, Xpred, x cutpoints for main prognostic covariates and trt covariates.
   //*****************************************************************************

   //---------------------------------------------
   // In-sample x
   //---------------------------------------------

   //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.

   // For all observations.
   size_t p = x.size()/n;           // Number of prognostic covariates

   //---------------------------------------------
   // Out-of-sample x, for predictions
   //---------------------------------------------

   // For all observations.
   vec xp = xpred;              // prognostic covariates

   // Check correct number of covariates.
   double npred = xp.size()/p;
   if(xp.size() != npred*p) Rcout << "error, wrong number of elements in prediction data set\n";

   //---------------------------------------------
   //x cutpoints
   //---------------------------------------------

   xinfo xi;

   // prognostic
   xi.resize(p);
   for(size_t i=0; i<p; ++i) {
      NumericVector tmp = xinfo_list[i];
      std::vector<double> tmp2;
      for(int j=0; j<tmp.size(); ++j) {
         tmp2.push_back(tmp[j]);
      }
      xi[i] = tmp2;
   }

   //*****************************************************************************
   // Setup for MCMC
   //*****************************************************************************

   //--------------------------------------------------
   //trees
   std::vector<tree> t(ntree);
   for(int i=0;i<ntree;i++) t[i].setm(ybar_vec/ntree); //if you sum the fit over the trees you get the fit.

   //--------------------------------------------------
   //prior and mcmc
   pinfo pi(tlen);

   pi.pbd = 1.0;      //prob of birth/death move
   pi.pb = .5;        //prob of birth given  birth/death

   pi.alpha = base_tree;    //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi.beta = power_tree;     //2 for bart means it is harder to build big trees.

   // Use override for sighat if provided by user.  Else use sd(y).
   if(sighat!=-999){
      pi.sigma = sighat;
   } else{
      pi.sigma = stddev(y);
   }

   // Use override for lambda if provided by user.  Else calculate using sighat and sigq.
   double qchi = 0;
   if(lambda==-999){
      qchi = R::qchisq(1-sigq, nu, true, false);
      lambda = (sighat * sighat * qchi) / nu;
   }

   // Initialize new time-related prior components.
   pi.mu0 = zeros(tlen);
   pi.ls = tlen / (PI * ecross);

   // If use_fscale=false, then var(f(x,t)) = con_sd^2 / m is variance of Cov matrix.
   // If use_fscale=true, then pi.var = 1/m, and C+ prior is induced with median con_sd via eta; sd of f(x,t) is con_sd*eta.
   pi.var = 1 / static_cast<double>(ntree);
   if(use_fscale==false) pi.var = pi.var * con_sd * con_sd;

   pi.Sigma0 = cov_se(tref, tref, pi.ls, pi.var); // Prior covariance matrix.
   pi.Prec0 = pi.Sigma0.i(); // Prior precision matrix.
   pi.eta = 1;
   pi.gamma = 1;

   //--------------------------------------------------
   // dinfo for in-sample

   // Initialize vector all fit; holds sum of fit of all trees for each observation.
   double* allfit = new double[n];

   // Initialize allfit to the ybar_t for each obs' corresponding time t value.
   for (size_t i=0; i<n; i++){
      idx.reset();
      idx = find(tref==tgt(i));
      allfit[i] = as_scalar(ybar_vec(idx));
   }

   // Initialize vectors to hold "residuals" (sum of all fits except tree j).
   // This is equal to: y - (allfit - temp)
   double* r = new double[n];

   // Initialize vector to hold fit of current tree j.
   double* ftemp = new double[n];

   // Set up data info.
   dinfo di;
   di.n=n; di.p=p; di.x = &x[0];
   di.y=r; //the "y" for each draw will be the residual

   // Set up data info time-point stuff.
   di.tlen = tlen;
   di.tref = tref;
   di.t = &tgt[0];

   // TROUBLESHOOTING:
   //  Rcout << "Completed di section" << endl;
   //  Rcout << "Checking pointer: di.t" << endl << di.t;
   //  Rcout << "Checking pointer: di.t value" << endl << di.t[0];

   //--------------------------------------------------
   // dinfo setup for out of sample

   dinfo dip; //data information for prediction
   dip.n = 0;

   if(np>0) {
      dip.n=np; dip.p=p; dip.x = &xp[0]; dip.y=0; //there are no y's!
      dip.tref = unique(tpred);  // Vector of unique time points, sorted in asc. order.
      dip.tlen = dip.tref.n_elem; // Number of unique time points (T).
      dip.t = &tpred[0];
   }

   //--------------------------------------------------
   //storage for ouput
   //--------------------------------------------------

   //for sigma draw
   double rss, restemp;
   NumericVector ssigma(nsim);

   //--------------------------------------------------
   // For eta draw
   // seta_sfy holds sum of (fit*y) and seta_sf2 holds sum of (fit^2)
   NumericVector seta(nsim);
   double seta_sfy; double seta_sf2;
   double seta_mean; double seta_var;

   //--------------------------------------------------
   // For gamma draw
   NumericVector sgamma(nsim);

   //--------------------------------------------------
   // For holding mcmc draws for in-sample and out-of-sample data.
   NumericMatrix sfit(nsim,n);
   NumericMatrix spred(nsim,dip.n);

   // For holding MH alphas for trees.
   NumericMatrix alpha (nsim+nburn,ntree);

   int thin = 1;

   //--------------------------------------------------
   //save stuff to tree file
   treef << xi << endl; //cutpoints
   treef << ntree << endl;  //number of trees
   treef << p << endl;  //dimension of x's
   treef << (int)(nsim/thin) << endl;

   //------------------------------------------------------
   // For rounding:
   // Save vector of indices for observations which are potentially rounded.
   // (For faster looping inside MCMC.)
   uvec round_idx;
   round_idx.resize(0);

   if(round_ind.size()>1){
      round_idx = find(round_ind==1);
   }

   //*****************************************************************************
   //* MCMC
   //*****************************************************************************

   if(silent_mode==false){ Rcout << "\nMCMC:\n"; }
   time_t tp;
   int time1 = time(&tp);

   //-------------------------------------------------
   // Loop through (nsim+nburn) MCMC iterations.
   //-------------------------------------------------
   for(int i=0;i<(nsim+nburn);i++) {

      // Progress report.
      if(silent_mode==false){ if(i%50==0) cout << "MCMC ITERATION: " << i << " of " << nsim+nburn <<  endl; }

      //-------------------------------------------------
      //draw trees
      //-------------------------------------------------

      for(int j=0;j<ntree;j++) {

         // In-sample tree fit.
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) {
            if(ftemp[k] != ftemp[k]) {
               Rcout << t[j] << endl;
               stop("nan in ftemp");
            }
            allfit[k] = allfit[k]-ftemp[k]*pi.eta; // allfits carrying correct scale: eta * treesum.
            r[k] = (y[k]-allfit[k])/pi.eta;
         }

         //birth/death
         alpha(i,j) = bd(t[j],xi,di,pi,gen);

         // draw mu's
         drmu(t[j],xi,di,pi,gen);

         // calculate fit, update allfit
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) allfit[k] += pi.eta*ftemp[k];

      } // End tree loop.

      //-------------------------------------------------
      //draw sigma.  Note: yhat = f(x) * eta.
      rss = 0.0;
      for(size_t k=0;k<n;k++) {
         restemp = y[k] - allfit[k];
         rss += restemp*restemp;
      }

      sighat = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));

      //-------------------------------------------------
      // draw eta

      if(use_fscale==true){

         seta_sfy = 0.0; seta_sf2 = 0.0;
         for(size_t k=0;k<n;k++) {
            seta_sfy += (allfit[k]/pi.eta) * y[k];
            seta_sf2 += (allfit[k]/pi.eta) * (allfit[k]/pi.eta);
         }

         // Calculate variance and mean for eta full conditional.
         seta_var = 1 / (1 / (pi.gamma*pi.gamma) + seta_sf2 / (sighat*sighat));
         //seta_mean = seta_var * (seta_sfy/(sighat*sighat));
         seta_mean = seta_var * (con_sd/(pi.gamma*pi.gamma) + seta_sfy/(pi.eta*sighat*sighat));

         // Draw normal full conditional for eta.
         double eta_old = pi.eta; // Save previous eta before drawing new one, for adjusting scaling.
         pi.eta = seta_mean + gen.normal(0.,1.) * sqrt(seta_var);

         // Update fits to have new pi.eta scaling.
         for(size_t k=0; k<n; ++k) {
            allfit[k] = allfit[k] * pi.eta / eta_old;
         }

      } else{
         pi.eta = 1;
      }

      // Update variance for backfitting with new sigma and eta draws: var is sig^2_y / pi.eta^2.
      pi.sigma = sighat/fabs(pi.eta);

      //-------------------------------------------------
      // draw gamma.  with eta ~ N(con_sd, gamma^2)
      double gamma_degf = 1;  // For C+(con_sd) prior.
      double gamma_a = (gamma_degf + 1) / 2;
      double gamma_b = (gamma_degf*con_sd*con_sd + pi.eta*pi.eta) / 2;

      pi.gamma = gen.gamma(gamma_a, gamma_b);

      //pi.gamma = sqrt((1 + ((pi.eta-con_sd) * (pi.eta-con_sd))) / gen.chi_square(2));

      //-------------------------------------------------
      // If some y's are rounded, draw updates for rounded y's.
      if(round_idx.size()>0){

        // Rcout << "drawing y's" << endl;

         for(int rnd = 0; rnd < round_idx.size(); rnd++){

            // Set up bounds for truncated normal.
            uword temp = round_idx[rnd];
            double temp_mean = allfit[temp];
            double temp_lb = temp_mean - round_c2;
            double temp_ub = temp_mean + round_c2;

            //cout << "Draw y truncated norm info for " << rnd << " of " << round_idx.size() << endl;
            //cout << "Mean: " << temp_mean << ", sd: " << pi.sigma << endl;
            //cout << "Lb: " << temp_lb << ", Ub: " << temp_ub << endl;

            // Draw the potentially-rounded y from a truncated normal.
            //cout << "y before: " << y[temp] << endl;
            y[temp] = rtnorm(temp_lb, temp_ub, temp_mean, pi.sigma);
            //scout << "y after: " << y[temp] << endl;
         }
      }

      //-------------------------------------------------
      // Save MCMC output.

      // Save only if exceeded number of iterations to nburn.
      if(i>=nburn) {

         // Save the ntree tree outputs.
         for(int j=0;j<ntree;j++) treef << t[j] << endl;

         // Save vector of sigma draws; (i-nburn) fills from the beginning of the vector.
         ssigma(i-nburn) = sighat;

         // Save vectors of eta and gamma draws.
         seta(i-nburn) = pi.eta;
         sgamma(i-nburn) = pi.gamma;

         //------------------------------
         // In-sample fits for matrices
         //------------------------------

         // Loop through trt=1 obs first.
         for(size_t k=0;k<n;k++) {

            if(monotone==false){
               // sfit saves a matrix of MCMC iterations (rows) for each obs (cols).
               sfit(i-nburn, k) = allfit[k] + alpha_ti[k]; // Allfit already carries mu scale. Just add back in overall means.
            } else{
               sfit(i-nburn, k) = fit_i_pava(k, t, xi, di, alpha_t, pi.eta, incr);
            }

         }

         //------------------------------
         // Matrix of out of sample fit values.
         //------------------------------

         if(dip.n){

            // Loop through trt=1 obs first.
            for(size_t k=0;k<dip.n;k++) {

               if(monotone==false){
                  spred(i-nburn, k) = fit_i(k, t, xi, dip, alpha_t, pi.eta);
               } else{
                  spred(i-nburn, k) = fit_i_pava(k, t, xi, dip, alpha_t, pi.eta, incr);
               }

            }

         } // End out of sample fit values.

      } // End saving MCMC outputs.

   } // End MCMC loop.

   //-------------------------------------------------
   // MCMC loop time keeping
   //-------------------------------------------------

   if(silent_mode==false){   Rcout << endl << "MCMC LOOP FINISHED SUCCESSFULLY" << endl; }
   int time2 = time(&tp);
   if(silent_mode==false){ Rcout << "time for loop: " << time2 - time1 << endl; }

   //*****************************************************************************
   // Rescale y vector.
   //*****************************************************************************
   y = y + alpha_ti;




   // //*****************************************************************************
   // //* Posterior predictive draws for each MCMC draw: N(f-hat, sighat)
   // //*****************************************************************************
   //
   // mat pp(nsim,n);
   // mat pp_pred(nsim,np);
   //
   // if(silent_mode==false){ Rcout << endl << "Begin posterior predictive draws." << endl; }
   //
   // // In-sample.  Do one row at a time; each row is an MCMC iteration.
   // for(int i=0; i<nsim; i++){
   //
   //    // Progress report.
   //    if(silent_mode==false){  if(i%50==0) cout << "Posterior predictive draw: " << i << " of " << nsim <<  endl; }
   //
   //    // In-sample.
   //    for(size_t j=0; j<n; j++){
   //       pp(i,j) = randn() * ssigma(i) + sfit(i,j);
   //    }
   //
   //    // Out-of-sample.
   //    for(size_t j=0; j<np; j++){
   //       pp_pred(i,j) = randn() * ssigma(i) + spred(i,j);
   //    }
   //
   // }
   //
   // if(silent_mode==false){ Rcout << endl << "Posterior predictive draws complete." << endl; }

   //*****************************************************************************
   //* Clean up and return results.
   //*****************************************************************************

   t.clear();
   delete[] allfit;
   delete[] r;
   delete[] ftemp;
   treef.close();

   // If save_tree = false, delete tree file.
   if(save_trees==false){
      std::remove(treef_name.c_str());
   }

   return(List::create(_["mcmcdraws"]      = sfit,     // Matrix of all mcmc draws for in-sample.  Rows = MCMC iter, cols = obs.
                       _["mcmcdraws_oos"]  = spred,    // Matrix of mcmc draws for out of sample.
                       _["sigma"]          = ssigma,   // Vector of sigma MCMC draws.
                       _["eta"]            = seta,     // Vector of eta MCMC draws.
                       _["gamma"]          = sgamma,   // Vector of gamma MCMC draws.
                       _["alpha"]          = alpha,    // matrix of MH alpha's for trees
                       _["y"]              = y,         // Vector of y values (may be updated with rounding draws),
                       _["round_idx"]      = round_idx
   ));
}
