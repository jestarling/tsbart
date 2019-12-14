# Fit tsBART.  Wrapper function for calling tsbartFit.cpp and tsbartProbit.cpp.


### For argument validation.
.ident <- function(...){
   # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
   args <- c(...)
   if( length( args ) > 2L ){
      #  recursively call ident()
      out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
   }else{
      out <- identical( args[1] , args[2] )
   }
   return( all( out ) )
}

### Rounds to nearest value.
round_to <- function(x, b) {
   round(x/b)*b
}

### For testing equality.
test_eq <- function(a,b){
   abs(a - b) < .001
}

tsbart <- function(y, tgt, x, tpred=NULL, xpred=NULL,
                   nburn=100, nsim=1000, ntree=200,
                   lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                   ecross=1, base_tree=.95, power_tree=2, sd_control = 2*sd(y),
                   use_fscale=TRUE,
                   probit=FALSE, yobs=NULL, verbose=T, mh=F, save_inputs=T,
                   monotone="no", binsize=NULL){

   ################################################################
   # Capture key arguments.
   ################################################################
   inputs = cbind.data.frame(
      'arg' = c('nburn',
                'nsim',
                'ntree',
                'lambda',
                'sigq',
                'sighat',
                'nu',
                'ecross',
                'base_tree',
                'power_tree',
                'sd_control',
                'use_fscale',
                'probit',
                'verbose',
                'mh',
                'monotone',
                'binsize'),
      'value' = c(ifelse(is.null(nburn),"NULL",nburn),
                  ifelse(is.null(nsim),"NULL",nsim),
                  ifelse(is.null(ntree),"NULL",ntree),
                  ifelse(is.null(lambda),"NULL",lambda),
                  ifelse(is.null(sigq),"NULL",sigq),
                  ifelse(is.null(sighat),"NULL",sighat),
                  ifelse(is.null(nu),"NULL",nu),
                  ifelse(is.null(ecross),"NULL",ecross),
                  ifelse(is.null(base_tree),"NULL",base_tree),
                  ifelse(is.null(power_tree),"NULL",power_tree),
                  ifelse(is.null(sd_control),"NULL",sd_control),
                  ifelse(is.null(use_fscale),"NULL",use_fscale),
                  ifelse(is.null(probit),"NULL",probit),
                  ifelse(is.null(verbose),"NULL",verbose),
                  ifelse(is.null(mh),"NULL",mh),
                  ifelse(is.null(monotone),"NULL",monotone),
                  ifelse(is.null(binsize),"NULL",binsize))
   )

   ################################################################
   # Validate inputs.
   ################################################################

   #---------------------------------------------------------------
   # If not predicting, set pred objects to first three obs.
   # These are not output.
   #---------------------------------------------------------------
   predict = 1

   if((is.null(tpred) & !is.null(xpred)) | (!is.null(tpred) & is.null(xpred))){
      stop('To predict, provide both tpred and xpred.')
   }
   if( is.null(tpred) & is.null(xpred) ){
      predict = 0
      tpred = tgt[1:min(3,length(tgt))]
      xpred = x[1:min(3,nrow(x)),,drop=F]
   }

   #---------------------------------------------------------------
   # Data size.
   #---------------------------------------------------------------

   # Check data size match.
   if( !.ident(length(y), length(tgt), nrow(x))){

      stop("Data size mismatch. The following should all be equal:
           length(y): ", length(y), "\n",
           "length(tgt): ", length(tgt), "\n",
           "nrow(x): ", nrow(x), "\n")
   }

   # Check out-of-sample data size match.
   if( !.ident(length(tpred), nrow(xpred))){

      stop("Data size mismatch. The following should all be equal:
           length(tpred): ", length(tpred), "\n",
           "nrow(xpred): ", nrow(xpred), "\n")
   }

   #---------------------------------------------------------------
   # Probit checks.
   #---------------------------------------------------------------
   if(probit==TRUE){

      # Data size match including yobs.
      if( !.ident(length(y), length(yobs), length(tgt), nrow(x))){
         stop("Data size mismatch. The following should all be equal:
              length(y): ", length(y), "\n",
              "length(yobs): ", length(yobs), "\n",
              "length(tgt): ", length(tgt), "\n",
              "nrow(x): ", nrow(x), "\n")
      }

      #Yobs must be only 0/1.  Y must not be only 0/1.
      if(length(unique(y))>2) warning("In probit case,
                                      y should contain initial latent variables,
                                      and yobs should contain observed binary values.")

      if(length(unique(yobs))>2) stop("In probit case,
                                      y should contain initial latent variables,
                                      and yobs should contain observed binary values.")
      if(is.null(yobs)) stop("yobs must be populated when probit=TRUE, and should contain observed binary values of 0/1.")

      # Warn user that manually input lambda and sighat values are ignored in probit case.
      if(!is.null(lambda) || !is.null(sighat)) warning("lambda and sighat inputs are ignored in probit case, as prior
                                                       tuning for sigma^2 is not applicable.")
   }

   if(probit==FALSE){
      if(!is.null(yobs)) stop("yobs is only for probit=TRUE case. Must be NULL when probit=FALSE.")
   }

   #---------------------------------------------------------------
   # Other inputs.
   #---------------------------------------------------------------

   if(any(is.na(y))) stop("Missing values in y")
   if(any(is.na(tgt))) stop("Missing values in tgt")
   if(any(is.na(tpred))) stop("Missing values in tpred")
   if(any(is.na(x))) stop("Missing values in x")
   if(any(is.na(xpred))) stop("Missing values in xpred")

   if(length(unique(y))<5 && probit==FALSE) warning("y appears to be discrete")
         # Ok for probit bc might be initializing 0/1 to two latent values.

   if(nburn<0) stop("nburn must be positive")
   if(nsim<0) stop("nsim must be positive")
   # if(nburn<100) warning("A low (<100) value for nburn was supplied")
   if(ecross<=0) stop("ecross must be positive")
   if(class(ecross)=="character" & ecross!="tune") stop("ecross must be a positive value or set to 'tune', case-sensitive.")
   #if(save_trees==T & is.null(tree_fname)) stop('If save_trees=TRUE, must spply a filename.')

   #---------------------------------------------------------------
   # Monotonicity and heaping bias for potentially rounded y's.
   #---------------------------------------------------------------

   # monotone must be NULL, or "incr", or "decr".
   if(!monotone %in% c('no','incr','decr')) stop("monotone must be 'no', 'incr', or 'decr'.")

   if(!is.null(binsize)){
      if((class(binsize)!="numeric") | (binsize<=0)) stop('binsize must be a positive numeric scalar, ex: 10,50, or 100.')
      if(length(binsize)>1) stop('binsize must be a positive numeric scalar, ex: 10,50, or 100.')
   }

   ################################################################
   # Order data frame in ascending order of t value.
   # (For monotonicity projection.)
   ################################################################
   perm = order(tgt, decreasing=FALSE)
   perm_oos = order(tpred, decreasing=FALSE)

   ################################################################
   # Calculate rounding half-bound-width and rounding indicators for binsize heaping
   # bias.
   ################################################################

   # Initialize rounding vectors.
   rounding_ind = 0
   rounding_halfWidth = 0

   # If binsize is provided, populate rounding details.
   if(!is.null(binsize)){

      # Vector of indicators for which observations may be rounded, based on cutoff binsize.
      rounding_ind = ifelse( test_eq( round_to(y,binsize), y), 1, 0)

      print(paste0('y potentially rounded to nearest ', binsize,'.'))
      print(paste0(sum(rounding_ind), ' potentially rounded observations (', round(sum(rounding_ind)/length(rounding_ind),3),' of obs)'))

      # Scalar c/2 bound half-widths for each observation. (Same for all obs, such as 100g/2.)
      rounding_halfWidth = binsize/2
   }

   ################################################################
   # Create model matrix and set up hyperparameters.
   ################################################################

   # Scale/center y.
   ybar = mean(y)
   ysd = sd(y)
   y = (y - ybar) / ysd

   # Scale rounding bound half-width c/2.
   if(!is.null(binsize)){
      rounding_halfWidth = rounding_halfWidth / ysd
   }

   # Model matrix.
   xx = tsbart::makeModelMatrix(x)
   xxpred = tsbart::makeModelMatrix(xpred)
   cutpoints = tsbart::makeCutpoints(xx)

   # Sighat and lambda calibration.
   if(is.null(sighat)){
      df = cbind.data.frame(y, tgt, xx)
      lmf = lm(y ~ ., data=df)
      sighat = sigma(lmf)
   }

   if(is.null(lambda)){
      qchi = qchisq(1-sigq, nu)
      lambda = (sighat * sighat * qchi) / nu
   }

   ################################################################
   # Set up probit parameters if probit=T.
   ################################################################
   offset=0

   if(probit==TRUE){
      phat = mean(unlist(yobs))
      offset = qnorm(phat)
   }

   ################################################################
   # Parameter tuning if necessary.
   ################################################################
   if(ecross=="tune"){
      tuned = tuneEcross(ecross_candidates = seq(.25,5,by=.25),
                         y[perm], tgt[perm], tpred[perm_oos], x[perm,], xpred[perm_oos,],
                         nburn=500, nsim=500, ntree=200,
                         lambda, sigq, sighat, nu,
                         base_tree, power_tree,
                         probit, yobs[perm], monotone=monotone)
      ecross = tuned$ecross_opt
   }


   ################################################################
   # Call tsbartFit.cpp or tsbartProbit.cpp if monotone=NULL,
   # or tsbartFitMono.cpp or tsbartProbitMono.cpp if monotone= "incr" or "decr"
   ################################################################
   out = NULL

   # Continuous response.
   if(probit==FALSE){
      out = tsbartFit(y=y[perm], tgt=tgt[perm], tpred=tpred[perm_oos],
                      x=t(xx[perm,]), xpred=t(xxpred[perm_oos,]), xinfo_list=cutpoints,
                      nburn=nburn, nsim=nsim, ntree=ntree,
                      lambda=lambda, sigq=sigq, sighat=sighat, nu=nu,
                      ecross=ecross, base_tree=base_tree, power_tree=power_tree,
                      treef_name_="tsb_trees.txt",
                      con_sd = ifelse(abs(2*ysd - sd_control)<1e-6, 2, sd_control/ysd),
                      use_fscale=use_fscale,
                      save_trees=FALSE,
                      silent_mode=!verbose,
                      monotone=ifelse(monotone=="no", FALSE, TRUE),
                      incr=ifelse(monotone=="incr", TRUE, FALSE),
                      round_ind=rounding_ind, round_c2=rounding_halfWidth)

   # Probit response, not monotone.
   } else{
      out = tsbartProbit(y=y[perm], yobs=yobs[perm], tgt=tgt[perm], tpred=tpred[perm_oos],
                         x=t(xx[perm,]), xpred=t(xxpred[perm_oos,]), xinfo_list=cutpoints,
                         nburn=nburn, nsim=nsim, offset=offset, ntree=ntree, ecross=ecross,
                         base_tree=base_tree, power_tree=power_tree,
                         treef_name_="tsb_trees.txt",
                         save_trees=FALSE,
                         silent_mode=!verbose)
   }

   ################################################################
   # Adjust alpha's and add accept/reject indicator.
   # Note: bd function returns:
   #     alpha in (0,1) for births which are accepted
   #     -alpha in (-1,0) for deaths which are accepted
   #     10+alpha for rejected births
   #     -alpha-10 for rejected deaths.
   ################################################################

   if(mh){
      # Con metropolis info.
      bd = ifelse(out$alpha<=0,0,1)             # 1 = birth, 0 = death
      accept = ifelse(abs(out$alpha)<10,1,0)   # 1 = accepted, 0 = rejected
      alpha = ifelse(accept==1, abs(out$alpha), abs(out$alpha)-10)

      # Assemble dataframes and convert bd to character.
      metrop = cbind.data.frame(
         'iter' = rep(1:(nburn+nsim), times=ntree),
         'tree' = rep(1:ntree, each=nburn+nsim),
         'accept' = as.numeric(accept),
         'alpha' = as.numeric(alpha),
         'bd' = as.numeric(bd)
      )

      metrop$bd = ifelse(metrop$bd==1,'birth','death')
   }

   ################################################################
   # Format output.
   ################################################################

   if(predict){
      output = list('mcmcdraws' = out$mcmcdraws[,order(perm)]*ysd + ybar,
                 'mcmcdraws_oos' = out$mcmcdraws_oos[,order(perm_oos)]*ysd + ybar,
                 'sigma' = out$sigma*ysd,
                 'treefit_sd' =  abs(sd_control * out$eta), # sd of BART fit f(x,t).
                 'eta' = out$eta,
                 'gamma' = out$gamma,
                 'ecross' = ecross)

   } else{
      output = list('mcmcdraws' = out$mcmcdraws[,order(perm)]*ysd + ybar,
                 'sigma' = out$sigma*ysd,
                 'treefit_sd' =  abs(sd_control * out$eta), # sd of BART fit f(x,t).
                 'eta' = out$eta,
                 'gamma' = out$gamma,
                 'ecross' = ecross)
   }

   # Include rounding info if provided.
   if(!is.null(binsize)){
      output$y_potentially_rounded_indicator = rounding_ind[order(perm)]
      output$y_imputed = as.numeric(out$y[order(perm)]) * ysd + ybar
   }

   # Include metropolis info if indicated.
   if(mh){
      output$metrop = metrop
   }

   # Include inputs if indicated.
   if(save_inputs){
      output$inputs = inputs
   }

   output$round_idx = out$round_idx
   # Return output.
   return(output)
   }
