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

tsbart <- function(y, tgt, tpred, x, xpred, nburn, nsim, ntree=200,
                   lambda=NULL, sigq=.9, sighat=NULL, nu=3,
                   ecross=1, base_tree=.95, power_tree=2, sd_control = 2*sd(y),
                   use_fscale=TRUE,
                   probit=FALSE, yobs=NULL, verbose=T){

   ################################################################
   # Validate inputs.
   ################################################################

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
      if(length(unique(y))<3) warning("y appears to be discrete.  In probit case,
                                      y should contain initial latent variables, and yobs should contain observed binary values.")

      if(length(unique(yobs))>2) stop("yobs appears to be discrete.  In probit case,
                                      y should contain initial latent variables, and yobs should contain observed binary values.")
      if(is.null(yobs)) stop("yobs must be populated when probit=TRUE, and should contain observed binary values.")

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

   if(length(unique(y))<5) warning("y appears to be discrete")

   if(nburn<0) stop("nburn must be positive")
   if(nsim<0) stop("nsim must be positive")
   # if(nburn<100) warning("A low (<100) value for nburn was supplied")
   if(ecross<=0) stop("ecross must be positive")
   if(class(ecross)=="character" & ecross!="tune") stop("ecross must be a positive value or set to 'tune', case-sensitive.")
   #if(save_trees==T & is.null(tree_fname)) stop('If save_trees=TRUE, must spply a filename.')

   ################################################################
   # Create model matrix and set up hyperparameters.
   ################################################################

   ybar = mean(y)
   ysd = sd(y)

   y = (y - ybar) / ysd

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
      phat = mean(unlist(y))
      offset = qnorm(phat)
   }

   ################################################################
   # Parameter tuning if necessary.
   ################################################################
   if(ecross=="tune"){
      tuned = tuneEcross(ecross_candidates = seq(.25,5,by=.25),
                         y, tgt, tpred, x, xpred, nburn=500, nsim=500, ntree=200,
                         lambda, sigq, sighat, nu,
                         base_tree, power_tree,
                         probit, yobs)
      ecross = tuned$ecross_opt
   }


   ################################################################
   # Call tsbartFit.cpp or tsbartProbit.cpp
   ################################################################
   out = NULL

   if(probit==FALSE){
      out = tsbartFit(y=y, tgt=tgt, tpred=tpred, x=t(xx), xpred=t(xxpred), xinfo_list=cutpoints,
                      nburn=nburn, nsim=nsim, ntree=ntree,
                      lambda=lambda, sigq=sigq, sighat=sighat, nu=nu,
                      ecross=ecross, base_tree=base_tree, power_tree=power_tree,
                      treef_name_="tsb_trees.txt",
                      con_sd = ifelse(abs(2*ysd - sd_control)<1e-6, 2, sd_control/ysd),
                      use_fscale=use_fscale,
                      save_trees=FALSE,
                      silent_mode=!verbose)
   } else{
      out = tsbartProbit(y=y, yobs=yobs, tgt=tgt, tpred=tpred, x=t(xx), xpred=t(xxpred), xinfo_list=cutpoints,
                         nburn=nburn, nsim=nsim, offset=offset, ntree=ntree, ecross=ecross,
                         base_tree=base_tree, power_tree=power_tree,
                         treef_name_="tsb_trees.txt",
                         save_trees=FALSE,
                         silent_mode=!verbose)
   }

   ################################################################
   # Format output.
   ################################################################

   return(list('mcmcdraws' = out$mcmcdraws*ysd + ybar,
               'mcmcdraws_oos' = out$mcmcdraws_oos*ysd + ybar,
               'sigma' = out$sigma*ysd,
               'eta' = out$eta,
               'gamma' = out$gamma,
               'ecross' = ecross))
   }
