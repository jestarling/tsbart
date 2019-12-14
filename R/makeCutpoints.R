# Creates cutpoint list.

makeCutpoints = function(X, gridlen = 10000){
   #---------------------------------------------------
   # FUNCTION: Creates cutpoints for a design matrix X.
   #---------------------------------------------------

   # Initialize empty list.
   cutpoints = list()

   # Loop through coluns of design matrix X.
   for(j in 1:ncol(X)){

      if(class(X[,j]) %in% c("character","factor")){
         # If discrete (factors or characters), use this.
         cutpoints[[j]] = sort(unique(X[,j]))

      }else if(sum(unique(X[,j]) %in% c(0,1))==2 & length(unique(X[,j]))==2){
         # If categorical (0/1), use this.
         cutpoints[[j]] = c(0,1)

      } else{
         min = min(X[,j])
         max = max(X[,j])
         cutpoints[[j]] = seq(min,max,length.out=gridlen)
      }
   }

   return(cutpoints)
}
