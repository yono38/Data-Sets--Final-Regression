# Final Project Mockup
require(lars)
require(lasso2)

main <- function(){
  # TODO: read in TF's, can use
  # we can access the rows in the ratios matrix by calling ratios[tf[i],]
  tfs <- scan(file="justTFS.txt", what="character", sep="\n")
  
  # TODO: get clusters
  # I'm gonna assume this comes in a matrix of the form
  clusters <- 0
  
  # List to store the results
  results <- list()
  
  for (i in 1:length(clusters)){
    y <- mean(clusters[,i])
    
    # Here we could also look for the best correlated TF's too
    x <- tfs
    
    # Do we want this to be tmp.l1ce or the final model?
    model <- bestFit(x,y)
    
    # add model to list of models for clusters
    results[[i]] <- model
  }
  
  # TODO: Save file into .rda for viz group?
  
  invisible()
  
  
  
}



# From here down is homework 5

bestFit <- function( x, y, kFolds=5, stepSize = .05, printSteps = FALSE){
  #require(lars)
  #require(lasso2)
  
  # Plotting the CV curve
  #if (printSteps == FALSE){
  #  cv.lars(x, y, K=kFolds)
  #}
  
  # Converting the step size to a range
  convertToRange = 1/stepSize
  bounds = (1:convertToRange)/convertToRange
  
  cv <- numeric(convertToRange)
  cv.err <- numeric(convertToRange)
  
  # Arbitrary Large Numbers to track best values so far
  mincv <- 10000
  mincv.err <- 10000
  bestT <- 0
  bestNumPredictors <- 0
  
  # Stepping through every shrinkage parameter
  for (i in 1:length(bounds)){ 
    l1ce.example  <- l1ce( y ~ x , sweep.out = ~ 1, standardize = TRUE,
                        bound = bounds[i], absolute.t = FALSE)
    pFromL1 <- which( abs(l1ce.example$coefficients) > 0)
    pFromL1 <- pFromL1[pFromL1 != 1] - 1
    
    # Making a model from the predictors found with the L1 Shrinkage
    cv.obj <- cv.lm( y, x, k = kFolds, p = pFromL1 )
    cv[i] <- cv.obj$cv
    cv.err[i] <- cv.obj$cv.err
    
    if (printSteps == TRUE){
      cat("Shrinkage Parameter is:", bounds[i], "\n")
      cat("CV is:", cv.obj$cv, "\n")
      cat("CV.err is:", cv.obj$cv.err, "\n")
    }
    
    # Testing for best shrinkage parameter
    if(cv.obj$cv < mincv){
      mincv = cv.obj$cv
      mincv.err = cv.obj$cv.err
      bestT = bounds[i]
      bestNumPredictors = length(pFromL1)
      bestPredictors = pFromL1
    }
  }

  if (printSteps == TRUE){
    plot.cv.lm( 1:convertToRange , cv, cv.err )
    cat("Best CV:",mincv, "\n")
    cat("Best t:", bestT, "\n")
    cat("The Best number of Predictors was:", bestNumPredictors, "\n")
    cat("The predictors were:", "\n")
    print(bestPredictors)
  }
      
  # This is the first shrinkage parameter that's within the  best cv+cv.err
  foundT <- 0
      
  for (i in 1:length(bounds)){
    if(cv[i] < (mincv + mincv.err)){
      foundT <- i*.05
      break
    }
  }
  
  cat("The shrinkage parameter t is:", foundT, "\n")
      
  l1ce.final  <- l1ce( y ~ x , sweep.out = ~ 1, standardize = TRUE,
                        bound = foundT, absolute.t = FALSE)
  predictors <- which( abs(l1ce.final$coefficients) > 0)
  predictors <- predictors[predictors != 1] - 1
      
  cat("The number of predictors is: ", length(predictors))

  if ( length(predictors) < dim(x)[2] ) { 
      x <- as.matrix( x[,predictors] )
  }
  
  # Refitting the Model with best parameters from l1ce()
  lm.final <- lm(y ~ x)
  plot( y, predict( lm.final ) )
  abline(0,1, col = 2, lwd = 3, lty = 2)
  summary(lm.final)
  
}


      
##### HELPER FUNCTIONS #####     
cv.folds <- function(n, folds = 10) {
  split(sample(1:n), rep(1:folds, length = n))
}

cv.lm <- function( y, x, k= 5, p = 1:dim(x)[2], method = "rss" ) {   
   if ( length(p) < dim(x)[2] ) { 
      x <- as.matrix( x[,p] )
   }   
   cv.subsets <- cv.folds( length(y), folds = k)
   cv.rss <- numeric(k)
   for (i in 1:k) {
   	  tmp.lm <- lm( y[ - cv.subsets[[i]] ] ~ x[ - cv.subsets[[i]],  ] )
   	  y.hat <- predict.from.lm( tmp.lm, x[cv.subsets[[i]],] )                                 
   	  cv.rss[i] <- mean( (y[ cv.subsets[[i]] ] - y.hat )**2 )  
   }	  
   return( list( cv = mean( cv.rss ), cv.err = sqrt( var( cv.rss ) / k ) ) ) 
}	


plot.cv.lm <- function (x, cv, cv.err ) {
    plot(x, cv, type = "b", ylim = range(cv, cv + cv.err, cv - cv.err))
    error.bars(x, cv + cv.err, cv - cv.err, width = 1/(1.5* length(x)) )
    invisible() 
}

    
predict.from.lm <- function( lm1, x) {
     x <- as.matrix( x )
     if (class(lm1) != "lm") {
       stop("input class is not lm")
     	return( FALSE )
     } 
     coeff <- lm1$coefficients
     n <- dim( x )[1]
     tmp.mat <- t( cbind( rep(1,n) , x) ) * coeff
     y.hat <-  apply( tmp.mat, 2, sum)
     invisible( y.hat ) 
}