# Final Project: Regression
# Patrick Grennan
# Jason Schapiro

require(lars)
require(lasso2)

main <- function(){
  load("baa.ratios.rda")
  
  # we can access the rows in the ratios matrix by calling ratios[tf[i],]
  tfNames <- scan(file="justTFS.txt", what="character", sep="\n")
  transFactors <- tfList(tfNames)
  
  # Since we don't have the clusters, this function is a placeholder
  clusters <- getClusters()
  
  # List to store the results
  results <- list()
  
  # this will probably be changed to length(clusters)
  for (i in 1:NROW(clusters)){
    # this should be ONE vector, because it should be the mean of the cluster
    y <- as.vector(clusters[i,])
    
    #str(transFactors)
    
    # Here we could also look for the best correlated TF's too
    x <- predictors(y,transFactors)
    
    # we have to put in the TRANSPOSE of x
    # Do we want this to be tmp.l1ce or the final model?
    model <- bestFit(t(x),y)
    
    # add model to list of models for clusters
    results[[i]] <- model
    #str(model)
    cat("SUCCESS\n")
  }
  
  # TODO: Save file into .rda for viz group?
  str(results)
  invisible()
  
  
  
}

getClusters <- function(){
  load("baa.ratios.rda")
  clusters <- as.vector(ratios[1,])
  clusters <- rbind(clusters,as.vector(ratios[1,]))
  invisible(clusters)
  
}


# A function that loads the data
tfList <- function(tfNames){
  load("baa.ratios.rda")
  miniList <- list()
  
  # This loop works
  for (i in 1:length(tfNames)){
    # cat(i, "\n")
    # remove as.vector for the names
    miniList[[i]] <- as.vector(ratios[tfNames[i],])
  }
  
  #for (q in 1:length(z[[1]])){
  #  rowTotal <- 0
  #  for (l in 1:length(z)){
  #    rowTotal = rowTotal+z[[l]][q]
  #  }
  #  cat("Mean is: ", (rowTotal/length(z)), "\n")
  #}
   
  # to return just the list
  #invisible(miniList)
  invisible(do.call(rbind,miniList))
  
}

  
predictors <- function(y, tfs) {
  cors <- double()
	for (i in 1:dim(tfs)[1]) {
		t = tfs[i,]
		cors <- append(cors, abs(cor(y, t, use="na.or.complete")))
	}
	
	invisible(tfs[sort(cors, decreasing = T, index.return = T)$ix[1:20],])
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
   
  # Return this?
  l1ce.final  <- l1ce( y ~ x , sweep.out = ~ 1, standardize = TRUE,
                        bound = foundT, absolute.t = FALSE)
  predictors <- which( abs(l1ce.final$coefficients) > 0)
  predictors <- predictors[predictors != 1] - 1
      
  cat("The number of predictors is: ", length(predictors), "\n")

  if ( length(predictors) < dim(x)[2] ) { 
      x <- as.matrix( x[,predictors] )
  }
  
  # Refitting the Model with best parameters from l1ce()
  lm.final <- lm(y ~ x)
  plot( y, predict( lm.final ) )
  abline(0,1, col = 2, lwd = 3, lty = 2)
  #summary(lm.final)
  
  invisible(l1ce.final)
  
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