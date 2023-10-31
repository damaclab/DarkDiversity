dark.pred.ca <-function( calib.dta, # Data frame (or matrix) with calibration dataset
                      pred.dta, # Data frame (or matrix) for which predictions are #made
                      # Following argument determines the cut-off threshold for predicted values
                      # This threshold is determined separately for each taxon (response variable)
                      # minpred: use smallest predicted value for real positive observations
                      # binminpred: as above, but used on data transformed to binary (0/1) form
                      # minobs: use smallest positive observed value
                      method = c("minpred","binminpred","minobs"),
                      var.expl = 0.6) # Minimum community variation explained by used #CA axes
{
  require( vegan)
  is.present <- function(x){ sum( x) > 1.0e-16 }
  n1 <- dim( calib.dta)[1]# Number of observations in calibration dataset
  n2 <- dim( pred.dta) [1]#Number of observations in dataset used for prediction
  nc <- dim( calib.dta)[2]# Number of columns (taxa) in calibration dataset
  
  # Append the prediction data frame to the end of calibration data frame
  merged.dta <- as.data.frame( rbind( as.matrix( calib.dta),
                                      matrix( 0, n2, nc)))
  
  # Make sure the column names match
  names.match <- match( names( pred.dta), names( merged.dta), 0)
  add.dta <- pred.dta[ , (names.match > 0)]
  add.cols <- names.match[ (names.match > 0)]
  merged.dta[(n1+1):(n1+n2), add.cols] <- add.dta
  
  
  # Adjust the names of appended rows
  dimnames( merged.dta)[[1]][(n1+1):(n1+n2)] <-
    paste("Pred", dimnames( add.dta)[[1]])
  # Omit columns which are empty in merged data
  merged.dta <- merged.dta[, sapply( merged.dta, is.present)]
  # Update count of columns (now in the combined data frame)
  nc <- dim( merged.dta)[2]
  # Calculate correspondence analysis, optionally binarise the data
  if(method=="binminpred")
  {
    bin.dta <- merged.dta;
    for(i in 1:nc) bin.dta[,i] <- as.numeric( (bin.dta[,i]>0));
    cca.1 <- cca( bin.dta)
  } else
    cca.1 <- cca( merged.dta)
  
  
  # Determine CA axes to use: first n.used.ax match required explained variation
  n.used.ax <- sum( cumsum(cca.1$CA$eig)/cca.1$tot.chi <= var.expl) + 1
  
  
  # Predict at the scale of original data using selected count of axes
  predicted.dta <- predict( cca.1, rank=n.used.ax, type="response")
  
  # Prepare data frame for predicted values (number of columns might differ)
  pred.pres <- as.data.frame( matrix( 0, n2, nc))
  # Copy over the matching row and column names
  dimnames(pred.pres)[[1]] <- dimnames(add.dta)[[1]]
  dimnames(pred.pres)[[2]] <- dimnames(merged.dta)[[2]]
  
  # Determine the cut-off thresholds for each column
  cut.off.vals <- rep( NA, nc)
  
  
  # Threshold determination for (bin)minpred method
  get.minpred <- function( real, pred)
  {
    pred.pos <- (pred > 0)
    real.pres <- (real > 0)
    cond.comb <- pred.pos & real.pres
    ifelse( sum(cond.comb) > 0, min( pred[cond.comb]), max( pred))
  }
  
  
  # Threshold determination for minobs method
  get.minobs <- function( real, pred)
  { min( real[real>0]) }
  
  for( i in 1:nc) # For each predicted column ...
  {
    spc.real <- merged.dta[,i]
    spc.pred <- predicted.dta[,i]# ... determine cut-off
    cut.off <- switch( method,
                       minpred = get.minpred( spc.real, spc.pred),
                       binminpred = get.minpred( spc.real, spc.pred),
                       minobs = get.minobs( spc.real, spc.pred),
                       NA)
    spc.pred <- spc.pred / cut.off
    spc.pred[spc.pred < 0.01] <- 0 # ... replace too low values with 0s
    predicted.dta[,i] <- spc.pred # copy back to predicted values ...
    spc.pred[spc.real > 0] <- 0 # Eliminate observed values in predicted cases
    pred.pres[,i] <- spc.pred[(n1+1):(n1+n2)]
    cut.off.vals[i] <- cut.off # Record cut-off value
  }
  spc.names <- dimnames(pred.pres)[[2]]
  row.names <- dimnames(pred.pres)[[1]]
  
  preds <- list()
  
  for( i in 1:n2) # For each case to be predicted, store in preds list the names
  { # and values of additional (dark-diversity) species
    x.case <- pred.pres[i,]
    x.ord <- rev(order(x.case))
    num.pos <- sum(x.case > 0)
    num.start <- ifelse( num.pos > 0, 1, 0)
    x.ord <- x.ord[num.start:num.pos]
    
    # Use name of predicted case as the name of 'preds' list' component (also list),
    # storing new (dark) species names and predicted values as two vectors8
    preds[[row.names[i]]] <- list( names = spc.names[x.ord],
                                   vals = as.numeric(x.case[x.ord]))
  }
  # Return value storing prediction results and info on adopted parameters of
  # the algorithm
  list( preds = preds, cut.offs = cut.off.vals, spc.names = spc.names,
        num.axes.tot = cca.1$CA$rank, num.axes.used = n.used.ax,
        method = method, n1=n1, n2=n2, pred.mat = predicted.dta[(n1+1):(n1+n2),])
}