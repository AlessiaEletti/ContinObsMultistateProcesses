my.msboot = function (theta, data, B = 5, id = "id", verbose = 0, ...) 
{
  if (!inherits(data, "msdata")) 
    stop("'data' must be a 'msdata' object")
  trans <- attr(data, "trans")
  ids <- unique(data[[id]])
  n <- length(ids)
  th <- theta(data, ...)
  
  res1 <- res2 <- res3 <- matrix(NA, nrow(th), B)
  
  for (b in 1:B) {
    if (verbose > 0) {
      cat("\nBootstrap replication", b, "\n")
      flush.console()
    }
    bootdata <- NULL
    bids <- sample(ids, replace = TRUE)
    bidxs <- unlist(sapply(bids, function(x) which(x == 
                                                     data[[id]])))
    bootdata <- data[bidxs, ]
    # if (verbose > 0) {
    #   print(date())
    #   print(events(bootdata))
    #   cat("applying theta ...")
    # }
    thstar <- try(theta(bootdata, ...))
    
    if(class(thstar) == 'try-error'){
      res1[, b] <- NA
      res2[, b] <- NA
      res3[, b] <- NA
    } else {
      res1[, b] <- thstar$pstate1
      res2[, b] <- thstar$pstate2
      res3[, b] <- thstar$pstate3
    }
    
    
    
  }
  if (verbose) 
    cat("\n")
  
  
  return(list(res1 = res1,
              res2 = res2,
              res3 = res3))
}