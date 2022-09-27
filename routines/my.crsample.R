my.crsample = function (Haz, tcond = 0, censtime = NULL) 
{
  if (is.null(censtime)) 
    fut <- Inf
  else fut <- censtime$time
  transs <- Haz$trans
  transun <- unique(transs)
  K <- length(transun)
  tt <- sort(unique(Haz$time))
  n <- length(tt)
  cim <- matrix(NA, n, 3 * K + 4)
  ci <- as.data.frame(cim)
  names(ci)[1] <- "time"
  names(ci)[2:(K + 1)] <- paste("Haz", as.character(1:K), 
                                sep = "")
  names(ci)[(K + 2):(2 * K + 1)] <- paste("haz", as.character(1:K), 
                                          sep = "")
  names(ci)[(2 * K + 2):(3 * K + 1)] <- paste("CI", as.character(1:K), 
                                              sep = "")
  names(ci)[3 * K + 2] <- "hazsum"
  names(ci)[3 * K + 3] <- "Hazsum"
  names(ci)[3 * K + 4] <- "S0"
  ci$time <- tt
  for (k in 1:K) {
    wh <- which(Haz$trans == transun[k])
    idx <- match(Haz$time[wh], tt)
    ci[, k + 1][idx] <- Haz$Haz[wh]
    ci[, k + 1] <- my.NAfix(ci[, k + 1], subst = 0)
    ci[, K + 1 + k] <- diff(c(0, ci[, k + 1]))
  }
  ci <- ci[ci$time >= tcond, ]
  n <- nrow(ci)
  for (k in 1:K) ci[, k + 1] <- cumsum(ci[, K + 1 + k])
  if (K == 1) 
    ci$hazsum <- ci[, 3]
  else ci$hazsum <- apply(ci[, ((K + 2):(2 * K + 1))], 1, 
                          sum)
  ci$S0 <- cumprod(1 - ci$hazsum)
  ci$Hazsum <- -log(ci$S0)
  nci <- nrow(ci)
  k <- NA
  tsample <- my.Hazsample(data.frame(time = ci$time, Haz = ci$Hazsum))
  if (fut < tsample) 
    crt <- fut
  else {
    crt <- tsample
    if (fut > tsample) {
      k <- sample(1:K, size = 1, prob = ci[which(ci$time == 
                                                   tsample), (K + 2):(2 * K + 1)])
    }
    else if (crt != Inf) {
      k <- sample(c(1:K, NA), size = 1, prob = c(ci[which(ci$time == 
                                                            tsample), (K + 2):(2 * K + 1)], censtime$jump))
    }
  }
  if (!is.na(k)) 
    trans <- unique(Haz$trans)[k]
  else trans <- NA
  return(list(t = crt, trans = trans))
}