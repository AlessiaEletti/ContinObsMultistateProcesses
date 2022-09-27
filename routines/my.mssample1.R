my.mssample1 = function (Haz, trans, history, beta.state, clock, output, tvec, 
                         cens) 
{
  if (!is.null(cens)) {
    pcens <- diff(c(0, 1 - cens$surv))
    idx <- sample(1:length(cens$time), size = 1, prob = pcens)
    fut <- cens$time[idx]
    censtime <- list(time = fut, jump = ifelse(idx > 1, 
                                               cens$Haz[idx] - cens$Haz[idx - 1], cens$Haz[idx]))
  }
  else censtime <- NULL
  K <- dim(trans)[1]
  trans2 <- to.trans2(trans)
  from <- to <- history$state
  tcond <- t0 <- Tstart <- history$time
  if (output == "state") 
    res <- matrix(0, length(tvec), K)
  else if (output == "path") {
    thepaths <- paths(trans)
    path <- c(to, rep(NA, ncol(thepaths) - 1))
    res <- matrix(0, length(tvec), nrow(thepaths))
  }
  else res <- NULL
  tstates <- history$tstate
  while (!is.na(to)) {
    from <- to
    nstates <- trans[from, ]
    transs <- nstates[!is.na(nstates)]
    allto <- which(!is.na(nstates))
    ntr <- length(transs)
    if (ntr != 0) {
      transnos <- transs
      for (tr in 1:ntr) Haz$Haz[Haz$trans == transnos[tr]] <- exp(sum(beta.state[, 
                                                                                 transnos[tr]] * tstates)) * Haz$Haz[Haz$trans == 
                                                                                                                       transnos[tr]]
      whh <- which(!is.na(match(Haz$trans, transnos)))
      if (clock == "forward") {
        crs <- my.crsample(Haz[whh, ], tcond, censtime)
        tcond <- Tstop <- crs$t
      }
      else {
        crs <- my.crsample(Haz[whh, ], t0, censtime)
        t0 <- 0
        tcond <- Tstop <- crs$t + tcond
      }
      transno <- crs$trans
      if (is.na(transno)) 
        to <- NA
      else {
        to <- trans2$to[transno]
        tstates[to] <- Tstop
      }
      if (output == "state") {
        res[((tvec >= Tstart) & (tvec < Tstop)), from] <- 1
        Tstart <- Tstop
      }
      else if (output == "path") {
        idx <- which(apply(thepaths, 1, function(x) identical(x, 
                                                              path)))
        res[((tvec >= Tstart) & (tvec < Tstop)), idx] <- 1
        path[which(is.na(path))[1]] <- to
        Tstart <- Tstop
      }
      else {
        res1 <- matrix(c(rep(NA, ntr), rep(Tstart, ntr), 
                         rep(Tstop, ntr), rep(Tstop - Tstart, ntr), 
                         rep(from, ntr), allto, rep(0, 2 * ntr)), ntr, 
                       8)
        res1[res1[, 6] == to, 7] <- 1
        res1[, 8] <- trans[from, allto]
        Tstart <- Tstop
        res <- rbind(res, res1)
      }
    }
    else {
      to <- NA
      if (output == "state") {
        res[tvec >= Tstart, from] <- 1
      }
      else if (output == "path") {
        idx <- which(apply(thepaths, 1, function(x) identical(x, 
                                                              path)))
        res[tvec >= Tstart, idx] <- 1
        path[which(is.na(path))[1]] <- to
      }
      else {
        res1 <- matrix(c(rep(NA, ntr), rep(Tstart, ntr), 
                         rep(Tstop, ntr), rep(Tstop - Tstart, ntr), 
                         rep(from, ntr), allto, rep(0, 2 * ntr)), ntr, 
                       8)
        res1[res1[, 6] == to, 7] <- 1
        res1[, 8] <- trans[from, allto]
        res <- rbind(res, res1)
      }
    }
  }
  return(res)
}