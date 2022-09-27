my.mssample = function (Haz, trans, history = list(state = 1, time = 0, tstate = NULL), 
                        beta.state = NULL, clock = c("forward", "reset"), output = c("state", 
                                                                                     "path", "data"), tvec, cens = NULL, M = 10, do.trace = NULL) 
{
  output <- match.arg(output)
  clock <- match.arg(clock)
  K <- dim(trans)[1]
  trans2 <- to.trans2(trans)
  ntrans <- nrow(trans2)
  if (length(history$state) == 1) 
    history$state <- rep(history$state, M)
  if (length(history$time) == 1) 
    history$time <- rep(history$time, M)
  if (length(history$state) != length(history$time)) 
    stop("lengths of history$state and history$time differ")
  if (!is.null(history$tstate)) {
    if (is.vector(history$tstate)) 
      if (length(history$tstate) != K) 
        stop("length of history$tstate should equal no of states")
    else history$tstate <- matrix(history$tstate, K, 
                                  M)
    if (is.null(beta.state)) 
      stop("beta.state should be specified when history$tstate not null")
  }
  if (!is.null(beta.state)) 
    if (any(dim(beta.state) != c(K, ntrans))) 
      stop("incorrect dimension of beta.state")
  if (output == "state") 
    res <- matrix(0, length(tvec), K)
  else if (output == "path") {
    thepaths <- paths(trans)
    L <- nrow(thepaths)
    res <- matrix(0, length(tvec), L)
  }
  else res <- NULL
  for (m in 1:M) {
    if(m %% 2500 == 0) print(paste('m =', m, 'paths out of M =', M, 'simulated.', sep = ' '))
    # if(m == 162) debug(my.mssample1)
    if (!is.null(history$tstate)) 
      res1 <- my.mssample1(Haz, trans, history = list(state = history$state[m], 
                                                   time = history$time[m], tstate = history$tstate[, 
                                                                                                   m]), beta.state = beta.state, clock = clock, 
                        output = output, tvec = tvec, cens = cens)
    else res1 <- my.mssample1(Haz, trans, history = list(state = history$state[m], 
                                                      time = history$time[m], tstate = rep(0, K)), beta.state = beta.state, 
                           clock = clock, output = output, tvec = tvec, cens = cens)
    if (output == "data") {
      res1[, 1] <- m
      res <- rbind(res, res1)
    }
    else res <- res + res1
    if (!is.null(do.trace)) 
      if (m%%do.trace == 0) {
        cat("Replication", m, "finished at", date(), 
            "\n")
        flush.console()
      }
  }
  if (output == "state") {
    res <- data.frame(cbind(tvec, res/M))
    names(res) <- c("time", paste("pstate", 1:K, sep = ""))
  }
  else if (output == "path") {
    res <- data.frame(cbind(tvec, res/M))
    names(res) <- c("time", paste("ppath", 1:L, sep = ""))
  }
  else if (output == "data") {
    res <- data.frame(res)
    names(res) <- c("id", "Tstart", "Tstop", "duration", 
                    "from", "to", "status", "trans")
    attr(res, "trans") <- trans
    class(res) <- "msdata"
  }
  return(res)
}