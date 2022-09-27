my.NAfix = function (x, subst = -Inf) 
{
  spec <- max(x[!is.na(x)]) + 1
  x <- c(spec, x)
  while (any(is.na(x))) x[is.na(x)] <- x[(1:length(x))[is.na(x)] - 
                                           1]
  x[x == spec] <- subst
  x <- x[-1]
  x
}