my.Hazsample = function (Haz, size = 1, replace = TRUE) 
{
  p <- diff(c(0, 1 - exp(-Haz$Haz))) # this is a CDF
  p <- c(p, exp(-Haz$Haz[nrow(Haz)]))
  # print(length(p))
  return(sample(c(Haz$time, Inf), size = size, prob = p, replace = replace))
}