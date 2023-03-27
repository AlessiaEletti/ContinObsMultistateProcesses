
# ************************* #
# DATA LOADING AND SETUP ####
# ************************* #
library(rstpm2)
library(readstata13)
library(ggplot2)
mex.1 <- read.dta13("http://fmwww.bc.edu/repec/bocode/m/multistate_example.dta")
transmat <- rbind("Post-surgery"=c(NA,1,2),
                  "Relapsed"=c(NA,NA,3),
                  "Died"=c(NA,NA,NA))
colnames(transmat) <- rownames(transmat)

# --- Dataset pre-processing ---
# transform observed state alive/deceased to 0/1
mex.2 <- transform(mex.1, osi = (osi=="deceased") + 0) 

# fix typo
levels(mex.2$size)[2] <- ">20-50 mm" 

# transforms wide format dataset to long format dataset
mex <- mstate::msprep(time = c(NA,"rf","os"),status=c(NA,"rfi","osi"),
                      data = mex.2,trans=transmat,id="pid",
                      keep = c("age","size","nodes","pr_1","hormon"))
# Comment: status indicates whether the transition represented by the row happened
# When individuals to not transition they will always have two rows, with stautus = 0
# for both transition 1->2 and 1->3 (note that transition 2->3 ONLY appears if 1->2 happened).

# creates dummy variables for size and hormon and coverts time Tstart and Tstop from months to years
# Note: 'time' is still in months and is given by time = Tstop - Tstart.
mex <- transform(mex,
                 size2 = (unclass(size) == 2)+0, # avoids issues with TRUE/FALSE
                 size3 = (unclass(size) == 3)+0,
                 hormon = (hormon == "yes")+0,
                 Tstart = Tstart/12,
                 Tstop = Tstop/12)


# **************** #
# MODEL FITTING ####
# **************** #
library(GJRM)

out.ar = gamlss(list(Tstop ~ s(log(Tstop), bs = 'mpi') + s(age) + size2 + size3 + s(nodes) + s(pr_1) + hormon +
                       ti(log(Tstop), pr_1)), surv = TRUE, margin = 'PH',
                data = mex[mex$trans == 1, ],
                type = 'R',
                cens = status)
conv.check(out.ar)
summary(out.ar)
AIC(out.ar)
BIC(out.ar)

out.ad = gamlss( list(Tstop ~ s(log(Tstop), bs = 'mpi') + s(age) + size2 + size3 + s(nodes) + s(pr_1) + hormon +
                        ti(log(Tstop), pr_1)),
                 surv = TRUE, margin = 'PH',
                 data = mex[mex$trans == 2, ],
                 type = 'R',
                 cens = status)
conv.check(out.ad)
summary(out.ad)
AIC(out.ad)
BIC(out.ad)


# define status indicator
status.factor = ifelse(mex$status == 0, 'R', 'U')
status.factor[mex$Tstart > 0] = paste(status.factor[mex$Tstart > 0], 'T', sep = '')
status.factor = as.factor(status.factor)

out.rd = gamlss( list(Tstop ~ s(log(Tstop), bs = 'mpi') + s(age) + size2 + size3 + s(nodes) + s(pr_1) + hormon +
                        ti(log(Tstop), pr_1)),
                 surv = TRUE, margin = 'PH',
                 data = mex[mex$trans == 3, ],
                 truncation.time = 'Tstart',
                 type = 'mixed',
                 cens = status.factor[mex$trans == 3])
conv.check(out.rd)
summary(out.rd)
AIC(out.rd)
BIC(out.rd)



# Smooths from 1-2 transition
plot(out.ar, eq = 1, pages = 1, pers = TRUE, phi = 10, theta = -70)

# Smooths from 1-3 transition
plot(out.ad, eq = 1, pages = 1, pers = TRUE, phi = 20, theta = -70)

# Smooths from 2-3 transition
plot(out.rd, eq = 1, pages = 1, pers = TRUE, phi = 20, theta = -70)



# ************************************************ #
# PREDICTING ESTIMATED TRANSITION PROBABILITIES ####
# ************************************************ #

# USING MSTATE TO OBTAIN TRANSITION PROBABILITIES FROM GJRM (SINCE WE CAN ESTIMATE CUMULATIVE HAZARD) 
library(mstate)

# For prediction of transition probabilities

par(mfrow = c(3,3),
    mar = c(5, 5.5, 4, 2) + 0.1)

times = seq(0, 15, length = 100) 

nodes.plots = sort(rep(c(0, 10, 20), 3))
size2.plots = rep(c(0, 1, 0), 3)
size3.plots = rep(c(0, 0, 1), 3)

# Modified mstate functions **********
source('routines/my.mssample.R')
source('routines/my.mssample1.R')
source('routines/my.crsample.R')
source('routines/my.NAfix.R')
source('routines/my.Hazsample.R')
# ************************************


for(i in 1:9){
  
  newdata = data.frame(age = 54, size2 = size2.plots[i], size3 = size3.plots[i], nodes = nodes.plots[i], pr_1 = 3, hormon = 1)
  
  pred.ar.test =  hazsurv.plot(out.ar, eq = 1, t.vec = times, newdata = newdata,
                               type = 'cumhaz', plot.out = F, n.sim = 200) 
  
  pred.ad.test =  hazsurv.plot(out.ad, eq = 1, t.vec = times, newdata = newdata,
                               type = 'cumhaz', plot.out = F, n.sim = 200)
  
  pred.rd.test =  hazsurv.plot(out.rd, eq = 1, t.vec = times, newdata = newdata,
                               type = 'cumhaz', plot.out = F, n.sim = 200)
  
  if(i == 1){ # save for later plot (of 95% confidence intervals)
    pred.ar.test.save = pred.ar.test
    pred.ad.test.save = pred.ad.test
    pred.rd.test.save = pred.rd.test
  }
  
  Hazprep = data.frame(time = rep(times, 3), 
                       Haz = c(pred.ar.test$ch, pred.ad.test$ch, pred.rd.test$ch), 
                       trans = c(rep(1, length(times)), # health -> relapse
                                 rep(2, length(times)), # health -> death
                                 rep(3, length(times)))) # relapse -> death
  
  
  set.seed(432)
  probs.test = my.mssample(Haz = Hazprep, trans = transmat, tvec = times, M = 10000)
  
  main = ''
  ylab = 'Probability'
  if(i == 1){ main = 'Size <= 20'; ylab = 'Nodes = 0\nProbability' }
  if(i == 2) main = '20 < Size < 50'
  if(i == 3) main = 'Size >= 50'
  if(i == 4) ylab = 'Nodes = 10\nProbability'
  if(i == 7) ylab = 'Nodes = 20\nProbability'
  
  plot(probs.test$time, probs.test$pstate1, type = 'l', lwd = 2, col = 'grey50', ylim = c(0,1),
       xlab = 'Follow-up time', main = main, ylab = ylab, cex.lab = 1.5)
  polygon(c(probs.test$time, rev(probs.test$time)), c(probs.test$pstate1+probs.test$pstate2, probs.test$pstate1+probs.test$pstate2+probs.test$pstate3),
          col = 'grey90', 
          border = NA)
  lines(probs.test$time, probs.test$pstate1+probs.test$pstate2, lwd = 2, col = 'grey70')
  polygon(c(probs.test$time, rev(probs.test$time)), c(probs.test$pstate1, rev(probs.test$pstate1+probs.test$pstate2)),
          col = 'grey70',
          border = NA)
  polygon(c(rev(probs.test$time), probs.test$time), c(rep(0, length(probs.test$time)), probs.test$pstate1),
          col = 'grey50',
          border = NA)
  lines(probs.test$time, probs.test$pstate1+probs.test$pstate2+probs.test$pstate3, lwd = 2, col = 'grey90')
  
}



# ***************************************************************** #
# OBTAINING CONFIDENCE INTERVALS FOR THE TRANSITION PROBABILTIES ####
# ***************************************************************** #

# Since hazsurv.plot already provides simulated cumulative intensity
# functions (to be used for computation of corresponding
# confidence interval) we can use these to obtain an equal number of probabilities and then get 
# quantiles of these.


# First fix the simulated cumulative transition-specific hazard functions
pred.ar.test.clean = pred.ar.test.save$s.sim[,!apply(apply(pred.ar.test.save$s.sim, 2, is.na), 2, any)]
pred.ad.test.clean = pred.ad.test.save$s.sim[,!apply(apply(pred.ad.test.save$s.sim, 2, is.na), 2, any)]
pred.rd.test.clean = pred.rd.test.save$s.sim[,!apply(apply(pred.rd.test.save$s.sim, 2, is.na), 2, any)]

acceptable.number.of.sims = 100 # to ensure we are using 100 simulations (this can be increased/decreased)
at.most = min(c(ncol(pred.ar.test.clean), ncol(pred.ad.test.clean), ncol(pred.rd.test.clean)), acceptable.number.of.sims)

pred.ar.test.clean = pred.ar.test.clean[,1:at.most]
pred.ad.test.clean = pred.ad.test.clean[,1:at.most]
pred.rd.test.clean = pred.rd.test.clean[,1:at.most]


set.seed(432)
start = Sys.time()
p.ar.sim = p.ad.sim = p.rd.sim = matrix(NA, nrow = nrow(pred.ar.test$s.sim), ncol = acceptable.number.of.sims)
for( i in 1:ncol(pred.ar.test.clean) ){ # this takes between 30 to 40 minutes to complete
  
  Hazprep = data.frame(time = rep(times, 3), 
                       Haz = c(pred.ar.test.clean[,i],
                               pred.ad.test.clean[,i],
                               pred.rd.test.clean[,i]), 
                       trans = c(rep(1, length(times)), 
                                 rep(2, length(times)),
                                 rep(3, length(times))))
  
  pred.sim = try(my.mssample(Haz = Hazprep, trans = transmat, tvec = times, M = 10000))
  if(class(pred.sim) != 'try-error'){ # just in case the computation of P fails - not the case anymore with modified mstate functions, but safety check
    p.ar.sim[,i] = pred.sim$pstate1
    p.ad.sim[,i] = pred.sim$pstate2
    p.rd.sim[,i] = pred.sim$pstate3
  }
  
  if(i %% 20 == 0) print(paste('P computed for', i, 'simulated cumulative hazards out of', ncol(pred.ar.test.clean), sep = ' '))
}
end = Sys.time()
rm(pred.sim)
end-start

p.ar.sim = p.ar.sim[,!apply(apply(p.ar.sim, 2, is.na), 2, any)] # none should be removed, just safety check
p.ar.CI = matrixStats::rowQuantiles(p.ar.sim, probs = c(0.025, 0.975))

p.ad.sim = p.ad.sim[,!apply(apply(p.ad.sim, 2, is.na), 2, any)]
p.ad.CI = matrixStats::rowQuantiles(p.ad.sim, probs = c(0.025, 0.975))

p.rd.sim = p.rd.sim[,!apply(apply(p.rd.sim, 2, is.na), 2, any)]
p.rd.CI = matrixStats::rowQuantiles(p.rd.sim, probs = c(0.025, 0.975))


# Plots of estimated transition probabilities with 95% confidence intervals
par(mfrow = c(1,3))
plot(probs.test$time, probs.test$pstate1, type = 'l', lwd = 1, ylim = c(0,1),
     ylab = 'Probabilities', xlab = 'Follow-up time', cex.lab = 1.5)
lines(probs.test$time, p.ar.CI[,1], lty = 2, lwd = 1)
lines(probs.test$time, p.ar.CI[,2], lty = 2, lwd = 1)

plot(probs.test$time, probs.test$pstate2, type = 'l', lwd = 1, ylim = c(0,1),
     ylab = 'Probabilities', xlab = 'Follow-up time', cex.lab = 1.5)
lines(probs.test$time, p.ad.CI[,1], lty = 2, lwd = 1)
lines(probs.test$time, p.ad.CI[,2], lty = 2, lwd = 1)

plot(probs.test$time, probs.test$pstate3, type = 'l', lwd = 1, ylim = c(0,1),
     ylab = 'Probabilities', xlab = 'Follow-up time', cex.lab = 1.5)
lines(probs.test$time, p.rd.CI[,1], lty = 2, lwd = 1)
lines(probs.test$time, p.rd.CI[,2], lty = 2, lwd = 1)


# ****************************************************************************************************************



# ************************************************ #
# PLOTS OF THE ESTIMATED TRANSITION INTENSITIES ####
# ************************************************ #

cex.lab = 1.5
cex.axis = 1.2

par(mfrow = c(1,1))
# Transition 1-2
hazsurv.plot(out.ar, eq = 1, t.vec = seq(min(mex$Tstop), max(mex$Tstop),
                                         length.out = 1000), newdata = data.frame(age = 54, size2 = 0, size3 = 1, nodes = 10, pr_1 = 3, hormon = 1), type = 'hazard',
             ylab = 'Health to Relapse trans. intensity', xlab = 'Time since surgery (years)', 
             cex.lab = cex.lab, cex.axis = cex.axis, n.sim = 1000) 
rug(mex$Tstop[mex$trans == 1])
abline(v = min(mex$Tstop[mex$trans == 1]), lty = 3)


# Transition 1-3
hazsurv.plot(out.ad, eq = 1, t.vec = seq(min(mex$Tstop), max(mex$Tstop),
                                         length.out = 1000), newdata = data.frame(age = 54, size2 = 0, size3 = 1, nodes = 10, pr_1 = 3, hormon = 1), type = 'hazard',
             ylab = 'Health to Death trans. intensity', xlab = 'Time since surgery (years)', 
             cex.lab = cex.lab, cex.axis = cex.axis)
rug(mex$Tstop[mex$trans == 2])
abline(v = min(mex$Tstop[mex$trans == 2]), lty = 3)


# Transition 2-3
hazsurv.plot(out.rd, eq = 1, t.vec = seq(min(mex$Tstop), max(mex$Tstop), length.out = 1000), newdata = data.frame(age = 54, size2 = 0, size3 = 1, nodes = 10, pr_1 = 3, hormon = 1), type = 'hazard',
             ylab = 'Relapse to Death trans. intensity', xlab = 'Time since surgery (years)', 
             cex.lab = cex.lab, cex.axis = cex.axis)
rug(mex$Tstop[mex$trans == 3])
abline(v = min(mex$Tstop[mex$trans == 3]), lty = 3)


# ****************************** #
# NO COVARIATES MODEL FITTING ####
# ****************************** #
library(GJRM)

# This is for Figure 6 in the Supplementary Material
c.ar.GJRM.noCov = gamlss( list( Tstop ~ s(log(Tstop), bs = 'mpi') ),
                          surv = TRUE, margin = 'PH',
                          data = mex[mex$trans == 1, ],
                          # truncation.time = 'Tstart',
                          type = 'R',
                          cens = status) # no problem when using type = 'mixed' and cens = status.factor (anyway no truncation occurs here)

c.ad.GJRM.noCov = gamlss( list(Tstop ~ s(log(Tstop), bs = 'mpi') ),
                          surv = TRUE, margin = 'PH',
                          data = mex[mex$trans == 2, ],
                          # truncation.time = 'Tstart',
                          type = 'R',
                          cens = status) # no problem when using type = 'mixed' and cens = status.factor (anyway no truncation occurs here)

# define status indicator
status.factor = ifelse(mex$status == 0, 'R', 'U')
status.factor[mex$Tstart > 0] = paste(status.factor[mex$Tstart > 0], 'T', sep = '')
status.factor = as.factor(status.factor)

c.rd.GJRM.noCov = gamlss( list(Tstop ~ s(log(Tstop), bs = 'mpi') ),
                          surv = TRUE, margin = 'PH',
                          data = mex[mex$trans == 3, ],
                          truncation.time = 'Tstart',
                          type = 'mixed',
                          cens = status.factor[mex$trans == 3])

# For plotting of estimated cumulative hazard functions (when no covariate were included in model fitting)
times = seq(0, max(mex$Tstop), length = 2000) 
par(mfrow = c(2,2))
pred.ar.test =  hazsurv.plot(c.ar.GJRM.noCov, eq = 1, t.vec = times, newdata = data.frame(age = 54, size2 = 0, size3 = 0, nodes = 10, pr_1 = 3,
                                                                                          hormon = 1), 
                             baseline = F,
                             type = 'cumhaz', pch = 19,
                             plot.out = T, xlab = 'Follow-up time (years since surgery)', ylab = 'Baseline Cumulative Hazard',
                             ylim = c(0, 4.5), intervals = F)
grid(nx = NA, ny = NULL)

pred.ad.test =  hazsurv.plot(c.ad.GJRM.noCov, eq = 1, t.vec = times, newdata = data.frame(age = 54, size2 = 0, size3 = 0, nodes = 10, pr_1 = 3,
                                                                                          hormon = 1), 
                             baseline = TRUE,
                             type = 'cumhaz', pch = 19,
                             plot.out = T, xlab = 'Follow-up time (years since surgery)', ylab = 'Baseline Cumulative Hazard',
                             ylim = c(0, 4.5), intervals = F)
grid(nx = NA, ny = NULL)

pred.rd.test =  hazsurv.plot(c.rd.GJRM.noCov, eq = 1, t.vec = times, newdata = data.frame(age = 54, size2 = 0, size3 = 0, nodes = 10, pr_1 = 3,
                                                                                          hormon = 1), 
                             baseline = TRUE,
                             type = 'cumhaz', pch = 19,
                             plot.out = T, xlab = 'Follow-up time (years since surgery)', ylab = 'Baseline Cumulative Hazard',
                             ylim = c(0,4.5), intervals = F)
grid(nx = NA, ny = NULL)











