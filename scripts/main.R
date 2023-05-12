
# prevalence analysis
# imperfect test and gold-standard
# with BUGS

library(rjags)
library(R2jags)

dat <- haven::read_dta("data/LVA and outcomes.dta")

# fake data

fdat <- data.frame(z = c(0,0,1,1,NA,NA),
                   x = c(0,1,0,1,0,1),
                   n = c(12,3,5,18,10,20))

offset <- c(0, cumsum(fdat$n)) + 1

x <- z <- NULL

for (j in 1:nrow(fdat)) {
  for (i in offset[j]:(offset[j+1]-1)) {
    x[i] <- fdat$x[j]
    z[i] <- fdat$z[j]
  }  
}

filein <- "BUGS/bugs_script.txt"
params <- c("phi", "psi")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)

dataJags <-
  list(x = x,
       z = z,
       n = length(x))

res_bugs <-
  jags(data = dataJags,
       inits = NULL,
       parameters.to.save = params,
       model.file = filein,
       n.chains = 1,
       n.iter,
       n.burnin,
       n.thin,
       DIC = TRUE)

# print(res_bugs)
# mcmcplot(res_bugs)
# plots <- traceplot(res_bugs)

R2WinBUGS::attach.bugs(res_bugs$BUGSoutput)

output <- res_bugs$BUGSoutput
