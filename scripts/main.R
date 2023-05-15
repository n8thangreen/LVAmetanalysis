
# prevalence analysis
# imperfect test and gold-standard
# with BUGS

library(rjags)
library(R2jags)
library(dplyr)

dat <- haven::read_dta("data/LVA and outcomes.dta")

# fake data

# hard-coded
fdat <- data.frame(z = c(0,0,1,1,NA,NA,0,1),
                   x = c(0,1,0,1,0,1,NA,NA),
                   n = c(12,3,5,18,10,20,5,10),
                   study = c(1,1,1,1,2,2,3,3))

offset <- c(0, cumsum(fdat$n)) + 1

x <- z <- NULL
study <- NULL

# create equivalent individual data
for (j in 1:nrow(fdat)) {
  for (i in offset[j]:(offset[j+1]-1)) {
    x[i] <- fdat$x[j]
    z[i] <- fdat$z[j]
    study[i] <- fdat$study[j]
  }  
}

dataJags <-
  list(x = x,                       # imperfect test
       z = z,                       # gold standard
       study = study,               # study id
       ns = length(unique(study)),  # number of studies
       n = length(x))               # number of patients

# from generative distribution
prev <- 0.7
indiv_dat <- 
  data.frame(study = rep(1:3, times = c(10,10,10))) |> 
  group_by(study) |> 
  mutate(id = 1:n(),
         sprev = boot::inv.logit(rnorm(1, car::logit(prev), 0.5))) |> 
  ungroup() |> 
  mutate(z = rbinom(n = n(), size = 1, prob = sprev),
         sens = 0.9,
         spec = 0.8) |> 
  rowwise() |> 
  mutate(x = ifelse(z == 1,
                    rbinom(n = 1, size = 1, prob = sens),
                    rbinom(n = 1, size = 1, prob = 1-spec)),
         x = ifelse(study == 2, NA, x),
         z = ifelse(study == 3, NA, z))

# priors
beta_pars_sens <- MoM_beta(0.9, 0.05)
beta_pars_1m_spec <- MoM_beta(1 - 0.8, 0.1)
mu_beta <- car::logit(prev)

dataJags <-
  list(x = indiv_dat$x,                       
       z = indiv_dat$z,                       
       study = indiv_dat$study,               
       ns = length(unique(indiv_dat$study)),  
       n = nrow(indiv_dat))

filein <- "BUGS/bugs_script.txt"
params <- c("phi", "psi", "psi0", "beta0")

n.iter <- 20000
n.burnin <- 1000
n.thin <- floor((n.iter - n.burnin)/500)

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
output
