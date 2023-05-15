
# prevalence analysis
# imperfect test and gold-standard
# with BUGS

library(rjags)
library(R2jags)
library(dplyr)

dat_raw <- haven::read_dta("data/LVA and outcomes.dta")

dat <-
  dat_raw |> 
  select(study, imaging, cohort, aneurysm) |> 
  rename(N = cohort,
         n = aneurysm) |> 
  mutate(
    # study = as.numeric(as.factor(study)),
    imaging = ifelse(imaging == "CMR+TTE",
                     "CMR", imaging))

maron2008 <- 
  tibble::tribble(
    ~study, ~n,	~z,	~x,
    "Maron 2008", 1271, 0, NA,	
    "Maron 2008", 12,	1,	0,
    "Maron 2008", 10,	1,	1,
    "Maron 2008", 6,	NA,	1)

dat_agg <- 
  rbind(
    dat |> 
      mutate(z = ifelse(imaging == "CMR", 1, NA),
             x = ifelse(imaging == "TTE", 1, NA)),
    dat |> 
      mutate(n = N - n,
             z = ifelse(imaging == "CMR", 0, NA),
             x = ifelse(imaging == "TTE", 0, NA))) |> 
  select(-imaging, -N) |> 
  arrange(study) |> 
  filter(study != "Maron 2008") |> 
  rbind(maron2008) |> 
  mutate(study_id = as.numeric(as.factor(study)))

offset <- c(0, cumsum(dat_agg$n)) + 1

x <- z <- NULL
study_id <- NULL

# create equivalent individual data
for (j in 1:nrow(dat_agg)) {
  for (i in offset[j]:(offset[j+1]-1)) {
    x[i] <- dat_agg$x[j]
    z[i] <- dat_agg$z[j]
    study_id[i] <- dat_agg$study_id[j]
  }  
}

prev <- 0.2

# priors
beta_pars_sens <- MoM_beta(0.8, 0.1)
beta_pars_1m_spec <- MoM_beta(1 - 0.8, 0.1)
mu_beta <- car::logit(prev)

dataJags <-
  list(x = x,                       
       z = z,                       
       study = study_id,               
       ns = length(unique(study_id)),  
       n = length(x))

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
