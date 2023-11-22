
# prevalence analysis for LVA data
# imperfect test and gold-standard
# with BUGS

library(rjags)
library(R2jags)
library(dplyr)

dat_raw <- haven::read_dta("data/LVA and outcomes.dta")

# marginal data
dat <-
  dat_raw |> 
  select(study, imaging, cohort, aneurysm) |> 
  rename(N = cohort,
         n = aneurysm) |> 
  mutate(
    # study = as.numeric(as.factor(study)),
    imaging = ifelse(imaging == "CMR+TTE",  # CMR as gold-standard
                     "CMR", imaging))

# cross test, joint data
# z: CMR
# x: TTE
maron2008 <- 
  tibble::tribble(
    ~study, ~n,	~z,	~x,
    "Maron 2008", 1271, 0, NA,	
    "Maron 2008", 12,	1,	0,
    "Maron 2008", 10,	1,	1,
    "Maron 2008", 6,	NA,	1)

# all data in joint format
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
# from count data
for (j in 1:nrow(dat_agg)) {
  for (i in offset[j]:(offset[j+1]-1)) {
    x[i] <- dat_agg$x[j]
    z[i] <- dat_agg$z[j]
    study_id[i] <- dat_agg$study_id[j]
  }  
}

# priors
prev <- 0.2
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

# fit model
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

# two-step approach
jags_mod <- jags.model(
  file = filein,
  data = dataJags,
  inits = NULL,
  n.chains = 1)

res_jags <- 
  rjags::coda.samples(jags_mod,
                      variable.names = params,
                      n.iter = n.iter,
                      thin = n.thin)

R2WinBUGS::attach.bugs(res_bugs$BUGSoutput)

output <- res_bugs$BUGSoutput
print(output, digits.summary = 4)

save(res_jags, file = "data/res_jags_main.RData")


#########
# plots #
#########

library(tidybayes)
library(ggplot2)

study_lup <- 
  dat_agg |> 
  select(study, study_id) |> 
  distinct()

dat_psi <- spread_draws(res_jags, psi[i]) |> 
  left_join(study_lup, by = join_by(i == study_id)) |> 
  ungroup() |> 
  select(-i)

dat_psi0 <- spread_draws(res_jags, psi0) |> 
  mutate(study = "Average") |> 
  rename(psi = psi0)

plot_dat <- rbind(dat_psi, dat_psi0)

dat_sum <- plot_dat |> 
  group_by(study) %>% 
  mean_qi(psi)

plot_dat %>%   
  ggplot(aes(x = psi, y = study)) +
  geom_vline(xintercept = mean(res_jags[[1]][,"psi0"]), linewidth = 0.25, lty = 2) +
  ggdist::stat_halfeye(.width = c(0.8, 0.95), fill = "dodgerblue") +
  xlab("Prevalence") +
  xlim(0, 0.2) +
  scale_x_continuous(labels = scales::percent)
  # # Add text labels
  # geom_text(
  # data = mutate_if(out_all_sum, is.numeric, round, 2),
  # aes(label = str_glue("{b_Intercept} [{.lower}, {.upper}]"), x = 0.2),
  # hjust = "inward") 
  
write.csv(plot_dat, file = "data/plot_dat_main.csv")
