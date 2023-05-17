
# meta-analysis
# of LVT data


library(meta)
library(dplyr)

dat_raw <- foreign::read.dta("data/LVA and outcomes.dta")

# frequentist

res <- metaprop(event = aneurysm, n = cohort, studlab = study, data = dat_raw)
res_imaging <- metaprop(event = aneurysm, n = cohort, studlab = study, byvar = imaging, data = dat_raw)
res_small <- metaprop(event = n_small, n = cohort, studlab = study, data = dat_raw)

forest(res, xlim = c(0, 0.1))
forest(res_imaging, xlim = c(0, 0.1))
forest(res_small, xlim = c(0, 0.1))

require(lme4)

modf <- glmer(cbind(aneurysm , cohort-aneurysm) ~ 1 + (1|study), family="binomial", data = dat_raw)
summary(modf)
ranef(modf)
fixef(modf)
exp(fixef(modf))

# Bayesian

library(rstanarm)
library(ggplot2)
library(brms)

modb <- stan_glmer(cbind(aneurysm , cohort-aneurysm) ~ 1 + (1|study),
                   data = dat_raw,
                   family = binomial(link = "logit"))

# average
exp(modb$coefficients[1])

mod_brm <- brm(aneurysm | trials(cohort) ~ 1 + (1|study),
               data = dat_raw,
               family = binomial())

save(mod_brm, file = "data/mod_brm_pooled.RData")

s#########
# plots #
#########

# https://mvuorre.github.io/posts/2016-09-29-bayesian-meta-analysis/
library(tidybayes)
library(ggdist)
library(forcats)
library(stringr)

# Study-specific effects are deviations + average
out_r <-
  spread_draws(mod_brm, r_study[study,term], b_Intercept) %>% 
  mutate(b_Intercept = exp(r_study + b_Intercept)) 

# Average effect
out_f <- spread_draws(mod_brm, b_Intercept) %>% 
  mutate(b_Intercept = exp(b_Intercept),
         study = "Average")

# Combine average and study-specific effects' data frames
out_all <- bind_rows(out_r, out_f) %>% 
  ungroup() %>%
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(study = fct_relevel(study, "Average")) %>% 
  # tidybayes garbles names so fix here
  mutate(study = str_replace_all(study, "\\.", " "))

# Data frame of summary numbers
out_all_sum <- group_by(out_all, study) %>% 
  mean_qi(b_Intercept)

out_all %>%   
  ggplot(aes(b_Intercept, study)) +
  geom_vline(xintercept = mean(out_f$b_Intercept), linewidth = 0.25, lty = 2) +
  stat_halfeye(.width = c(0.8, 0.95), fill = "dodgerblue") +
  scale_x_continuous(labels = scales::percent) +
  # Add text labels
  geom_text(
    data = mutate_if(out_all_sum, is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept} [{.lower}, {.upper}]"), x = 0.2),
    hjust = "inward") +
  xlab("Prevalence") #+
  # # Observed as empty points
  # geom_point(
  #   data = dat %>% mutate(study = str_replace_all(study, "\\.", " ")), 
  #   aes(x=yi), position = position_nudge(y = -.2), shape = 1)

write.csv(out_all, file = "data/plot_dat_pooled.csv")
