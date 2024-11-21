
# frequentist and Bayesian meta-analysis
# of aneurysm (LVA) data
# for outcomes:
#  stroke, LV thrombus, SCD, imaging, size


library(meta)
library(dplyr)
library(lme4)

filename <- "LVA and outcomes_2023 11 30.dta"
# filename <- "LVA and outcomes.dta"

dat_raw <- foreign::read.dta(glue::glue("data/{filename}"))
# write.csv(dat_raw, file = "data/LVA-and-outcomes.csv")

dat_raw$nsvt_aneu_n <- as.numeric(dat_raw$nsvt_aneu_n)

##############
# frequentist
##############

## prevalence in total cohort

res_aneurysm <-
  metaprop(event = aneurysm, n = cohort, studlab = study, data = dat_raw)

res_scd <-
  metaprop(event = nscd, n = cohort, studlab = study, data = dat_raw)

res_imaging <-
  metaprop(event = aneurysm, n = cohort, studlab = study, byvar = imaging, data = dat_raw)

res_small <-
  metaprop(event = n_small, n = cohort, studlab = study, data = dat_raw)

res_medium <-
  metaprop(event = n_medium, n = cohort, studlab = study, data = dat_raw)

res_large <-
  metaprop(event = n_large, n = cohort, studlab = study, data = dat_raw)

dat_size <- dat_raw |> 
  reshape2:::melt.data.frame(measure.vars = c("n_small", "n_medium", "n_large"),
                             variable.name = "size")
res_size <-
  metaprop(event = value, n = cohort, studlab = study, byvar = size, data = dat_size)

resbind_size <-
  metabind(res_small, res_medium, res_large,
           outclab = "", pooled = "common", backtransf = FALSE)
##TODO: error
# forest(resbind_size, print.I2 = FALSE, print.pval.Q = FALSE, print.subgroup.labels = FALSE)

## prevalence in aneurysm group

res_stroke <-
  metaprop(event = ncva, n = aneurysm, studlab = study, data = dat_raw)

res_lvthrombus <-
  metaprop(event = nlvthrombus, n = aneurysm, studlab = study, data = dat_raw)

res_svt_aneu <-
  metaprop(event = nsvt_aneu_n, n = aneurysm, studlab = study, data = dat_raw)

res_scd_per_aneurysm <-
  metaprop(event = nscd, n = aneurysm, studlab = study, byvar = atpy_n, data = dat_raw)

res_small_per_aneurysm <-
  metaprop(event = n_small, n = aneurysm, studlab = study, data = dat_raw)

res_medium_per_aneurysm <-
  metaprop(event = n_medium, n = aneurysm, studlab = study, data = dat_raw)

res_large_per_aneurysm <-
  metaprop(event = n_large, n = aneurysm, studlab = study, data = dat_raw)


# The same comment applies for sudden death; can the forest plot show n_scd/aneurysm 
# (rather than n_scd/cohort? Could this be stratified with variable atpy_n (which I have included to the attached dataset).
#   ATP is a surrogate for SCD which erroneously inflates SCD outcomes.


# SCD in patients with small and big LVA

dat_raw$nscd_big <- dat_raw$nscd_medium + dat_raw$nscd_large
dat_raw$n_big <- dat_raw$n_medium + dat_raw$n_large
dat_raw$n_big <- ifelse(is.na(dat_raw$n_big), 100, dat_raw$n_big)

dat_size <- dat_raw |> 
  reshape2:::melt.data.frame(measure.vars = c("n_small", "n_big"),
                             variable.name = "size_label", id.vars = "study")
dat_scdsize <- dat_raw |> 
  reshape2:::melt.data.frame(measure.vars = c("nscd_small", "nscd_big"),
                             variable.name = "size_label", id.vars = "study") |> 
  mutate(size_label = gsub(replacement = "", "scd", size_label))

dat_scdsize <- merge(dat_scdsize, dat_size, by = c("study", "size_label")) |> 
  rename(nscd = value.x, n = value.y) |> 
  filter(!is.na(nscd))

res_nscd_big <-
  metaprop(event = nscd_big, n = n_big, studlab = study, data = dat_raw)

res_nscd_size <-
  metaprop(event = nscd, n = n, studlab = study, byvar = size_label, data = dat_scdsize)

#########
# plots #
#########

# custom plot
forest_plot <- function(x, save = FALSE, ...) {
  
  if (save) {
    var_name <- deparse(substitute(x))
    png(glue::glue("plots/{var_name}.png"), height = 500, width = 550)
    on.exit(dev.off())
  }  
  
  weight <- x$n / max(x$n)              # linear
  # weight <- exp(1 + x$n / max(x$n))    # exponential 
  # weight <- log(1 + x$n)               # logarithmic
  
  x$w.random <- weight
  meta::forest(x, weight.study = "random", ...) #, xlim = c(0, 0.1))
}


forest_plot(res_aneurysm)
forest_plot(res_stroke)
forest_plot(res_lvthrombus)
forest_plot(res_svt_aneu)
forest_plot(res_scd)
forest_plot(res_imaging)
forest_plot(res_small)
forest_plot(res_medium)
forest_plot(res_large)

# don't think that this plot is strictly correct because overall pooling is double counting
# so should remove this
forest_plot(res_size)

forest_plot(res_small_per_aneurysm)
forest_plot(res_medium_per_aneurysm)
forest_plot(res_large_per_aneurysm)
forest_plot(res_scd_per_aneurysm)
# forest(res_nscd_big)
forest_plot(res_nscd_size)


#####################
# fit model directly

modf <- lme4::glmer(cbind(aneurysm , cohort-aneurysm) ~ 1 + (1|study),
                    family = "binomial",
                    data = dat_raw)
summary(modf)
ranef(modf)
fixef(modf)
exp(fixef(modf))


###########
# Bayesian
###########

library(rstanarm)
library(brms)

modb <- stan_glmer(cbind(aneurysm , cohort-aneurysm) ~ 1 + (1|study),
                   data = dat_raw,
                   family = binomial(link = "logit"))

# average
exp(modb$coefficients[1])

# same model in brms
mod_brm <- brm(aneurysm | trials(cohort) ~ 1 + (1|study),
               data = dat_raw,
               family = binomial())

save(mod_brm, file = "data/mod_brm_pooled.RData")

#########
# plots #
#########

# https://mvuorre.github.io/posts/2016-09-29-bayesian-meta-analysis/
library(ggplot2)
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

ggsave(filename = "plots/forest_plot_posterior_aneurysm.png",
       width = 20, height = 20, units = "cm", dpi = 640, device = "png")

write.csv(out_all, file = "data/plot_dat_pooled.csv")
