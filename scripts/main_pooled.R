
# frequentist and Bayesian meta-analysis
# of aneurysm (LVA) data
# for outcomes:
#  stroke, LV thrombus, SCD, imaging, size
#
# using Freeman-Tukey double arcsine transformation

library(meta)
library(dplyr)
library(lme4)
library(tidyr)

filename <- "LVA and outcomes_size 2024 01 28.dta"
# filename <- "LVA and outcomes_2023 11 30.dta"
# filename <- "LVA and outcomes.dta"

dat_raw <- foreign::read.dta(here::here(glue::glue("data/{filename}")))
# write.csv(dat_raw, file = "data/LVA-and-outcomes.csv")

dat_raw$nsvt_aneu_n <- as.numeric(dat_raw$nsvt_aneu_n)

# remove duplicate rows
dat_raw <- dat_raw[!duplicated(dat_raw$study), ]

# # replace NAs with 0
# dat_raw <- dat_raw |> 
#   mutate(across(
#     c(aneurysm, ncva, nlvthrombus, nsvt_aneu_n, nscd,
#       n_small, n_medium, n_large, nscd_small, nscd_medium, nscd_large, ncva_small, ncva_big, nthrombus_small, nthrombus_big),
#     ~ replace_na(.x, 0)))


##############
# frequentist
##############

## prevalence in _total_ cohort

trans_method <- "PFT"
# trans_method <- "PLOGIT"

# remove studies with NAs
res_aneurysm <-
  dat_raw[!is.na(dat_raw$aneurysm), ] |> 
  metaprop(event = aneurysm, n = cohort, studlab = study, sm = trans_method, method.tau = "REML",
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_scd <-
  dat_raw[!is.na(dat_raw$nscd), ] |> 
  metaprop(event = nscd, n = cohort, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_imaging <-
  dat_raw[!is.na(dat_raw$aneurysm), ] |> 
  metaprop(event = aneurysm, n = cohort, studlab = study, sm = trans_method, subgroup = imaging, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_small <-
  dat_raw[!is.na(dat_raw$n_small), ] |> 
  metaprop(event = n_small, n = cohort, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_medium <-
  dat_raw[!is.na(dat_raw$n_medium), ] |> 
  metaprop(event = n_medium, n = cohort, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_large <-
  dat_raw[!is.na(dat_raw$n_large), ] |> 
  metaprop(event = n_large, n = cohort, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_size <-
  dat_raw |> 
  reshape2:::melt.data.frame(measure.vars = c("n_small", "n_medium", "n_large"),
                             variable.name = "size") |> 
  metaprop(event = value, n = cohort, studlab = study, sm = trans_method, subgroup = size, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

resbind_size <-
  metabind(res_small, res_medium, res_large,
           outclab = "", pooled = "common", backtransf = FALSE)
##TODO: error
# forest(resbind_size, print.I2 = FALSE, print.pval.Q = FALSE, print.subgroup.labels = FALSE)

## prevalence in _aneurysm_ group

res_stroke <-
  dat_raw[!is.na(dat_raw$ncva), ] |> 
  metaprop(event = ncva, n = aneurysm, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_lvthrombus <-
  dat_raw[!is.na(dat_raw$nlvthrombus), ] |> 
  metaprop(event = nlvthrombus, n = aneurysm, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_svt_aneu <-
  dat_raw[!is.na(dat_raw$nsvt_aneu_n), ] |> 
  metaprop(event = nsvt_aneu_n, n = aneurysm, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_scd_in_lvaa <-
  dat_raw[!is.na(dat_raw$nscd), ] |> 
  metaprop(event = nscd, n = aneurysm, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)           

res_scd_per_aneurysm <-
  dat_raw[!is.na(dat_raw$nscd), ] |> 
  metaprop(event = nscd, n = aneurysm, studlab = study, sm = trans_method, subgroup = atpy_n, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_small_per_aneurysm <-
  dat_raw[!is.na(dat_raw$n_small), ] |> 
  metaprop(event = n_small, n = aneurysm, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_medium_per_aneurysm <-
  dat_raw[!is.na(dat_raw$n_medium), ] |> 
  metaprop(event = n_medium, n = aneurysm, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_large_per_aneurysm <-
  dat_raw[!is.na(dat_raw$n_large), ] |> 
  metaprop(event = n_large, n = aneurysm, studlab = study, sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

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

non_zero_studies <- dat_raw$n_big != 0 & !is.na(dat_raw$n_big)
res_nscd_big <-
  metaprop(event = nscd_big, n = n_big, studlab = study, data = dat_raw[non_zero_studies, ],
           sm = trans_method, method.tau = "REML",
           common = FALSE,
           backtransf = TRUE,  # proportions
           method.ci = "CP"   # exact binomial confidence intervals
  )

non_zero_studies <- dat_scdsize$n != 0 & !is.na(dat_scdsize$n)
res_nscd_size <-
  metaprop(event = nscd, n = n, studlab = study, subgroup = size_label, data = dat_scdsize[non_zero_studies, ],
           sm = trans_method, method.tau = "REML",
           common = FALSE,
           backtransf = TRUE,  # proportions
           method.ci = "CP"   # exact binomial confidence intervals
  )

###################
# pool odds-ratios

# scd
non_zero_studies <-
  !(dat_raw$nscd_small == 0 & dat_raw$nscd_big == 0) &
  !(is.na(dat_raw$nscd_small) & is.na(dat_raw$nscd_big))

res_scd_size_or <- metabin(
  event.e = nscd_big,        # Events in the treatment group
  n.e = n_big,               # Total number in treatment group
  event.c = nscd_small,      # Events in the control group
  n.c = n_small,             # Total number in control group
  studlab = study,           # Study labels
  sm = "OR",                 # Summary measure: odds ratio (OR)
  method = "MH",             # Mantel-Haenszel method for pooling
  data = dat_raw[non_zero_studies, ],
  method.tau = "REML",
  common = FALSE
)
res_scd_size_or$label.e <- "Big"
res_scd_size_or$label.c <- "Small" 

# cva
non_zero_studies <-
  !(dat_raw$ncva_small == 0 & dat_raw$ncva_big == 0) &
  !(is.na(dat_raw$ncva_small) & is.na(dat_raw$ncva_big))

res_cva_size_or <- metabin(
  event.e = ncva_big,
  n.e = n_big,       
  event.c = ncva_small,  
  n.c = n_small,     
  studlab = study,     
  sm = "OR",           
  method = "MH",
  data = dat_raw[non_zero_studies, ],
  method.tau = "REML",
  common = FALSE
)
res_cva_size_or$label.e <- "Big"
res_cva_size_or$label.c <- "Small" 

# thrombi
non_zero_studies <-
  !(dat_raw$nthrombus_small == 0 & dat_raw$nthrombus_big == 0) &
  !(is.na(dat_raw$nthrombus_small) & is.na(dat_raw$nthrombus_big))

res_thrombi_size_or <- metabin(
  event.e = nthrombus_big,
  n.e = n_big,       
  event.c = nthrombus_small,  
  n.c = n_small,        
  studlab = study,     
  sm = "OR",           
  method = "MH",
  data = dat_raw[non_zero_studies, ],
  method.tau = "REML",
  common = FALSE
)
res_thrombi_size_or$label.e <- "Big"
res_thrombi_size_or$label.c <- "Small" 

###################
# pool odds ratios using Peto's method

# scd
non_zero_studies <- 
  !(dat_raw$nscd_small == 0 & dat_raw$nscd_big == 0) & 
  !(is.na(dat_raw$nscd_small) & is.na(dat_raw$nscd_big))

res_scd_size_or_peto <- metabin(
  event.e = nscd_big,
  n.e = n_big,
  event.c = nscd_small,
  n.c = n_small,
  studlab = study,
  sm = "OR",
  method = "Peto",
  data = dat_raw[non_zero_studies, ],
  common = TRUE,
  random = FALSE
)

res_scd_size_or_peto$label.e <- "Big"
res_scd_size_or_peto$label.c <- "Small"

# cva
non_zero_studies <- 
  !(dat_raw$ncva_small == 0 & dat_raw$ncva_big == 0) & 
  !(is.na(dat_raw$ncva_small) & is.na(dat_raw$ncva_big))

res_cva_size_or_peto <- metabin(
  event.e = ncva_big,  
  n.e = n_big,         
  event.c = ncva_small,
  n.c = n_small,       
  studlab = study,     
  sm = "OR",
  method = "Peto",
  data = dat_raw[non_zero_studies, ],
  common = TRUE,
  random = FALSE
)
res_cva_size_or_peto$label.e <- "Big"
res_cva_size_or_peto$label.c <- "Small"

# thrombi
non_zero_studies <- 
  !(dat_raw$nthrombus_small == 0 & dat_raw$nthrombus_big == 0) & 
  !(is.na(dat_raw$nthrombus_small) & is.na(dat_raw$nthrombus_big))

res_thrombi_size_or_peto <- metabin(
  event.e = nthrombus_big,  
  n.e = n_big,         
  event.c = nthrombus_small,
  n.c = n_small,       
  studlab = study,     
  sm = "OR",
  method = "Peto",
  data = dat_raw[non_zero_studies, ],
  common = TRUE,
  random = FALSE
)
res_thrombi_size_or_peto$label.e <- "Big"
res_thrombi_size_or_peto$label.c <- "Small"

#########
# plots #
#########

# custom plot
forest_plot <- function(x,
                        save = FALSE,
                        filetxt = "",
                        colvars = c("effect", "ci", "w.random"),  #, "Var"),
                        rhs_text = "Treatment",
                        lhs_text = "Control", 
                        exactCI = FALSE,
                        pooledCP = FALSE, ...) {
  
  inv_logit <- function(x) {
    1 / (1 + exp(-x))
  }
  
  if (save) {
    var_name <- deparse(substitute(x)) 
    png(glue::glue("plots/{var_name}{filetxt}.png"), height = 500, width = 650)
    on.exit(dev.off())
  }  
  
  x$Var <- x$seTE^2
  
  if (x$sm == "OR") {
    weight <- 1 / x$Var
    
    if (exactCI) {
      # calculate exact CI
      tab <- array(c(x$event.e, x$n.e - x$event.e,
                     x$event.c, x$n.c - x$event.c),
                   dim = c(x$k, 2, 2))
      
      for (k in 1:x$k) {
        studytab <- t(tab[k, , ])
        
        if (any(studytab == 0)) {
          studytab <- studytab + 1  # has to be integer so can't use 0.5
        }
        
        exact <- exact2x2::exact2x2(studytab)
        x$lower[k] <- log(exact$conf.int[1]) 
        x$upper[k] <- log(exact$conf.int[2])
      }
    }
  } else {
    weight <- x$n / max(x$n)  # linear
    
    if (pooledCP) {
      # clopper-pearson for pooled rate
      alpha <- 0.05
      total_successes <- sum(x$event)
      total_trials <- sum(x$n)
      pooled_prop <- total_successes / total_trials
      
      # TE.random <- qbeta(0.5, total_successes, total_trials - total_successes)
      # lower.random <- qbeta(alpha / 2, total_successes, total_trials - total_successes + 1)
      # upper.random <- qbeta(1 - alpha / 2, total_successes + 1, total_trials - total_successes)
      
      # Clopper-Pearson confidence interval for pooled proportion
      ci <- binom::binom.confint(total_successes, total_trials, method = "exact")
      
      TE.random <- ci["mean"]
      lower.random <- ci["lower"]
      upper.random <- ci["upper"]
      
      x$TE.random <- p2asin(TE.random)
      x$lower.random <- p2asin(lower.random)
      x$upper.random <- p2asin(upper.random)
    }
  }
  # weight <- exp(1 + x$n / max(x$n))    # exponential 
  # weight <- log(1 + x$n)               # logarithmic
  
  sd <- backtrans_delta_PFT(x$TE.random, x$tau2)^0.5
    
  x$w.random <- weight
  meta:::forest.meta(
    x,
    weight.study = "random",
    label.left = glue::glue("Favours {lhs_text}"),
    label.right = glue::glue("Favours {rhs_text}"),
    # text.add = paste("Variance:", labels_var),
    # prediction = TRUE,
    rightcols = colvars,
    text.addline1 = paste0("sd = ", sd),
    ...) #, xlim = c(0, 0.1))
}

#
backtrans_delta_PFT <- function(ft_value, var_ft) {
  g <- function(x) sin(x / 2)^2        # Back-transform to proportion
  g_prime <- function(x) sin(x) / 2    # Derivative of g(x)
  
  # Variance on original scale
  (g_prime(ft_value)^2) * var_ft
}


forest_plot(res_aneurysm, filetxt = trans_method)
forest_plot(res_aneurysm, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
# grid::grid.text("Favours Control", x = 0.3, y = 0.1)
# grid::grid.text("Favours Treatment", x = 0.3, y = 0.1)

forest_plot(res_scd_in_lvaa, filetxt = trans_method)
forest_plot(res_stroke, filetxt = trans_method)
forest_plot(res_lvthrombus, filetxt = trans_method)
forest_plot(res_svt_aneu, filetxt = trans_method)
forest_plot(res_scd, filetxt = trans_method)
forest_plot(res_imaging, filetxt = trans_method)
forest_plot(res_small, filetxt = trans_method)
forest_plot(res_medium, filetxt = trans_method)
forest_plot(res_large, filetxt = trans_method)

forest_plot(res_scd_in_lvaa, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_stroke, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_lvthrombus, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_svt_aneu, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_scd, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_imaging, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_small, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_medium, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)
forest_plot(res_large, filetxt = paste0(trans_method, "_pooledCP"), pooledCP = TRUE)

## odds-ratios

forest_plot(res_scd_size_or, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small")
forest_plot(res_cva_size_or, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small")
forest_plot(res_thrombi_size_or, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small")

forest_plot(res_scd_size_or_peto, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small")
forest_plot(res_cva_size_or_peto, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small")
forest_plot(res_thrombi_size_or_peto, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small")

# exact study ci
forest_plot(res_scd_size_or, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small",
            exactCI = TRUE, filetxt = "_exactCI")
forest_plot(res_cva_size_or, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small",
            exactCI = TRUE, filetxt = "_exactCI")
forest_plot(res_thrombi_size_or, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small",
            exactCI = TRUE, filetxt = "_exactCI")

forest_plot(res_scd_size_or_peto, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small",
            exactCI = TRUE, filetxt = "_exactCI")
forest_plot(res_cva_size_or_peto, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small",
            exactCI = TRUE, filetxt = "_exactCI")
forest_plot(res_thrombi_size_or_peto, colvars = c("effect", "ci", "Var"),
            plotwidth = "3cm", lhs_text = "Big", rhs_text = "Small",
            exactCI = TRUE, filetxt = "_exactCI")

# don't think that this plot is strictly correct because overall pooling is double counting
# so should remove this
forest_plot(res_size, filetxt = trans_method)

forest_plot(res_small_per_aneurysm, filetxt = trans_method)
forest_plot(res_medium_per_aneurysm, filetxt = trans_method)
forest_plot(res_large_per_aneurysm, filetxt = trans_method)
forest_plot(res_scd_in_lvaa, filetxt = trans_method)
forest_plot(res_scd_per_aneurysm, filetxt = trans_method)
# forest(res_nscd_big)
forest_plot(res_nscd_size, filetxt = trans_method)


#####################
# fit model directly

modf <- lme4::glmer(cbind(aneurysm , cohort-aneurysm) ~ 1 + (1|study), family="binomial", data = dat_raw)
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

###########
# ggplots #
###########

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
