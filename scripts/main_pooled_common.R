# frequentist of aneurysm (LVA) data
# for outcomes:
#  stroke, LV thrombus, SCD, imaging, size
#
# using Freeman-Tukey (PFT) double arcsine transformation
# instead of logit transformation
#
# COMMOM EFFECT models for supplementary material
#
# confidence intervals for individual study results
# Clopper-Pearson interval ('exact' binomial interval)
# method.ci = "CP"

library(meta)
library(dplyr)
library(lme4)
library(tidyr)
library(stringr)

filename <- "LVA and outcomes_size 2024 01 28.dta"

dat_raw <- foreign::read.dta(here::here(glue::glue("data/{filename}")))

dat_raw$nsvt_aneu_n <- as.numeric(dat_raw$nsvt_aneu_n)

# remove duplicate rows
dat_raw <- dat_raw[!duplicated(dat_raw$study), ]

# reformat study names for paper
desired_ids <- c(17, 9, 18, 8, 16, 20, 15, 14, 7, 19) 

dat_raw <- dat_raw %>%
  mutate(
    # Extract author name (everything before the last space)
    author = str_extract(study, "^.*(?=\\s\\d{4})"),
    # Extract year (the four digits at the end)
    year = str_extract(study, "\\d{4}$"),
    # Add the predefined IDs
    id = desired_ids, # Use this if you have predefined IDs
    # Format the new reference string
    study = paste0(author, " (", year, ") [", id, "]")
  ) |> 
  select(-author)

##############
# frequentist
##############

## prevalence in _total_ cohort

trans_method <- "PFT"
# trans_method <- "PLOGIT"

# remove studies with NAs
res_aneurysm <-
  dat_raw[!is.na(dat_raw$aneurysm), ] |> 
  metaprop(event = aneurysm, n = cohort, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_scd <-
  dat_raw[!is.na(dat_raw$nscd), ] |> 
  metaprop(event = nscd, n = cohort, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_imaging <-
  dat_raw[!is.na(dat_raw$aneurysm), ] |> 
  metaprop(event = aneurysm, n = cohort, studlab = study,
           sm = trans_method, subgroup = imaging, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_small <-
  dat_raw[!is.na(dat_raw$n_small), ] |> 
  metaprop(event = n_small, n = cohort, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_medium <-
  dat_raw[!is.na(dat_raw$n_medium), ] |> 
  metaprop(event = n_medium, n = cohort, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_large <-
  dat_raw[!is.na(dat_raw$n_large), ] |> 
  metaprop(event = n_large, n = cohort, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_size <-
  dat_raw |> 
  reshape2:::melt.data.frame(measure.vars = c("n_small", "n_medium", "n_large"),
                             variable.name = "size") |> 
  metaprop(event = value, n = cohort, studlab = study,
           sm = trans_method, subgroup = size, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

resbind_size <-
  metabind(res_small, res_medium, res_large,
           outclab = "",
           common = TRUE,
           random = FALSE,
           backtransf = FALSE)

## prevalence in _aneurysm_ group

res_stroke <-
  dat_raw[!is.na(dat_raw$ncva), ] |> 
  metaprop(event = ncva, n = aneurysm, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_lvthrombus <-
  dat_raw[!is.na(dat_raw$nlvthrombus), ] |> 
  metaprop(event = nlvthrombus, n = aneurysm, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_svt_aneu <-
  dat_raw[!is.na(dat_raw$nsvt_aneu_n), ] |> 
  metaprop(event = nsvt_aneu_n, n = aneurysm, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_scd_in_lvaa <-
  dat_raw[!is.na(dat_raw$nscd), ] |> 
  metaprop(event = nscd, n = aneurysm, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)           

res_scd_per_aneurysm <-
  dat_raw[!is.na(dat_raw$nscd), ] |> 
  metaprop(event = nscd, n = aneurysm, studlab = study,
           sm = trans_method, subgroup = atpy_n, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_small_per_aneurysm <-
  dat_raw[!is.na(dat_raw$n_small), ] |> 
  metaprop(event = n_small, n = aneurysm, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_medium_per_aneurysm <-
  dat_raw[!is.na(dat_raw$n_medium), ] |> 
  metaprop(event = n_medium, n = aneurysm, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

res_large_per_aneurysm <-
  dat_raw[!is.na(dat_raw$n_large), ] |> 
  metaprop(event = n_large, n = aneurysm, studlab = study,
           sm = trans_method, method.tau = "REML", 
           backtransf = TRUE,  # proportions
           common = TRUE,
           random = FALSE,
           method.ci = "CP",   # exact binomial confidence intervals
           data = _)

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

dat_scdsize <- dat_scdsize |> 
  merge(dat_size, by = c("study", "size_label")) |> 
  rename(nscd = value.x, n = value.y) |> 
  filter(!is.na(nscd))

non_zero_studies <- dat_raw$n_big != 0 & !is.na(dat_raw$n_big)
res_nscd_big <-
  metaprop(event = nscd_big, n = n_big, studlab = study,
           data = dat_raw[non_zero_studies, ],
           sm = trans_method,
           method.tau = "REML",
           common = TRUE,
           random = FALSE,
           backtransf = TRUE,  # proportions
           method.ci = "CP"    # exact binomial confidence intervals
  )

non_zero_studies <- dat_scdsize$n != 0 & !is.na(dat_scdsize$n)
res_nscd_size <-
  metaprop(event = nscd, n = n, studlab = study, subgroup = size_label,
           data = dat_scdsize[non_zero_studies, ],
           sm = trans_method,
           method.tau = "REML",
           common = TRUE,
           random = FALSE,
           backtransf = TRUE,  # proportions
           method.ci = "CP"    # exact binomial confidence intervals
  )

#########
# plots #
#########

# custom plot
forest_plot <- function(x,
                        save = TRUE,
                        filetxt = "",
                        colvars = c("effect", "ci", "w.random"),
                        rhs_text = "Treatment",
                        lhs_text = "Control", ...) {
  
  if (save) {
    var_name <- deparse(substitute(x)) 
    # png(glue::glue("plots/{var_name}{filetxt}_common.png"), height = 500, width = 650)     # standard resolution
    png(glue::glue("plots/{var_name}{filetxt}_common.png"), height = 1400, width = 2500, res = 300)   # high resolution
    on.exit(dev.off())
  }  

  x$Var <- x$seTE^2
  
  # pooled on natural scale
  sd <- backtrans_delta_PFT(x$TE.common, x$tau2)^0.5
  text.addline1 <- paste0("\u03C3 = ", sprintf("%.4f", sd))
  
  if (x$sm == "OR") {
    weight <- 1 / x$Var
  } else {
    weight <- x$n / max(x$n)  # linear
  }
  
  ## Clopper-Pearson for pooled rate
  alpha <- 0.05
  total_successes <- sum(x$event)
  total_trials <- sum(x$n)
  pooled_prop <- total_successes / total_trials
  
  # Clopper-Pearson confidence interval for pooled rate/proportion
  ci <- binom::binom.confint(total_successes, total_trials, method = "exact")
  
  TE.common <- ci["mean"]
  lower.common <- ci["lower"]
  upper.common <- ci["upper"]
  
  # back-transform
  # proportions to arcsine transformed proportions
  x$TE.common <- p2asin(TE.common)
  x$lower.common <- p2asin(lower.common)
  x$upper.common <- p2asin(upper.common)
  
  x$w.common <- weight
  
  meta::forest(
    x,
    weight.study = "common",
    label.left = glue::glue("Favours {lhs_text}"),
    label.right = glue::glue("Favours {rhs_text}"),
    rightcols = colvars,
    text.addline1 = text.addline1,
    ...)
}

# transform from linear scale to rate
backtrans_delta_PFT <- function(ft_value, var_ft) {
  g <- function(x) sin(x / 2)^2        # Back-transform to proportion
  g_prime <- function(x) sin(x) / 2    # Derivative of g(x)
  
  # Variance on original scale
  (g_prime(ft_value)^2) * var_ft
}

forest_plot(res_aneurysm, filetxt = trans_method)
forest_plot(res_scd_in_lvaa, filetxt = trans_method)
forest_plot(res_stroke, filetxt = trans_method)
forest_plot(res_lvthrombus, filetxt = trans_method)
forest_plot(res_svt_aneu, filetxt = trans_method)
forest_plot(res_scd, filetxt = trans_method)
forest_plot(res_imaging, filetxt = trans_method)
forest_plot(res_small, filetxt = trans_method)
forest_plot(res_medium, filetxt = trans_method)
forest_plot(res_large, filetxt = trans_method)

forest_plot(res_small_per_aneurysm, filetxt = trans_method)
forest_plot(res_medium_per_aneurysm, filetxt = trans_method)
forest_plot(res_large_per_aneurysm, filetxt = trans_method)
forest_plot(res_scd_in_lvaa, filetxt = trans_method)
forest_plot(res_scd_per_aneurysm, filetxt = trans_method)
# forest(res_nscd_big)
forest_plot(res_nscd_size, filetxt = trans_method)

