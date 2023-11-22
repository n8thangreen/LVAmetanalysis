
# forest plots combining two test uncerainty and basic model
#
# https://mjskay.github.io/ggdist/articles/slabinterval.html


library(dplyr)

plot_dat_pooled <- read.csv("data/plot_dat_pooled.csv")
plot_dat_main <- read.csv("data/plot_dat_main.csv")

plot_dat_pooled <- 
  plot_dat_pooled |> 
  select(study, b_Intercept) |> 
  rename(psi = b_Intercept) |> 
  mutate(model = "Pooled")

plot_dat_main <-
  plot_dat_main |> 
  select(study, psi) |> 
  mutate(model = "Two test")

plot_dat <- 
  rbind(plot_dat_main, plot_dat_pooled)

#######
# plot

plot_dat %>%   
  ggplot(aes(x = psi, y = study,  group = model, fill = model)) +
  # geom_vline(xintercept = mean(res_jags[[1]][,"psi0"]), linewidth = 0.25, lty = 2) +
  # ggdist::stat_halfeye(.width = c(0.8, 0.95), position = "dodge") +
  stat_gradientinterval(position = "dodge") +
  xlab("Prevalence") +
  xlim(0, 0.2) +
  scale_x_continuous(labels = scales::percent)
# # Add text labels
# geom_text(
# data = mutate_if(out_all_sum, is.numeric, round, 2),
# aes(label = str_glue("{b_Intercept} [{.lower}, {.upper}]"), x = 0.2),
# hjust = "inward") 

ggsave(filename = "plots/forest_plot_combined.png",
       width = 20, height = 20, units = "cm", dpi = 640, device = "png")

