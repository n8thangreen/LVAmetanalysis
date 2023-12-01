
# forestplot package
# version of plots


base_data <- tibble::tibble(mean  = exp(res_scd$TE),
                            lower = res_scd$lower,
                            upper = res_scd$upper,
                            study = res_scd$studlab,
                            event = res_scd$event,
                            prop = round(exp(res_scd$TE),2),
                            CI = paste0("[", round(res_scd$lower,2), ", ", round(res_scd$upper,2), "]"),
                            n = res_scd$n) |> 
  mutate(CI = ifelse(CI == "[NA, NA]", "", CI))
         
library(forestplot)
library(ggplot2)

png("plots/forestplot_version.png", height = 700, width = 750)
base_data |>
  forestplot(labeltext = c(study, event, n, prop, CI),
             boxsize = 0.2,
             clip = c(0, 0.025)) |> 
  # xlog = TRUE) |> 
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(study = c("Study"),
                event = c("Events"),
                n = c("Total"),
                prop = c("Proportion"),
                CI = c("95% CI")) |>
  fp_append_row(mean  = exp(res_scd$TE.common),
                lower = exp(res_scd$lower.common),
                upper = exp(res_scd$upper.common),
                study = "Common effect model",
                is.summary = TRUE) |>
  fp_append_row(mean  = exp(res_scd$TE.random),
                lower = exp(res_scd$lower.random),
                upper = exp(res_scd$upper.random),
                study = "Random effect model",
                is.summary = TRUE) |>
  fp_decorate_graph(graph.pos = 4) |>
  fp_set_zebra_style("#EFEFEF") |> 
  fp_decorate_graph(grid = structure(c(exp(res_scd$TE.common),
                                       exp(res_scd$TE.random)), 
                                     gp = gpar(lty = 2, col = "#CCCCFF")))

dev.off()
