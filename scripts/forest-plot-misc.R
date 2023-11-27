# some sample code for making prettier forest plots
# from Anthony Hatswell

p_dialysis$consolidated <-
  mcmc_areas(
    temp3,
    point_est = "mean",
    area_method = c("equal height"),
    
    pars = c(LETTERS[1:nrow(d_dialysis)], "fma", "rma", "bpp"),
    
    prob_outer = 0.95
  ) +
  
  scale_x_continuous(
    limits = c(0, 1.30),
    breaks = seq(0, 1, 0.1),
    minor_breaks = seq(0, 1, .05),
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(labels = rev(
    c(
      d_dialysis_names,
      "Fixed Effects",
      "Random Effects",
      "Bayesian Power Prior"
    )
  ),
  limits = rev,
  expand = expansion(add = c(0.5, 1.1))) +
  
  theme_light() +
  
  ggtitle("Hemodialysis") +
  
  xlab("Health State Utility Value") +
  
  theme(
    axis.text.x = element_text(
      size = 8,
      color = "black",
      family = "sans",
      margin = margin(
        t = 0,
        r = 0,
        b = 5,
        l = 0
      )
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black",
      family = "sans",
      
      face = c(
        rep('bold', 3),
        rep('plain', nrow(d_dr) - r_dr$selected),
        'italic',
        rep('plain', r_dr$selected - 1)
      )
    )
  ) +
  
  theme(panel.spacing = unit(1, "lines")) +
  
  
  
  annotate(
    geom = "text",
    
    label = "FE Weight:",
    
    x = 1.05,
    
    y = 2,
    
    size = 3,
    
    fontface = 'italic',
    
    hjust = "center",
    
    angle = 90
  ) +
  
  annotate(
    geom = "label",
    
    label = paste0(sprintf("%.1f", maweights[, 1]), "%"),
    
    x = 1.05,
    
    y = c((nrow(maweights) + 3):4),
    
    size = 2,
    
    hjust = "center"
  ) +
  
  annotate(
    geom = "text",
    
    label = "RE Weight:",
    
    x = 1.15,
    
    y = 2,
    
    size = 3,
    
    fontface = 'italic',
    
    hjust = "center",
    
    angle = 90
  ) +
  
  annotate(
    geom = "label",
    
    label =  paste0(sprintf("%.1f", maweights[, 2]), "%"),
    
    x = 1.15,
    
    y = c((nrow(maweights) + 3):4),
    
    size = 2,
    
    hjust = "center"
  ) +
  
  annotate(
    geom = "text",
    
    label = "BPP Weight:",
    
    x = 1.25,
    
    y = 2,
    
    size = 3,
    
    fontface = 'italic',
    
    hjust = "center",
    
    angle = 90
  ) +
  
  annotate(
    geom = "label",
    
    label =  paste0(sprintf("%.1f", d_dialysis[, 3])),
    
    x = 1.25,
    
    y = c((nrow(maweights) + 3):4),
    
    size = 2,
    
    hjust = "center"
  )

ggsave(
  file = paste0(s_path, "/dialysis/consolidated.png"),
  p_dialysis$consolidated,
  height = s_pheight,
  width = s_pwidth,
  units = "px",
  device = "png"
)