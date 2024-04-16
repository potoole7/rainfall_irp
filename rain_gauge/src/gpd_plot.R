#### Explore Generalised Pareto Disitribution for some parameter values ####

library(SpatialExtremes)
library(dplyr)
library(ggplot2)
library(latex2exp)

sigmas <- c(1, 2)
xis <- c(1, 5, 20)
# xi <- c(0, 1, 2)
n <- 1000

pars <- tidyr::crossing("sigma" = sigmas, "xi" = xis) %>% 
  arrange(xi, sigma)

densities <- data.frame(
  "sigma"   = rep(pars$sigma, n),
  "xi"      = rep(pars$xi, n),
  "density" = NA
) %>% 
  mutate(param = paste0("sigma = ", sigma, ", xi = ", xi)) %>% 
  arrange(xi, sigma)

for (i in seq_len(nrow(pars))) {
  # print(pars$sigma[i])
  # print(pars$xi[i])
  set.seed(123)
  densities[
    densities$sigma == pars$sigma[i] & densities$xi == pars$xi[i], "density"
  ] <- 
    # rgpd(n, 0, pars$sigma[i], pars$xi[i])
    dgpd(seq(0, 5, len = 1000), 0, pars$sigma[i], pars$xi[i])
}

p <- densities %>% 
  mutate(
    param = factor(
      param, 
      c(
        "sigma = 1, xi = 1", 
        "sigma = 2, xi = 1", 
        "sigma = 1, xi = 5", 
        "sigma = 2, xi = 5", 
        "sigma = 1, xi = 20", 
        "sigma = 2, xi = 20"
      )
    )
  ) %>% 
  ggplot(aes(
    x = rep(seq(0, 5, len = 1000), nrow(pars)), 
    y = density, 
    # colour = factor(param),
    colour = factor(param),
    # linetype = factor(param),
    linetype = factor(param),
  )) + 
  geom_line(linewidth = 1.3, alpha = 0.8) + 
  scale_x_continuous(
    expand = c(0, 0), 
    limits = c(0, 5), 
    n.breaks = 5,
  ) + 
  scale_y_continuous(
    expand = c(0, 0), 
    limits = c(0, 1), 
    n.breaks = 5,
  ) + 
  labs(
    x = "", 
    y = "Density", 
    colour = expression(~ sigma * "," ~ xi), 
    linetype = expression(~ sigma * "," ~ xi), 
    title = expression(
      "Generalised Pareto distribution density for some" ~ sigma * "," ~ xi
    )
  ) + 
  scale_color_manual(
    "", 
    # values = sort(rep(c("red", "blue", "green"), 2)),
    values = sort(rep(RColorBrewer::brewer.pal(3, "Dark2"), 2)),
    labels = rep(unname(TeX(paste0("$\\sigma =$", pars$sigma,  ", $xi =$", pars$xi))), 2)
  ) + 
  scale_linetype_manual(
    "", 
    # values = sort(rep(c(1, 4), 3)), 
    values = rep(c(1, 4), 3), 
    labels = rep(unname(TeX(paste0("$\\sigma =$", pars$sigma,  ", $xi =$", pars$xi))), 2)
  ) + 
  theme_bw(base_size = 8) + 
  # Altering plot text size
  theme(
    axis.text.x         = element_text(size = rel(1.5), colour = "black"),
    axis.text.y         = element_text(size = rel(1.5), colour = "black"),
    axis.title.y         = element_text(size = rel(1.5), colour = "black"),
    legend.text         = element_text(size = rel(1.6), colour = "black"),
    legend.position     = c(0.88, 0.85),
    # remove white box behind legend
    legend.background   = element_rect(colour = NA, fill = NA),
    legend.key          = element_rect(colour = NA, fill = NA),
    plot.title          = element_text(
      size = rel(1.6), hjust = 0.5, vjust = -1.5
    ), 
    # panel spacing and right-hand margin so x-axis labels fit & don't touch
    # panel.spacing       = unit(0.65, units = "cm"), 
    # panel.border        = element_blank(),
    # panel.grid.major.y  = element_blank(),
    # plot.margin         = unit(c(0, 0.5, 0.2, 0), "cm"),
    strip.background    = element_rect(fill = NA, colour = "white"),
    panel.background    = element_rect(fill = NA, color = "black")
  )

# dev.new(width = 6.3, height = 6, noRStudioGD = TRUE)
# p
# dev.off()
ggsave("gdp_plot.png", p, width = 6.3, height = 6)
