#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al: https://www.medrxiv.org/content/10.1101/2021.05.26.21257519v1
#
# Overlayed trajectories plot
#

library(tidyverse)
library(patchwork)
library(splines)
library(sitar)

overlayed_plot_fun <- function(data1, data2, data3){ 
  ggplot(plot_d(data1), aes(.x, .y)) + geom_line(aes(col = "SITAR Model")) + theme_classic() +
    geom_line(aes(x = age, y = tblh_bmc_, col = "Linear Spline LME Model"), data = data2) + 
    geom_line(aes(x = age, y = tblh_bmc_, col = "Natural Cubic Spline LME Model"), data = data3) +
    scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) +
    scale_y_continuous(breaks=seq(500,3500,1000), limits=c(500, 3500)) +
    labs(x = 'Age - years', y = 'BMC - grams') + theme(
      legend.title = element_blank(), legend.direction = "vertical") +
    guides(colour = guide_legend(override.aes = list(size = 1))) 
} 

overlayed_plot_alsp_f <- overlayed_plot_fun(
  data1 = sitar_alsp_f_best, data2 = alsp_ls_f_pred, 
  data3 = alsp_ns_f_pred) + ggtitle("(a) ALSPAC Females")

overlayed_plot_alsp_m <- overlayed_plot_fun(
  data1 = sitar_alsp_m_best, data2 = alsp_ls_m_pred, 
  data3 = alsp_ns_m_pred) + ggtitle("(b) ALSPAC Males")

overlayed_plot_bmdcs_f <- overlayed_plot_fun(
  data1 = sitar_bmdcs_f_best, data2 = bmdcs_ls_f_pred, 
  data3 = bmdcs_ns_f_pred) + ggtitle("(c) BMDCS Females")

overlayed_plot_bmdcs_m <- overlayed_plot_fun(
  data1 = sitar_bmdcs_m_best, data2 = bmdcs_ls_m_pred, 
  data3 = bmdcs_ns_m_pred) + ggtitle("(d) BMDCS Males")

overlayed_plot_pbmas_f <- overlayed_plot_fun(
  data1 = sitar_pbmas_f_best, data2 = pbmas_ls_f_pred, 
  data3 = pbmas_ns_f_pred) + ggtitle("(e) PBMAS Females")

overlayed_plot_pbmas_m <- overlayed_plot_fun(
  data1 = sitar_pbmas_m_best, data2 = pbmas_ls_m_pred, 
  data3 = pbmas_ns_m_pred) + ggtitle("(f) PBMAS Males")

graphics.off()
pdf("results/Fig6_overlayed_traj.pdf",
    width = 5, height = 6.5)
((overlayed_plot_alsp_f | overlayed_plot_alsp_m) / 
    (overlayed_plot_bmdcs_f | overlayed_plot_bmdcs_m) / 
    (overlayed_plot_pbmas_f | overlayed_plot_pbmas_m)) + 
  plot_layout(guides = "collect") & theme(legend.title = element_blank(), legend.position = 'bottom')
dev.off()

rm(overlayed_plot_fun, overlayed_plot_alsp_f, overlayed_plot_alsp_m, overlayed_plot_bmdcs_f,
   overlayed_plot_bmdcs_m, overlayed_plot_pbmas_f, overlayed_plot_pbmas_m)
