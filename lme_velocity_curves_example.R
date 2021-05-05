#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al.
#
# EFFECT OF INCREASING RE KNOTS ON INDIVUDAL VELOCITY CURVES in LME MODEL
#

library(lme4)
library(spluti)
library(splines)
library(tidyverse)
library(patchwork)

nslme_re_ex <- map(c(3, 5, 7), possibly(~ { eval(parse(text = paste0(
   "lmer(tblh_bmc_ ~ ns(age, df = 7) + ( ns(age, df = ", .x, ") | id ), data = pbmas_m)"
   )))}, otherwise = NA_real_))

# k = 2

(nslme_re_ex_k2_mean <- BSpl2PP(pbmas_m$age, fun = "ns", df = 7, coef = fixef(nslme_re_ex[[1]])[-1]))
(nslme_re_ex_k2_mean <- AddInterceptPP(nslme_re_ex_k2_mean, fixef(nslme_re_ex[[1]])[1]))

(nslme_re_ex_k2_re <- as.matrix(ranef(nslme_re_ex[[1]])$id))
(nslme_re_ex_k2_spl_random <- BSpl2PP(pbmas_m$age, fun = "ns", df = 3, coef = t(nslme_re_ex_k2_re[, -1, drop = FALSE])))
(nslme_re_ex_k2_spl_random <- AddInterceptPP(nslme_re_ex_k2_spl_random, nslme_re_ex_k2_re[, 1]))
(nslme_re_ex_k2_spl_subject <- AddPP(nslme_re_ex_k2_mean, nslme_re_ex_k2_spl_random))

# k = 4

(nslme_re_ex_k4_mean <- BSpl2PP(pbmas_m$age, fun = "ns", df = 7, coef = fixef(nslme_re_ex[[2]])[-1]))
(nslme_re_ex_k4_mean <- AddInterceptPP(nslme_re_ex_k4_mean, fixef(nslme_re_ex[[2]])[1]))

(nslme_re_ex_k4_re <- as.matrix(ranef(nslme_re_ex[[2]])$id))
(nslme_re_ex_k4_spl_random <- BSpl2PP(pbmas_m$age, fun = "ns", df = 5, coef = t(nslme_re_ex_k4_re[, -1, drop = FALSE])))
(nslme_re_ex_k4_spl_random <- AddInterceptPP(nslme_re_ex_k4_spl_random, nslme_re_ex_k4_re[, 1]))
(nslme_re_ex_k4_spl_subject <- AddPP(nslme_re_ex_k4_mean, nslme_re_ex_k4_spl_random))

# k = 6

(nslme_re_ex_k6_mean <- BSpl2PP(pbmas_m$age, fun = "ns", df = 7, coef = fixef(nslme_re_ex[[3]])[-1]))
(nslme_re_ex_k6_mean <- AddInterceptPP(nslme_re_ex_k6_mean, fixef(nslme_re_ex[[3]])[1]))

(nslme_re_ex_k6_re <- as.matrix(ranef(nslme_re_ex[[3]])$id))
(nslme_re_ex_k6_spl_random <- BSpl2PP(pbmas_m$age, fun = "ns", df = 7, coef = t(nslme_re_ex_k6_re[, -1, drop = FALSE])))
(nslme_re_ex_k6_spl_random <- AddInterceptPP(nslme_re_ex_k6_spl_random, nslme_re_ex_k6_re[, 1]))
(nslme_re_ex_k6_spl_subject <- AddPP(nslme_re_ex_k6_mean, nslme_re_ex_k6_spl_random))

#### FIGURE 8 ####

(k2_re_plot <- as.data.frame(plot(nslme_re_ex_k2_spl_subject, deriv = 1)) %>% select(x, y.50:y.54) %>% 
   pivot_longer(cols = starts_with("y."), names_to = "id", values_to = "y") %>% mutate(
     knots = "2 knots"))

(k4_re_plot <- as.data.frame(plot(nslme_re_ex_k4_spl_subject, deriv = 1)) %>% select(x, y.50:y.54) %>% 
    pivot_longer(cols = starts_with("y."), names_to = "id", values_to = "y") %>% mutate(
      knots = "4 knots"))

(k6_re_plot <- as.data.frame(plot(nslme_re_ex_k6_spl_subject, deriv = 1)) %>% select(x, y.50:y.54) %>% 
    pivot_longer(cols = starts_with("y."), names_to = "id", values_to = "y") %>% mutate(
      knots = "6 knots"))

graphics.off()
pdf("results/Fig8_lme_velocity.pdf", width = 3, height = 3.5)

k2_re_plot %>% full_join(k4_re_plot) %>% full_join(k6_re_plot) %>% ggplot(
  data = ., aes(x = x, y = y, col = id)) + geom_line(size = 0.5) + theme_bw() + 
  facet_grid(knots ~ .) + scale_y_continuous(breaks = seq(0,600,200), limits = c(-50,500)) +
  theme(legend.position = "none") + labs(x = 'Age - years', y = 'BMC growth velocity - grams per year')

dev.off()