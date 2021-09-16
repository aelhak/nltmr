#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al: https://www.medrxiv.org/content/10.1101/2021.05.26.21257519v1
#
# GAMM EXAMPLE IN PBMAS COHORT
#

library(mgcv)
library(gamm4)
library(gratia)
library(splines)
library(tidyverse)
library(patchwork)

#### PBMAS FEMALES: ####

gamm_pbmas_f <- gamm4(
  tblh_bmc_ ~ s(age, k = 8), random = ~(ns(age, df = 4)|id),
  data = pbmas_f)

(gamm_pbmas_f_pred <- as.data.frame(
  predict(gamm_pbmas_f$gam, data.frame(
    age = seq(min(pbmas_f$age), max(pbmas_f$age), length = 100)), se.fit = TRUE)))
gamm_pbmas_f_pred$age <- seq(min(pbmas_f$age), max(pbmas_f$age), length = 100)
gamm_pbmas_f_pred

#### PBMAS MALES: ####

gamm_pbmas_m <- gamm4(
  tblh_bmc_ ~ s(age, k = 8), random = ~(ns(age, df = 4)|id),
  data = pbmas_m)

(gamm_pbmas_m_pred <- as.data.frame(
  predict(gamm_pbmas_m$gam, data.frame(
    age = seq(min(pbmas_m$age), max(pbmas_m$age), length = 100)), se.fit = TRUE)))
gamm_pbmas_m_pred$age <- seq(min(pbmas_m$age), max(pbmas_m$age), length = 100)
gamm_pbmas_m_pred

# TRAJECTORY CURVES
gamm_curve_pbmas <- ggplot(
  data = gamm_pbmas_m_pred, aes(x = age, y = fit, ymin = fit-(1.96*se.fit), ymax = fit+(1.96*se.fit))) +
  geom_line(aes(col = "Males"), data = gamm_pbmas_m_pred) + 
  geom_line(aes(col = "Females"), data = gamm_pbmas_f_pred) +
  geom_ribbon(data = gamm_pbmas_m_pred, alpha=0.1) + 
  geom_ribbon(data = gamm_pbmas_f_pred, alpha=0.1) + theme_classic() +
  scale_colour_manual(name="legend", values=c("red2", "black")) + 
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.82),
        legend.background = element_blank(),
        legend.box.background = element_rect(
          fill = "white", colour = "black")) + 
  scale_y_continuous(breaks=seq(500,3500,500), limits=c(500, 3500)) +
  labs(x = 'Age - years', y = 'BMC - grams') + ggtitle("(a) GAMM BMC Trajectory in PBMAS") 

# VELOCITY CURVES
gamm_pbmas_m_d1 <- derivatives(gamm_pbmas_m$gam, order = 1, type = "central")
gamm_pbmas_f_d1 <- derivatives(gamm_pbmas_f$gam, order = 1, type = "central")

gamm_pbmas_vel_plot <- ggplot(
  data = gamm_pbmas_m_d1, aes(x = data, y = derivative)) +
  geom_line(aes(col = "Males"), data = gamm_pbmas_m_d1) + 
  geom_line(aes(col = "Females"), data = gamm_pbmas_f_d1) + 
  scale_colour_manual(name="legend", values=c("red2", "black")) + 
  geom_vline(xintercept = 13.703975, linetype="dotted", color = "black") +
  geom_vline(xintercept = 12.218093, linetype="dotted", color = "red2") + theme_classic() +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  theme(legend.title = element_blank(), legend.position = c(0.82, 0.82),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "black")) + 
  scale_y_continuous(breaks=seq(0,350,50), limits=c(-10, 350)) +
  labs(x = 'Age - years', y = 'BMC velocity - grams per year') +
  ggtitle("(b) Estimated BMC Velocity") +
  annotate("text", x = 15.35, y= -7.3, label= "13.7y", col = "black") + 
  annotate("text", x = 10.4, y= -7.3, label = "12.2y", col = "red2")

#### FIGURE 8: MEAN TRAJECTORY & VELOCITY PLOTS ####

graphics.off()
pdf("results/Fig8_gamm.pdf", width = 5, height = 6)

gamm_curve_pbmas / gamm_pbmas_vel_plot

dev.off()
