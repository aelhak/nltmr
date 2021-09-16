#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al: https://www.medrxiv.org/content/10.1101/2021.05.26.21257519v1
#
# NATURAL CUBIC SPLINE LME MODELS
#

library(lme4)
library(spluti) # https://github.com/ZheyuanLi/spluti
library(splines)
library(tidyverse)
library(patchwork)

# ALSPAC FEMALES ####

nslme_alsp_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ ns(age, df = ", .x, ") + ( ns(age, df = 2) | id ), 
  REML = F, data = alsp_f)")))}, otherwise = NA_real_))

(nslme_alsp_f <- Filter(Negate(anyNA), nslme_alsp_f))
(nslme_alsp_f_bic <- unlist(map(nslme_alsp_f, BIC)))
(nslme_alsp_f_best <- nslme_alsp_f[[which.min(nslme_alsp_f_bic)]]) # k = 3

nslme_alsp_f_resid <- augment(nslme_alsp_f_best) %>% full_join(alsp_f)
nslme_alsp_f_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/nslme_alsp_f_resid.csv") 

(alsp_ns_f_pred <- data.frame(age = seq(min(alsp_f$age), max(alsp_f$age), length = 500)))
(alsp_ns_f_pred$tblh_bmc_ <- predict(nslme_alsp_f_best, alsp_ns_f_pred, re.form = NA))
(alsp_ns_f_mm <- model.matrix(terms(nslme_alsp_f_best), alsp_ns_f_pred))
(alsp_ns_f_pred$se <- sqrt(diag(alsp_ns_f_mm %*% vcov(nslme_alsp_f_best) %*% t(alsp_ns_f_mm))))

(alsp_ns_f_spl_mean <- BSpl2PP(alsp_f$age, fun = "ns", df = 4, coef = fixef(nslme_alsp_f_best)[-1]))
(alsp_ns_f_spl_mean <- AddInterceptPP(alsp_ns_f_spl_mean, fixef(nslme_alsp_f_best)[1]))
(alsp_ns_f_vel_mean <- as.data.frame(plot(alsp_ns_f_spl_mean, deriv = 1)))

(alsp_f_apv_mean <- solve(alsp_ns_f_spl_mean, b = 0, deriv = 2)) 
(alsp_f_pv_mean <- predict(alsp_ns_f_spl_mean, alsp_f_apv_mean, deriv = 1))
(alsp_f_apv_pv_mean <- as.data.frame(cbind(alsp_f_apv_mean, alsp_f_pv_mean)))
(alsp_f_apv_pv_mean <- alsp_f_apv_pv_mean %>% top_n(1, alsp_f_pv_mean))

rm(nslme_alsp_f, nslme_alsp_f_bic, alsp_ns_f_mm, alsp_ns_f_spl_mean, nslme_alsp_f_resid)

#### ALSPAC MALES: ####

nslme_alsp_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ ns(age, df = ", .x, ") + ( ns(age, df = 2) | id ), 
  REML = F, data = alsp_m)")))}, otherwise = NA_real_))

(nslme_alsp_m <- Filter(Negate(anyNA), nslme_alsp_m))
(nslme_alsp_m_bic <- unlist(map(nslme_alsp_m, BIC)))
(nslme_alsp_m_best <- nslme_alsp_m[[which.min(nslme_alsp_m_bic)]]) # k = 6

nslme_alsp_m_resid <- augment(nslme_alsp_m_best) %>% full_join(alsp_m)
nslme_alsp_m_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/nslme_alsp_m_resid.csv")  

(alsp_ns_m_pred <- data.frame(age = seq(min(alsp_m$age), max(alsp_m$age), length = 500)))
(alsp_ns_m_pred$tblh_bmc_ <- predict(nslme_alsp_m_best, alsp_ns_m_pred, re.form = NA))
(alsp_ns_m_mm <- model.matrix(terms(nslme_alsp_m_best), alsp_ns_m_pred))
(alsp_ns_m_pred$se <- sqrt(diag(alsp_ns_m_mm %*% vcov(nslme_alsp_m_best) %*% t(alsp_ns_m_mm))))

(alsp_ns_m_spl_mean <- BSpl2PP(alsp_m$age, fun = "ns", df = 7, coef = fixef(nslme_alsp_m_best)[-1]))
(alsp_ns_m_spl_mean <- AddInterceptPP(alsp_ns_m_spl_mean, fixef(nslme_alsp_m_best)[1]))
(alsp_ns_m_vel_mean <- as.data.frame(plot(alsp_ns_m_spl_mean, deriv = 1)))

(alsp_m_apv_mean <- solve(alsp_ns_m_spl_mean, b = 0, deriv = 2)) 
(alsp_m_pv_mean <- predict(alsp_ns_m_spl_mean, alsp_m_apv_mean, deriv = 1))
(alsp_m_apv_pv_mean <- as.data.frame(cbind(alsp_m_apv_mean, alsp_m_pv_mean)))
(alsp_m_apv_pv_mean <- alsp_m_apv_pv_mean %>% top_n(1, alsp_m_pv_mean))

rm(nslme_alsp_m, nslme_alsp_m_bic, alsp_ns_m_mm, alsp_ns_m_spl_mean, nslme_alsp_m_resid)

#### BMDCS FEMALES: ####

nslme_bmdcs_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ ns(age, df = ", .x, ") + ( ns(age, df = 2) | id ), 
  REML = F, data = bmdcs_f)")))}, otherwise = NA_real_))

(nslme_bmdcs_f <- Filter(Negate(anyNA), nslme_bmdcs_f))
(nslme_bmdcs_f_bic <- unlist(map(nslme_bmdcs_f, BIC)))
(nslme_bmdcs_f_best <- nslme_bmdcs_f[[which.min(nslme_bmdcs_f_bic)]]) # k = 5

nslme_bmdcs_f_resid <- augment(nslme_bmdcs_f_best) %>% full_join(bmdcs_f)
nslme_bmdcs_f_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/nslme_bmdcs_f_resid.csv")   

(bmdcs_ns_f_pred <- data.frame(age = seq(min(bmdcs_f$age), max(bmdcs_f$age), length = 500)))
(bmdcs_ns_f_pred$tblh_bmc_ <- predict(nslme_bmdcs_f_best, bmdcs_ns_f_pred, re.form = NA))
(bmdcs_ns_f_mm <- model.matrix(terms(nslme_bmdcs_f_best), bmdcs_ns_f_pred))
(bmdcs_ns_f_pred$se <- sqrt(diag(bmdcs_ns_f_mm %*% vcov(nslme_bmdcs_f_best) %*% t(bmdcs_ns_f_mm))))

(bmdcs_ns_f_spl_mean <- BSpl2PP(bmdcs_f$age, fun = "ns", df = 6, coef = fixef(nslme_bmdcs_f_best)[-1]))
(bmdcs_ns_f_spl_mean <- AddInterceptPP(bmdcs_ns_f_spl_mean, fixef(nslme_bmdcs_f_best)[1]))
(bmdcs_ns_f_vel_mean <- as.data.frame(plot(bmdcs_ns_f_spl_mean, deriv = 1)))

(bmdcs_f_apv_mean <- solve(bmdcs_ns_f_spl_mean, b = 0, deriv = 2)) 
(bmdcs_f_pv_mean <- predict(bmdcs_ns_f_spl_mean, bmdcs_f_apv_mean, deriv = 1))
(bmdcs_f_apv_pv_mean <- as.data.frame(cbind(bmdcs_f_apv_mean, bmdcs_f_pv_mean)))
(bmdcs_f_apv_pv_mean <- bmdcs_f_apv_pv_mean %>% top_n(1, bmdcs_f_pv_mean))

rm(nslme_bmdcs_f, nslme_bmdcs_f_bic, bmdcs_ns_f_mm, bmdcs_ns_f_spl_mean, nslme_bmdcs_f_resid)

#### BMDCS MALES: ####

nslme_bmdcs_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ ns(age, df = ", .x, ") + ( ns(age, df = 2) | id ), 
  REML = F, data = bmdcs_m)")))}, otherwise = NA_real_))

(nslme_bmdcs_m <- Filter(Negate(anyNA), nslme_bmdcs_m))
(nslme_bmdcs_m_bic <- unlist(map(nslme_bmdcs_m, BIC)))
(nslme_bmdcs_m_best <- nslme_bmdcs_m[[which.min(nslme_bmdcs_m_bic)]]) # k = 5

nslme_bmdcs_m_resid <- augment(nslme_bmdcs_m_best) %>% full_join(bmdcs_m)
nslme_bmdcs_m_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/nslme_bmdcs_m_resid.csv")  

(bmdcs_ns_m_pred <- data.frame(age = seq(min(bmdcs_m$age), max(bmdcs_m$age), length = 500)))
(bmdcs_ns_m_pred$tblh_bmc_ <- predict(nslme_bmdcs_m_best, bmdcs_ns_m_pred, re.form = NA))
(bmdcs_ns_m_mm <- model.matrix(terms(nslme_bmdcs_m_best), bmdcs_ns_m_pred))
(bmdcs_ns_m_pred$se <- sqrt(diag(bmdcs_ns_m_mm %*% vcov(nslme_bmdcs_m_best) %*% t(bmdcs_ns_m_mm))))

(bmdcs_ns_m_spl_mean <- BSpl2PP(bmdcs_m$age, fun = "ns", df = 6, coef = fixef(nslme_bmdcs_m_best)[-1]))
(bmdcs_ns_m_spl_mean <- AddInterceptPP(bmdcs_ns_m_spl_mean, fixef(nslme_bmdcs_m_best)[1]))
(bmdcs_ns_m_vel_mean <- as.data.frame(plot(bmdcs_ns_m_spl_mean, deriv = 1)))

(bmdcs_m_apv_mean <- solve(bmdcs_ns_m_spl_mean, b = 0, deriv = 2)) 
(bmdcs_m_pv_mean <- predict(bmdcs_ns_m_spl_mean, bmdcs_m_apv_mean, deriv = 1))
(bmdcs_m_apv_pv_mean <- as.data.frame(cbind(bmdcs_m_apv_mean, bmdcs_m_pv_mean)))
(bmdcs_m_apv_pv_mean <- bmdcs_m_apv_pv_mean %>% top_n(1, bmdcs_m_pv_mean))

rm(nslme_bmdcs_m, nslme_bmdcs_m_bic, bmdcs_ns_m_mm, bmdcs_ns_m_spl_mean, nslme_bmdcs_m_resid)

#### PBMAS FEMALES: ####

nslme_pbmas_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ ns(age, df = ", .x, ") + ( ns(age, df = 2) | id ), 
  REML = F, data = pbmas_f)")))}, otherwise = NA_real_))

(nslme_pbmas_f <- Filter(Negate(anyNA), nslme_pbmas_f))
(nslme_pbmas_f_bic <- unlist(map(nslme_pbmas_f, BIC)))
(nslme_pbmas_f_best <- nslme_pbmas_f[[which.min(nslme_pbmas_f_bic)]]) # k = 6

nslme_pbmas_f_resid <- augment(nslme_pbmas_f_best) %>% full_join(pbmas_f)
nslme_pbmas_f_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/nslme_pbmas_f_resid.csv")   

(pbmas_ns_f_pred <- data.frame(age = seq(min(pbmas_f$age), max(pbmas_f$age), length = 500)))
(pbmas_ns_f_pred$tblh_bmc_ <- predict(nslme_pbmas_f_best, pbmas_ns_f_pred, re.form = NA))
(pbmas_ns_f_mm <- model.matrix(terms(nslme_pbmas_f_best), pbmas_ns_f_pred))
(pbmas_ns_f_pred$se <- sqrt(diag(pbmas_ns_f_mm %*% vcov(nslme_pbmas_f_best) %*% t(pbmas_ns_f_mm))))

(pbmas_ns_f_spl_mean <- BSpl2PP(pbmas_f$age, fun = "ns", df = 7, coef = fixef(nslme_pbmas_f_best)[-1]))
(pbmas_ns_f_spl_mean <- AddInterceptPP(pbmas_ns_f_spl_mean, fixef(nslme_pbmas_f_best)[1]))
(pbmas_ns_f_vel_mean <- as.data.frame(plot(pbmas_ns_f_spl_mean, deriv = 1)))

(pbmas_f_apv_mean <- solve(pbmas_ns_f_spl_mean, b = 0, deriv = 2)) 
(pbmas_f_pv_mean <- predict(pbmas_ns_f_spl_mean, pbmas_f_apv_mean, deriv = 1))
(pbmas_f_apv_pv_mean <- as.data.frame(cbind(pbmas_f_apv_mean, pbmas_f_pv_mean)))
(pbmas_f_apv_pv_mean <- pbmas_f_apv_pv_mean %>% top_n(1, pbmas_f_pv_mean))

rm(nslme_pbmas_f, nslme_pbmas_f_bic, pbmas_ns_f_mm, pbmas_ns_f_spl_mean, nslme_pbmas_f_resid)

#### PBMAS MALES: ####

nslme_pbmas_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ ns(age, df = ", .x, ") + ( ns(age, df = 2) | id ), 
  REML = F, data = pbmas_m)")))}, otherwise = NA_real_))

(nslme_pbmas_m <- Filter(Negate(anyNA), nslme_pbmas_m))
(nslme_pbmas_m_bic <- unlist(map(nslme_pbmas_m, BIC)))
(nslme_pbmas_m_best <- nslme_pbmas_m[[which.min(nslme_pbmas_m_bic)]]) # k = 6

nslme_pbmas_m_resid <- augment(nslme_pbmas_m_best) %>% full_join(pbmas_m)
nslme_pbmas_m_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/nslme_pbmas_m_resid.csv")   

(pbmas_ns_m_pred <- data.frame(age = seq(min(pbmas_m$age), max(pbmas_m$age), length = 500)))
(pbmas_ns_m_pred$tblh_bmc_ <- predict(nslme_pbmas_m_best, pbmas_ns_m_pred, re.form = NA))
(pbmas_ns_m_mm <- model.matrix(terms(nslme_pbmas_m_best), pbmas_ns_m_pred))
(pbmas_ns_m_pred$se <- sqrt(diag(pbmas_ns_m_mm %*% vcov(nslme_pbmas_m_best) %*% t(pbmas_ns_m_mm))))

(pbmas_ns_m_spl_mean <- BSpl2PP(pbmas_m$age, fun = "ns", df = 7, coef = fixef(nslme_pbmas_m_best)[-1]))
(pbmas_ns_m_spl_mean <- AddInterceptPP(pbmas_ns_m_spl_mean, fixef(nslme_pbmas_m_best)[1]))
(pbmas_ns_m_vel_mean <- as.data.frame(plot(pbmas_ns_m_spl_mean, deriv = 1)))

(pbmas_m_apv_mean <- solve(pbmas_ns_m_spl_mean, b = 0, deriv = 2)) 
(pbmas_m_pv_mean <- predict(pbmas_ns_m_spl_mean, pbmas_m_apv_mean, deriv = 1))
(pbmas_m_apv_pv_mean <- as.data.frame(cbind(pbmas_m_apv_mean, pbmas_m_pv_mean)))
(pbmas_m_apv_pv_mean <- pbmas_m_apv_pv_mean %>% top_n(1, pbmas_m_pv_mean))

rm(nslme_pbmas_m, nslme_pbmas_m_bic, pbmas_ns_m_mm, pbmas_ns_m_spl_mean, nslme_pbmas_m_resid)

#### FIGURE 4: NATURAL SPLINE MEAN TRAJECTORY & VELOCITY PLOTS ####

nspline_plot_fun <- function(data1, data2){ ggplot(
  data = data1, aes(x = age, y = tblh_bmc_, ymin = tblh_bmc_-(1.96*se), ymax = tblh_bmc_+(1.96*se))) +
    geom_line(aes(col = "Males"), data = data1) + geom_line(aes(col = "Females"), data = data2) +
    geom_ribbon(data = data1, alpha=0.1) + geom_ribbon(data = data2, alpha=0.1) + theme_classic() +
    scale_colour_manual(name="legend", values=c("red2", "black")) + 
    scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
    theme(legend.title = element_blank(), legend.position = c(0.15, 0.82),
          legend.background = element_blank(),
          legend.box.background = element_rect(fill = "white", colour = "black")) + 
    scale_y_continuous(breaks=seq(500,3500,500), limits=c(500, 3500)) +
    labs(x = 'Age - years', y = 'BMC - grams')
} 

nspline_curve_alsp <- nspline_plot_fun(data1 = alsp_ns_m_pred, data2 = alsp_ns_f_pred) + 
  ggtitle("(a) Natural Cubic Spline BMC Trajectory: ALSPAC")

nspline_curve_bmdcs <- nspline_plot_fun(data1 = bmdcs_ns_m_pred, data2 = bmdcs_ns_f_pred) + 
  ggtitle("(c) Natural Cubic Spline BMC Trajectory: BMDCS")

nspline_curve_pbmas <- nspline_plot_fun(data1 = pbmas_ns_m_pred, data2 = pbmas_ns_f_pred) + 
  ggtitle("(e) Natural Cubic Spline BMC Trajectory: PBMAS")

nspline_vel_plot_fun <- function(data1, data2, apv1, apv2){ ggplot(
  data = data1, aes(x = x, y = y)) +
    geom_line(aes(col = "Males"), data = data1) + geom_line(aes(col = "Females"), data = data2) + 
    geom_vline(xintercept = apv1, linetype="dotted", color = "black") +
    geom_vline(xintercept = apv2, linetype="dotted", color = "red2") + theme_classic() +
    scale_colour_manual(name="legend", values=c("red2", "black")) + 
    scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
    theme(legend.title = element_blank(), legend.position = c(0.82, 0.82),
          legend.background = element_blank(),
          legend.box.background = element_rect(fill = "white", colour = "black")) + 
    scale_y_continuous(breaks=seq(0,350,50), limits=c(-10, 350)) +
    labs(x = 'Age - years', y = 'BMC velocity - grams per year')
}

nspline_vel_alsp <- nspline_vel_plot_fun(
  data1 = alsp_ns_m_vel_mean, data2 = alsp_ns_f_vel_mean, 
  apv1 = alsp_m_apv_pv_mean$alsp_m_apv_mean, apv2 = alsp_f_apv_pv_mean$alsp_f_apv_mean) + 
  ggtitle("(b) Estimated BMC Velocity: ALSPAC")

nspline_vel_bmdcs <- nspline_vel_plot_fun(
  data1 = bmdcs_ns_m_vel_mean, data2 = bmdcs_ns_f_vel_mean, 
  apv1 = bmdcs_m_apv_pv_mean$bmdcs_m_apv_mean, apv2 = bmdcs_f_apv_pv_mean$bmdcs_f_apv_mean) + 
  ggtitle("(d) Estimated BMC Velocity: BMDCS")

nspline_vel_pbmas <- nspline_vel_plot_fun(
  data1 = pbmas_ns_m_vel_mean, data2 = pbmas_ns_f_vel_mean, 
  apv1 = pbmas_m_apv_pv_mean$pbmas_m_apv_mean, apv2 = pbmas_f_apv_pv_mean$pbmas_f_apv_mean) + 
  ggtitle("(f) Estimated BMC Velocity: PBMAS")

graphics.off()
pdf("results/Fig4_nsplines.pdf",
    width = 8.6, height = 9.5)

((nspline_curve_alsp / nspline_curve_bmdcs / nspline_curve_pbmas) |
  (nspline_vel_alsp / nspline_vel_bmdcs / nspline_vel_pbmas)) + plot_layout(
    widths = c(1.2, 1)) 

dev.off()

rm(nspline_plot_fun, nspline_curve_alsp, nspline_curve_bmdcs, nspline_curve_pbmas,
   nspline_vel_plot_fun, nspline_vel_alsp, nspline_vel_bmdcs, nspline_vel_pbmas, 
   alsp_ns_f_vel_mean, alsp_ns_m_vel_mean, bmdcs_ns_f_vel_mean, bmdcs_ns_m_vel_mean, 
   pbmas_ns_f_vel_mean, pbmas_ns_m_vel_mean)
