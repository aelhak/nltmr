#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al
#
# LINEAR SPLINE LME MODELS
#

library(lme4)
library(lspline)
library(tidyverse)
library(patchwork)
library(broom.mixed)

#### ALSPAC FEMALES ####

lslme_alsp_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ elspline(age, ", .x, ") + (elspline(age, 2) | id), 
  REML = F, data = alsp_f)")))}, otherwise = NA_real_))

(lslme_alsp_f <- Filter(Negate(anyNA), lslme_alsp_f))
(lslme_alsp_f_bic <- unlist(map(lslme_alsp_f, BIC)))
(lslme_alsp_f_best <- lslme_alsp_f[[which.min(lslme_alsp_f_bic)]]) # k = 6

lslme_alsp_f_resid <- augment(lslme_alsp_f_best) %>% full_join(alsp_f)
lslme_alsp_f_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/lslme_alsp_f_resid.csv") 

print(c(min(alsp_f$age), attr(elspline(alsp_f$age, 7), "knots"), max(alsp_f$age)))
as.data.frame(tidy(lslme_alsp_f_best, conf.int = T))

(alsp_ls_f_pred <- data.frame(age = seq(min(alsp_f$age), max(alsp_f$age), length = 50)))
(alsp_ls_f_pred$tblh_bmc_ <- predict(lslme_alsp_f_best, alsp_ls_f_pred, re.form = NA))
(alsp_ls_f_mm <- model.matrix(terms(lslme_alsp_f_best), alsp_ls_f_pred))
(alsp_ls_f_pred$se <- sqrt(diag(alsp_ls_f_mm %*% vcov(lslme_alsp_f_best) %*% t(alsp_ls_f_mm))))

rm(lslme_alsp_f, lslme_alsp_f_bic, alsp_ls_f_mm, lslme_alsp_f_resid)

#### ALSPAC MALES ####

lslme_alsp_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ elspline(age, ", .x, ") + (elspline(age, 2) | id), 
  REML = F, data = alsp_m)")))}, otherwise = NA_real_))

(lslme_alsp_m <- Filter(Negate(anyNA), lslme_alsp_m))
(lslme_alsp_m_bic <- unlist(map(lslme_alsp_m, BIC)))
(lslme_alsp_m_best <- lslme_alsp_m[[which.min(lslme_alsp_m_bic)]]) # k = 4

lslme_alsp_m_resid <- augment(lslme_alsp_m_best) %>% full_join(alsp_m)
lslme_alsp_m_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/lslme_alsp_m_resid.csv") 

print(c(min(alsp_m$age), attr(elspline(alsp_m$age, 5), "knots"), max(alsp_m$age)))
as.data.frame(tidy(lslme_alsp_m_best, conf.int = T))

(alsp_ls_m_pred <- data.frame(age = seq(min(alsp_m$age), max(alsp_m$age), length = 50)))
(alsp_ls_m_pred$tblh_bmc_ <- predict(lslme_alsp_m_best, alsp_ls_m_pred, re.form = NA))
(alsp_ls_m_mm <- model.matrix(terms(lslme_alsp_m_best), alsp_ls_m_pred))
(alsp_ls_m_pred$se <- sqrt(diag(alsp_ls_m_mm %*% vcov(lslme_alsp_m_best) %*% t(alsp_ls_m_mm))))

rm(lslme_alsp_m, lslme_alsp_m_bic, alsp_ls_m_mm, lslme_alsp_m_resid)

#### BMDCS FEMALES ####

lslme_bmdcs_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ elspline(age, ", .x, ") + (elspline(age, 2) | id), 
  REML = F, data = bmdcs_f)")))}, otherwise = NA_real_))

(lslme_bmdcs_f <- Filter(Negate(anyNA), lslme_bmdcs_f))
(lslme_bmdcs_f_bic <- unlist(map(lslme_bmdcs_f, BIC)))
(lslme_bmdcs_f_best <- lslme_bmdcs_f[[which.min(lslme_bmdcs_f_bic)]]) # k = 5

lslme_bmdcs_f_resid <- augment(lslme_bmdcs_f_best) %>% full_join(bmdcs_f)
lslme_bmdcs_f_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/lslme_bmdcs_f_resid.csv")  

print(c(min(bmdcs_f$age), attr(elspline(bmdcs_f$age, 6), "knots"), max(bmdcs_f$age)))
as.data.frame(tidy(lslme_bmdcs_f_best, conf.int = T))

(bmdcs_ls_f_pred <- data.frame(age = seq(min(bmdcs_f$age), max(bmdcs_f$age), length = 50)))
(bmdcs_ls_f_pred$tblh_bmc_ <- predict(lslme_bmdcs_f_best, bmdcs_ls_f_pred, re.form = NA))
(bmdcs_ls_f_mm <- model.matrix(terms(lslme_bmdcs_f_best), bmdcs_ls_f_pred))
(bmdcs_ls_f_pred$se <- sqrt(diag(bmdcs_ls_f_mm %*% vcov(lslme_bmdcs_f_best) %*% t(bmdcs_ls_f_mm))))

rm(lslme_bmdcs_f, lslme_bmdcs_f_bic, bmdcs_ls_f_mm, lslme_bmdcs_f_resid)

#### BMDCS MALES ####

lslme_bmdcs_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ elspline(age, ", .x, ") + (elspline(age, 2) | id), 
  REML = F, data = bmdcs_m)")))}, otherwise = NA_real_))

(lslme_bmdcs_m <- Filter(Negate(anyNA), lslme_bmdcs_m))
(lslme_bmdcs_m_bic <- unlist(map(lslme_bmdcs_m, BIC)))
(lslme_bmdcs_m_best <- lslme_bmdcs_m[[which.min(lslme_bmdcs_m_bic)]]) # k = 4

lslme_bmdcs_m_resid <- augment(lslme_bmdcs_m_best) %>% full_join(bmdcs_m)
lslme_bmdcs_m_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/lslme_bmdcs_m_resid.csv") 

print(c(min(bmdcs_m$age), attr(elspline(bmdcs_m$age, 5), "knots"), max(bmdcs_m$age)))
as.data.frame(tidy(lslme_bmdcs_m_best, conf.int = T))

(bmdcs_ls_m_pred <- data.frame(age = seq(min(bmdcs_m$age), max(bmdcs_m$age), length = 50)))
(bmdcs_ls_m_pred$tblh_bmc_ <- predict(lslme_bmdcs_m_best, bmdcs_ls_m_pred, re.form = NA))
(bmdcs_ls_m_mm <- model.matrix(terms(lslme_bmdcs_m_best), bmdcs_ls_m_pred))
(bmdcs_ls_m_pred$se <- sqrt(diag(bmdcs_ls_m_mm %*% vcov(lslme_bmdcs_m_best) %*% t(bmdcs_ls_m_mm))))

rm(lslme_bmdcs_m, lslme_bmdcs_m_bic, bmdcs_ls_m_mm, lslme_bmdcs_m_resid)

#### PBMAS FEMALES ####

lslme_pbmas_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ elspline(age, ", .x, ") + (elspline(age, 2) | id), 
  REML = F, data = pbmas_f)")))}, otherwise = NA_real_))

(lslme_pbmas_f <- Filter(Negate(anyNA), lslme_pbmas_f))
(lslme_pbmas_f_bic <- unlist(map(lslme_pbmas_f, BIC)))
(lslme_pbmas_f_best <- lslme_pbmas_f[[which.min(lslme_pbmas_f_bic)]]) # k = 4

lslme_pbmas_f_resid <- augment(lslme_pbmas_f_best) %>% full_join(pbmas_f)
lslme_pbmas_f_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/lslme_pbmas_f_resid.csv")  

print(c(min(pbmas_f$age), attr(elspline(pbmas_f$age, 5), "knots"), max(pbmas_f$age)))
as.data.frame(tidy(lslme_pbmas_f_best, conf.int = T))

(pbmas_ls_f_pred <- data.frame(age = seq(min(pbmas_f$age), max(pbmas_f$age), length = 50)))
(pbmas_ls_f_pred$tblh_bmc_ <- predict(lslme_pbmas_f_best, pbmas_ls_f_pred, re.form = NA))
(pbmas_ls_f_mm <- model.matrix(terms(lslme_pbmas_f_best), pbmas_ls_f_pred))
(pbmas_ls_f_pred$se <- sqrt(diag(pbmas_ls_f_mm %*% vcov(lslme_pbmas_f_best) %*% t(pbmas_ls_f_mm))))

rm(lslme_pbmas_f, lslme_pbmas_f_bic, pbmas_ls_f_mm, lslme_pbmas_f_resid)

#### PBMAS MALES ####

lslme_pbmas_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "lmer(tblh_bmc_ ~ elspline(age, ", .x, ") + (elspline(age, 2) | id), 
  REML = F, data = pbmas_m)")))}, otherwise = NA_real_))

(lslme_pbmas_m <- Filter(Negate(anyNA), lslme_pbmas_m))
(lslme_pbmas_m_bic <- unlist(map(lslme_pbmas_m, BIC)))
(lslme_pbmas_m_best <- lslme_pbmas_m[[which.min(lslme_pbmas_m_bic)]]) # k = 6

lslme_pbmas_m_resid <- augment(lslme_pbmas_m_best) %>% full_join(pbmas_m)
lslme_pbmas_m_resid %>% group_by(visit) %>% select(tblh_bmc_, .fitted, .resid) %>% 
  summarise_all(list(mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1))) %>% 
  write_csv("results/lslme_pbmas_m_resid.csv")  

print(c(min(pbmas_m$age), attr(elspline(pbmas_m$age, 7), "knots"), max(pbmas_m$age)))
as.data.frame(tidy(lslme_pbmas_m_best, conf.int = T))

(pbmas_ls_m_pred <- data.frame(age = seq(min(pbmas_m$age), max(pbmas_m$age), length = 50)))
(pbmas_ls_m_pred$tblh_bmc_ <- predict(lslme_pbmas_m_best, pbmas_ls_m_pred, re.form = NA))
(pbmas_ls_m_mm <- model.matrix(terms(lslme_pbmas_m_best), pbmas_ls_m_pred))
(pbmas_ls_m_pred$se <- sqrt(diag(pbmas_ls_m_mm %*% vcov(lslme_pbmas_m_best) %*% t(pbmas_ls_m_mm))))

rm(lslme_pbmas_m, lslme_pbmas_m_bic, pbmas_ls_m_mm, lslme_pbmas_m_resid)

#### FIGURE 3 ####

lspline_plot_fun <- function(data1, data2){ ggplot(
  data = data1, aes(x = age, y = tblh_bmc_, ymin = tblh_bmc_-(1.96*se), ymax = tblh_bmc_+(1.96*se))) +
    geom_line(aes(col = "Males"), data = data1) + geom_line(aes(col = "Females"), data = data2) + 
    geom_ribbon(data = data1, alpha = 0.1) + geom_ribbon(data = data2, alpha = 0.1) + theme_bw() +
    scale_colour_manual(name = "legend", values=c("red2", "black")) +
    scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
    theme(legend.title = element_blank(), legend.position = c(0.17, 0.84),
          legend.background = element_blank(),
          legend.box.background = element_rect(fill = "white", colour = "black")) + 
    scale_y_continuous(breaks=seq(500,3500,500), limits=c(500, 3500)) +
    labs(x = 'Age - years', y = 'BMC - grams')
} 

lspline_curve_alsp <- lspline_plot_fun(data1 = alsp_ls_m_pred, data2 = alsp_ls_f_pred) + 
  ggtitle("(a) Linear Spline BMC Trajectory: ALSPAC")

lspline_curve_bmdcs <- lspline_plot_fun(data1 = bmdcs_ls_m_pred, data2 = bmdcs_ls_f_pred) + 
  ggtitle("(b) Linear Spline BMC Trajectory: BMDCS")

lspline_curve_pbmas <- lspline_plot_fun(data1 = pbmas_ls_m_pred, data2 = pbmas_ls_f_pred) + 
  ggtitle("(c) Linear Spline BMC Trajectory: PBMAS")

graphics.off()
pdf("results/Fig3_lsplines.pdf",
    width = 14, height = 4)
(lspline_curve_alsp | lspline_curve_bmdcs | lspline_curve_pbmas) 
dev.off()

rm(lspline_plot_fun, lspline_curve_alsp, lspline_curve_bmdcs, lspline_curve_pbmas)
