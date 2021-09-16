#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al: https://www.medrxiv.org/content/10.1101/2021.05.26.21257519v1
#
# LATENT TRAJECTORY MODELS (GROWTH MIXTURE MODELS)
#

library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(splines)
library(lattice)
library(lcmm)

set.seed(2020)

#### ALSPAC FEMALES ####

# 1-5 CLASS MODELS & FIT STATISTICS 

ltm_alsp_f_k1 <- hlme(
  fixed = tblh_bmc_ ~ 1 + ns(age, 4), random = ~ 1 + age,
  idiag = F, subject = 'id', ng = 1, data = alsp_f) 

ltm_alsp_f_k2k5 <- map(2:5, ~ { eval(parse(text = paste0(
    "gridsearch(hlme(fixed = tblh_bmc_ ~ 1 + ns(age, 4), random = ~ 1 + age, idiag = F, 
    subject = 'id', mixture = ~ ns(age, 4), ng = ", .x, ", nwg = T, data = alsp_f), 
    minit = ltm_alsp_f_k1, maxiter = 10, rep = 100)")))})

(ltm_alsp_f_results <- c(list(ltm_alsp_f_k1),ltm_alsp_f_k2k5))
(ltm_alsp_f_results <- ltm_alsp_f_results[-5])

map_df(ltm_alsp_f_results, ~ {return(data.frame(summarytable(.x, which = c(
  "G", "npm", "BIC", "entropy", "%class"))))})  

# PLOT TRAJECTORIES FROM 1-5 CLASS MODELS

graphics.off()
pdf("results/ltm_alsp_f_k1k4.pdf",
    width = 7, height = 7)
par(
  mfrow    = c(2,2),
  mar      = c(5, 5, 2, 2)
)

iwalk(map(ltm_alsp_f_results, ~ predictY(., data.frame(age = seq(min(alsp_f$age), max(alsp_f$age), length = 100)), var.time ="age")), ~{
  plot(.x, main = "", legend.loc = "topleft", ylim=c(500,4000), xlim=c(5,40),xlab="Age - years", ylab="BMC - grams", lty=1)
  mtext("BMC Latent Trajectories from 1 to 4 Class Models in ALSPAC Females", outer = T, line = -1.5)})

dev.off()

# PLOT TRAJECTORIES FROM SELECTED MODEL & PPROB

(ltm_alsp_f_best <- ltm_alsp_f_results[[3]])
ltm_alsp_f_pred <- data.frame(age = seq(min(alsp_f$age), max(alsp_f$age), length = 50))
ltm_alsp_f_pred <- predictY(ltm_alsp_f_best, ltm_alsp_f_pred, var.time ="age", draws = T)
(ltm_alsp_f_pred <- cbind(ltm_alsp_f_pred$times, ltm_alsp_f_pred$pred))
summarytable(ltm_alsp_f_best, which = "%class")
postprob(ltm_alsp_f_best)

ltm_alsp_f_plot <- ggplot(
  data = ltm_alsp_f_pred, aes(x = age)) +
  geom_line(aes(y = Ypred_class1, col = "Class 1 (58.3%)")) +
  geom_line(aes(y = Ypred_class2, col = "Class 2 (13.3%)")) + 
  geom_line(aes(y = Ypred_class3, col = "Class 3 (28.4%)")) + theme_classic() +
  geom_ribbon(aes(ymin = lower.Ypred_class1, ymax = upper.Ypred_class1), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class2, ymax = upper.Ypred_class2), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class3, ymax = upper.Ypred_class3), alpha=0.1) +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  scale_y_continuous(breaks=seq(500, 3500, 500), limits=c(500, 3500)) + 
  labs(x = 'Age - years', y = 'BMC - grams') + scale_color_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.20, 0.80),
        legend.background = element_blank(),
        legend.box.background = element_rect(
          fill = "white", colour = "black")) +
  ggtitle("(a) BMC Latent Trajectories: ALSPAC Females")

#### ALSPAC MALES ####

# 1-5 CLASS MODELS & FIT STATISTICS 

ltm_alsp_m_k1 <- hlme(
  fixed = tblh_bmc_ ~ 1 + ns(age, 7), random = ~ 1 + age,
  idiag = F, subject = 'id', ng = 1, data = alsp_m)

ltm_alsp_m_k2k5 <- map(2:5, ~ { eval(parse(text = paste0(
  "gridsearch(hlme(fixed = tblh_bmc_ ~ 1 + ns(age, 7), random = ~ 1 + age, idiag = F, 
    subject = 'id', mixture = ~ ns(age, 7), ng = ", .x, ", nwg = T, data = alsp_m), 
    minit = ltm_alsp_m_k1, maxiter = 10, rep = 100)")))})

(ltm_alsp_m_results <- c(list(ltm_alsp_m_k1),ltm_alsp_m_k2k5))

map_df(ltm_alsp_m_results, ~ {return(data.frame(summarytable(.x, which = c(
  "G", "npm", "BIC", "entropy", "%class"))))})  

# PLOT TRAJECTORIES FROM 1-5 CLASS MODELS

graphics.off()
pdf("results/ltm_alsp_m_k1k5.pdf",
    width = 7, height = 8)
par(
  mfrow    = c(3,2),
  mar      = c(5, 5, 2, 2)
)

iwalk(map(ltm_alsp_m_results, ~ predictY(., data.frame(age = seq(min(alsp_m$age), max(alsp_m$age), length = 100)), var.time ="age")), ~{
  plot(.x, main = "", legend.loc = "topleft", ylim=c(500,4000), xlim=c(5,40),xlab="Age - years", ylab="BMC - grams", lty=1)
  mtext("BMC Latent Trajectories from 1 to 5 Class Models in ALSPAC Males", outer = T, line = -1.5)})

dev.off()

# PLOT TRAJECTORIES FROM SELECTED MODEL & PPROB

(ltm_alsp_m_best <- ltm_alsp_m_results[[2]])
ltm_alsp_m_pred <- data.frame(age = seq(min(alsp_m$age), max(alsp_m$age), length = 50))
ltm_alsp_m_pred <- predictY(ltm_alsp_m_best, ltm_alsp_m_pred, var.time ="age", draws = T)
ltm_alsp_m_pred <- cbind(ltm_alsp_m_pred$times, ltm_alsp_m_pred$pred)
summarytable(ltm_alsp_m_best, which = "%class")
postprob(ltm_alsp_m_best)

ltm_alsp_m_plot <- ggplot(
  data = ltm_alsp_m_pred, aes(x = age)) +
  geom_line(aes(y = Ypred_class1, col = "Class 1: (34.4%)")) +
  geom_line(aes(y = Ypred_class2, col = "Class 2: (65.6%)")) + theme_classic() +
  geom_ribbon(aes(ymin = lower.Ypred_class1, ymax = upper.Ypred_class1), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class2, ymax = upper.Ypred_class2), alpha=0.1) +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  scale_y_continuous(breaks=seq(500, 3500, 500), limits=c(500, 3500)) + 
  labs(x = 'Age - years', y = 'BMC - grams') + scale_color_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.20, 0.86),
        legend.background = element_blank(),
        legend.box.background = element_rect(
          fill = "white", colour = "black")) +
  ggtitle("(b) BMC Latent Trajectories: ALSPAC Males")

#### BMDCS FEMALES ####

# 1-5 CLASS MODELS & FIT STATISTICS 

ltm_bmdcs_f_k1 <- hlme(
  fixed = tblh_bmc_ ~ 1 + ns(age, 6), random = ~ 1 + age,
  idiag = F, subject = 'id', ng = 1, data = bmdcs_f)

ltm_bmdcs_f_k2k5 <- map(2:5, ~ { eval(parse(text = paste0(
    "gridsearch(hlme(fixed = tblh_bmc_ ~ 1 + ns(age, 6),random = ~ 1 + age, idiag = F, 
    subject = 'id', mixture = ~ ns(age, 6), ng = ", .x, ", nwg = T, data = bmdcs_f), 
    minit = ltm_bmdcs_f_k1, maxiter = 10, rep = 100)")))})

(ltm_bmdcs_f_results <- c(list(ltm_bmdcs_f_k1),ltm_bmdcs_f_k2k5))

map_df(ltm_bmdcs_f_results, ~ {return(data.frame(summarytable(.x, which = c(
  "G", "npm", "BIC", "entropy", "%class"))))})  

# PLOT TRAJECTORIES FROM 1-5 CLASS MODELS

graphics.off()
pdf("results/ltm_bmdcs_f_k1k5.pdf",
    width = 7, height = 8)
par(
  mfrow    = c(3,2),
  mar      = c(5, 5, 2, 2)
)

iwalk(map(ltm_bmdcs_f_results, ~ predictY(., data.frame(age = seq(min(bmdcs_f$age), max(bmdcs_f$age), length = 100)), var.time ="age")), ~{
  plot(.x, main = "", legend.loc = "topleft", ylim=c(500,4000), xlim=c(5,40),xlab="Age - years", ylab="BMC - grams", lty=1)
  mtext("BMC Latent Trajectories from 1 to 5 Class Models in BMDCS Females", outer = T, line = -1.5)})

dev.off()

# PLOT TRAJECTORIES FROM SELECTED MODEL & PPROB

ltm_bmdcs_f_best <- ltm_bmdcs_f_results[[3]]
ltm_bmdcs_f_pred <- data.frame(age = seq(min(bmdcs_f$age), max(bmdcs_f$age), length = 50))
ltm_bmdcs_f_pred <- predictY(ltm_bmdcs_f_best, ltm_bmdcs_f_pred, var.time ="age", draws = T)
(ltm_bmdcs_f_pred <- cbind(ltm_bmdcs_f_pred$times, ltm_bmdcs_f_pred$pred))
summarytable(ltm_bmdcs_f_best, which = "%class")
postprob(ltm_bmdcs_f_best)

ltm_bmdcs_f_plot <- ggplot(
  data = ltm_bmdcs_f_pred, aes(x = age)) +
  geom_line(aes(y = Ypred_class1, col = "Class 1 (20.7%)")) +
  geom_line(aes(y = Ypred_class2, col = "Class 2 (19.1%)")) +
  geom_line(aes(y = Ypred_class3, col = "Class 3 (60.2%)")) + theme_classic() +
  geom_ribbon(aes(ymin = lower.Ypred_class1, ymax = upper.Ypred_class1), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class2, ymax = upper.Ypred_class2), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class3, ymax = upper.Ypred_class3), alpha=0.1) +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  scale_y_continuous(breaks=seq(500, 3500, 500), limits=c(500, 3500)) + 
  labs(x = 'Age - years', y = 'BMC - grams') + scale_color_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.20, 0.80),
        legend.background = element_blank(),
        legend.box.background = element_rect(
          fill = "white", colour = "black")) +
  ggtitle("(c) BMC Latent Trajectories: BMDCS Females")

#### BMDCS MALES ####

# 1-5 CLASS MODELS & FIT STATISTICS 

ltm_bmdcs_m_k1 <- hlme(
  fixed = tblh_bmc_ ~ 1 + ns(age, 6), random = ~ 1 + age,
  idiag = F, subject = 'id', ng = 1, data = bmdcs_m) 

ltm_bmdcs_m_k2k5 <- map(2:5, ~ { eval(parse(text = paste0(
    "gridsearch(hlme(fixed = tblh_bmc_ ~ 1 + ns(age, 6), random = ~ 1 + age, idiag = F, 
    subject = 'id', mixture = ~ ns(age, 6), ng = ", .x, ", nwg = T, data = bmdcs_m), 
    minit = ltm_bmdcs_m_k1, maxiter = 10, rep = 100)")))})

(ltm_bmdcs_m_results <- c(list(ltm_bmdcs_m_k1),ltm_bmdcs_m_k2k5))

map_df(ltm_bmdcs_m_results, ~ {return(data.frame(summarytable(.x, which = c(
  "G", "npm", "BIC", "entropy", "%class"))))})  

# PLOT TRAJECTORIES FROM 1-5 CLASS MODELS

graphics.off()
pdf("results/ltm_bmdcs_m_k1k5.pdf",
    width = 7, height = 8)
par(
  mfrow    = c(3,2),
  mar      = c(5, 5, 2, 2)
)

iwalk(map(ltm_bmdcs_m_results, ~ predictY(., data.frame(age = seq(min(bmdcs_m$age), max(bmdcs_m$age), length = 100)), var.time ="age")), ~{
  plot(.x, main = "", legend.loc = "topleft", ylim=c(500,4000), xlim=c(5,40),xlab="Age - years", ylab="BMC - grams", lty=1)
  mtext("BMC Latent Trajectories from 1 to 5 Class Models in BMDCS Males", outer = T, line = -1.5)})

dev.off()

# PLOT TRAJECTORIES FROM SELECTED MODEL & PPROB

ltm_bmdcs_m_best <- ltm_bmdcs_m_results[[2]]
ltm_bmdcs_m_pred <- data.frame(age = seq(min(bmdcs_m$age), max(bmdcs_m$age), length = 50))
ltm_bmdcs_m_pred <- predictY(ltm_bmdcs_m_best, ltm_bmdcs_m_pred, var.time ="age", draws = T)
(ltm_bmdcs_m_pred <- cbind(ltm_bmdcs_m_pred$times, ltm_bmdcs_m_pred$pred))
summarytable(ltm_bmdcs_m_best, which = "%class")
postprob(ltm_bmdcs_m_best)

ltm_bmdcs_m_plot <- ggplot(
  data = ltm_bmdcs_m_pred, aes(x = age)) +
  geom_line(aes(y = Ypred_class1, col = "Class 1 (37.8%)")) +
  geom_line(aes(y = Ypred_class2, col = "Class 2 (62.2%)")) + theme_classic() +
  geom_ribbon(aes(ymin = lower.Ypred_class1, ymax = upper.Ypred_class1), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class2, ymax = upper.Ypred_class2), alpha=0.1) +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  scale_y_continuous(breaks=seq(500, 3500, 500), limits=c(500, 3500)) + 
  labs(x = 'Age - years', y = 'BMC - grams') + scale_color_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.20, 0.85),
        legend.background = element_blank(),
        legend.box.background = element_rect(
          fill = "white", colour = "black")) +
  ggtitle("(d) BMC Latent Trajectories: BMDCS Males")

#### PBMAS FEMALES ####

# 1-5 CLASS MODELS & FIT STATISTICS 

ltm_pbmas_f_k1 <- hlme(
  fixed = tblh_bmc_ ~ 1 + ns(age, 7), random = ~ 1 + age,
  idiag = F, subject = 'id', ng = 1, data = pbmas_f) 

ltm_pbmas_f_k2k5 <- map(2:5, ~ { eval(parse(text = paste0(
    "gridsearch(hlme(fixed = tblh_bmc_ ~ 1 + ns(age, 7), random = ~ 1 + age, idiag = F, 
    subject = 'id', mixture = ~ ns(age, 7), ng = ", .x, ", nwg = T, data = pbmas_f), 
    minit = ltm_pbmas_f_k1, maxiter = 10, rep = 100)")))})

(ltm_pbmas_f_results <- c(list(ltm_pbmas_f_k1),ltm_pbmas_f_k2k5))

map_df(ltm_pbmas_f_results, ~ {return(data.frame(summarytable(.x, which = c(
  "G", "npm", "BIC", "entropy", "%class"))))})  

# PLOT TRAJECTORIES FROM 1-5 CLASS MODELS

graphics.off()
pdf("results/ltm_pbmas_f_k1k5.pdf",
    width = 7, height = 8)
par(
  mfrow    = c(3,2),
  mar      = c(5, 5, 2, 2)
)

iwalk(map(ltm_pbmas_f_results, ~ predictY(., data.frame(age = seq(min(pbmas_f$age), max(pbmas_f$age), length = 100)), var.time ="age")), ~{
  plot(.x, main = "", legend.loc = "topleft", ylim=c(500,4000), xlim=c(5,40),xlab="Age - years", ylab="BMC - grams", lty=1)
  mtext("BMC Latent Trajectories from 1 to 5 Class Models in PBMAS Females", outer = T, line = -1.5)})

dev.off()

# PLOT TRAJECTORIES FROM SELECTED MODEL & PPROB

ltm_pbmas_f_best <- ltm_pbmas_f_results[[4]]
ltm_pbmas_f_pred <- data.frame(age = seq(min(pbmas_f$age), max(pbmas_f$age), length = 50))
ltm_pbmas_f_pred <- predictY(ltm_pbmas_f_best, ltm_pbmas_f_pred, var.time ="age", draws = T)
ltm_pbmas_f_pred <- cbind(ltm_pbmas_f_pred$times, ltm_pbmas_f_pred$pred)
summarytable(ltm_pbmas_f_best, which = "%class")
postprob(ltm_pbmas_f_best)

ltm_pbmas_f_plot <- ggplot(
  data = ltm_pbmas_f_pred, aes(x = age)) +
  geom_line(aes(y = Ypred_class1, col = "Class 1 (39.4%)")) +
  geom_line(aes(y = Ypred_class2, col = "Class 2 (15.0%)")) +
  geom_line(aes(y = Ypred_class3, col = "Class 3 (32.3%)")) +
  geom_line(aes(y = Ypred_class4, col = "Class 4 (13.4%)")) + theme_classic() +
  geom_ribbon(aes(ymin = lower.Ypred_class1, ymax = upper.Ypred_class1), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class2, ymax = upper.Ypred_class2), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class3, ymax = upper.Ypred_class3), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class4, ymax = upper.Ypred_class4), alpha=0.1) +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  scale_y_continuous(breaks=seq(500, 3500, 500), limits=c(500, 3500)) + 
  labs(x = 'Age - years', y = 'BMC - grams') + scale_color_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.20, 0.75),
        legend.background = element_blank(),
        legend.box.background = element_rect(
          fill = "white", colour = "black")) +
  ggtitle("(e) BMC Latent Trajectories: PBMAS Females")

#### PBMAS MALES ####

# 1-5 CLASS MODELS & FIT STATISTICS 

ltm_pbmas_m_k1 <-  hlme(
  fixed = tblh_bmc_ ~ 1 + ns(age, 7), random = ~ 1 + age,
  idiag = F, subject = 'id', ng = 1, data = pbmas_m)

ltm_pbmas_m_k2k5 <- map(2:5, ~ { eval(parse(text = paste0(
    "gridsearch(hlme(fixed = tblh_bmc_ ~ 1 + ns(age, 7), random = ~ 1 + age, idiag = F, 
    subject = 'id', mixture = ~ ns(age, 7), ng = ", .x, ", nwg = T, data = pbmas_m), 
    minit = ltm_pbmas_m_k1, maxiter = 10, rep = 100)")))})

(ltm_pbmas_m_results <- c(list(ltm_pbmas_m_k1),ltm_pbmas_m_k2k5))

map_df(ltm_pbmas_m_results, ~ {return(data.frame(summarytable(.x, which = c(
  "G", "npm", "BIC", "entropy", "%class"))))})  

# PLOT TRAJECTORIES FROM 1-5 CLASS MODELS

graphics.off()
pdf("results/ltm_pbmas_m_k1k5.pdf",
    width = 7, height = 8)
par(
  mfrow    = c(3,2),
  mar      = c(5, 5, 2, 2)
)

iwalk(map(ltm_pbmas_m_results, ~ predictY(., data.frame(age = seq(min(pbmas_m$age), max(pbmas_m$age), length = 100)), var.time ="age")), ~{
  plot(.x, main = "", legend.loc = "topleft", ylim=c(500,4000), xlim=c(5,40),xlab="Age - years", ylab="BMC - grams", lty=1)
  mtext("BMC Latent Trajectories from 1 to 5 Class Models in PBMAS Males", outer = T, line = -1.5)})

dev.off()

# PLOT TRAJECTORIES FROM SELECTED MODEL & PPROB

ltm_pbmas_m_best <- ltm_pbmas_m_results[[3]]
ltm_pbmas_m_pred <- data.frame(age = seq(min(pbmas_m$age), max(pbmas_m$age), length = 50))
ltm_pbmas_m_pred <- predictY(ltm_pbmas_m_best, ltm_pbmas_m_pred, var.time ="age", draws = T)
ltm_pbmas_m_pred <- cbind(ltm_pbmas_m_pred$times, ltm_pbmas_m_pred$pred)
summarytable(ltm_pbmas_m_best, which = "%class")
postprob(ltm_pbmas_m_best)

ltm_pbmas_m_plot <- ggplot(
  data = ltm_pbmas_m_pred, aes(x = age)) +
  geom_line(aes(y = Ypred_class1, col = "Class 1 (37.5%)")) +
  geom_line(aes(y = Ypred_class2, col = "Class 2 (35.7%)")) +
  geom_line(aes(y = Ypred_class3, col = "Class 3 (26.8%)")) + theme_classic() +
  geom_ribbon(aes(ymin = lower.Ypred_class1, ymax = upper.Ypred_class1), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class2, ymax = upper.Ypred_class2), alpha=0.1) +
  geom_ribbon(aes(ymin = lower.Ypred_class3, ymax = upper.Ypred_class3), alpha=0.1) +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), limits=c(4, 40)) + 
  scale_y_continuous(breaks=seq(500, 3500, 500), limits=c(500, 3500)) + 
  labs(x = 'Age - years', y = 'BMC - grams') + scale_color_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.20, 0.80),
        legend.background = element_blank(),
        legend.box.background = element_rect(
          fill = "white", colour = "black")) +
  ggtitle("(f) BMC Latent Trajectories: PBMAS Males")

##### FIGURE 7 ####

graphics.off()
pdf("results/Fig7_gmm.pdf",
    width = 9, height = 10)
((ltm_alsp_f_plot | ltm_alsp_m_plot) / 
    (ltm_bmdcs_f_plot | ltm_bmdcs_m_plot) / 
    (ltm_pbmas_f_plot | ltm_pbmas_m_plot)) 
dev.off()
