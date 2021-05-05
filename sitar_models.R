#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al
#
# SITAR MODELS
#

library(sitar)
library(splines)
library(tidyverse)

#### ALSPAC FEMALES ####

sitar_alsp_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "sitar(x = age, y = tblh_bmc_, id = id, data = alsp_f, df = ", .x, ")")))
}, otherwise = NA_real_))

(sitar_alsp_f <- Filter(Negate(anyNA), sitar_alsp_f))
(sitar_alsp_f_bic <- unlist(map(sitar_alsp_f, BICadj)))
(sitar_alsp_f_best <- sitar_alsp_f[[which.min(sitar_alsp_f_bic)]])  # k = 2
rm(sitar_alsp_f, sitar_alsp_f_bic)
varexp(sitar_alsp_f_best)

(alsp_f_mean_pv = plot(sitar_alsp_f_best, 'v', apv = TRUE)$apv[2])
(alsp_f_mean_apv = plot(sitar_alsp_f_best, 'v', apv = TRUE)$apv[1])

#### ALSPAC MALES ####

sitar_alsp_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "sitar(x = age, y = tblh_bmc_, id = id, data = alsp_m, df = ", .x, ")")))
}, otherwise = NA_real_))

(sitar_alsp_m <- Filter(Negate(anyNA), sitar_alsp_m))
(sitar_alsp_m_bic <- unlist(map(sitar_alsp_m, BICadj)))
(sitar_alsp_m_best <- sitar_alsp_m[[which.min(sitar_alsp_m_bic)]])  # k = 5
rm(sitar_alsp_m, sitar_alsp_m_bic)
varexp(sitar_alsp_m_best)

(alsp_m_mean_pv = plot(sitar_alsp_m_best, 'v', apv = TRUE)$apv[2])
(alsp_m_mean_apv = plot(sitar_alsp_m_best, 'v', apv = TRUE)$apv[1])

#### BMDCS FEMALES ####

sitar_bmdcs_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "sitar(x = age, y = tblh_bmc_, id = id, data = bmdcs_f, df = ", .x, ")")))
}, otherwise = NA_real_))

(sitar_bmdcs_f <- Filter(Negate(anyNA), sitar_bmdcs_f))
(sitar_bmdcs_f_bic <- unlist(map(sitar_bmdcs_f, BICadj)))
(sitar_bmdcs_f_best <- sitar_bmdcs_f[[which.min(sitar_bmdcs_f_bic)]])  # k = 5
rm(sitar_bmdcs_f, sitar_bmdcs_f_bic)
varexp(sitar_bmdcs_f_best)

(bmdcs_f_mean_pv = plot(sitar_bmdcs_f_best, 'v', apv = TRUE)$apv[2])
(bmdcs_f_mean_apv = plot(sitar_bmdcs_f_best, 'v', apv = TRUE)$apv[1])

#### BMDCS MALES ####

sitar_bmdcs_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "sitar(x = age, y = tblh_bmc_, id = id, data = bmdcs_m, df = ", .x, ")")))
}, otherwise = NA_real_))

(sitar_bmdcs_m <- Filter(Negate(anyNA), sitar_bmdcs_m))
(sitar_bmdcs_m_bic <- unlist(map(sitar_bmdcs_m, BICadj)))
(sitar_bmdcs_m_best <- sitar_bmdcs_m[[which.min(sitar_bmdcs_m_bic)]])  # k = 4
rm(sitar_bmdcs_m, sitar_bmdcs_m_bic)
varexp(sitar_bmdcs_m_best)

(bmdcs_m_mean_pv = plot(sitar_bmdcs_m_best, 'v', apv = TRUE)$apv[2])
(bmdcs_m_mean_apv = plot(sitar_bmdcs_m_best, 'v', apv = TRUE)$apv[1])

#### PBMAS FEMALES ####

sitar_pbmas_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "sitar(x = age, y = tblh_bmc_, id = id, data = pbmas_f, df = ", .x, ")")))
}, otherwise = NA_real_))

(sitar_pbmas_f <- Filter(Negate(anyNA), sitar_pbmas_f))
(sitar_pbmas_f_bic <- unlist(map(sitar_pbmas_f, BICadj)))
(sitar_pbmas_f_best <- sitar_pbmas_f[[which.min(sitar_pbmas_f_bic)]])  # k = 2
rm(sitar_pbmas_f, sitar_pbmas_f_bic)
varexp(sitar_pbmas_f_best)

(pbmas_f_mean_pv = plot(sitar_pbmas_f_best, 'v', apv = TRUE)$apv[2])
(pbmas_f_mean_apv = plot(sitar_pbmas_f_best, 'v', apv = TRUE)$apv[1])

#### PBMAS MALES ####

sitar_pbmas_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
    "sitar(x = age, y = tblh_bmc_, id = id, data = pbmas_m, df = ", .x, ")")))
  }, otherwise = NA_real_))

(sitar_pbmas_m <- Filter(Negate(anyNA), sitar_pbmas_m))
(sitar_pbmas_m_bic <- unlist(map(sitar_pbmas_m, BICadj)))
(sitar_pbmas_m_best <- sitar_pbmas_m[[which.min(sitar_pbmas_m_bic)]])  # k = 2
rm(sitar_pbmas_m, sitar_pbmas_m_bic)
varexp(sitar_pbmas_m_best)

(pbmas_m_mean_pv = plot(sitar_pbmas_m_best, 'v', apv = TRUE)$apv[2])
(pbmas_m_mean_apv = plot(sitar_pbmas_m_best, 'v', apv = TRUE)$apv[1])

#### FIGURE 5: SITAR PLOTS FROM BEST MODELS ####

graphics.off()
pdf("results/Fig5_sitar.pdf",
    width=11.1, height=11.2)
par(
  mfrow    = c(3,2),
  mar      = c(5, 5, 3, 2)
)

plot(sitar_alsp_f_best, opt = c('d', 'v'), las = 1, apv = F, 
     legend = NULL, ylim = c(500,3500), xlim = c(5,40), vlim = c(0,400),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main="(a) SITAR BMC Trajectory and Velocity: ALSPAC Females", 
     font.main = 1, adj  = 0, cex.main=1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams", side = 2, line = 3.25)
mtext("BMC velocity - grams / year", side = 4, line = 2.6)
abline(v = alsp_f_mean_apv, lwd = 2, lty = 3, col = 'red')

plot(sitar_alsp_m_best, opt = c('d', 'v'), las = 1, apv = F,
     legend = NULL, ylim = c(500,3500), xlim = c(5,40), vlim = c(0,400),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main = "(b) SITAR BMC Trajectory and Velocity: ALSPAC Males", 
     font.main = 1, adj  = 0, cex.main = 1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams", side = 2, line = 3.25)
mtext("BMC velocity - grams / year", side = 4, line = 2.6)
abline(v = alsp_m_mean_apv, lwd = 2, lty = 3, col = 'red')

plot(sitar_bmdcs_f_best, opt = c('d', 'v'), las = 1, apv = F, 
     legend = NULL, ylim = c(500,3500), xlim = c(5,40), vlim = c(0,400),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main="(c) SITAR BMC Trajectory and Velocity: BMDCS Females", 
     font.main = 1, adj  = 0, cex.main=1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams", side = 2, line = 3.25)
mtext("BMC velocity - grams / year", side = 4, line = 2.6)
abline(v = bmdcs_f_mean_apv, lwd = 2, lty = 3, col = 'red')

plot(sitar_bmdcs_m_best, opt = c('d', 'v'), las = 1, apv = F, 
     legend = NULL, ylim = c(500,3500), xlim = c(5,40), vlim = c(0,400),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main="(d) SITAR BMC Trajectory and Velocity: BMDCS Males", 
     font.main = 1, adj  = 0, cex.main=1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams", side = 2, line = 3.25)
mtext("BMC velocity - grams / year", side = 4, line = 2.6)
abline(v = bmdcs_m_mean_apv, lwd = 2, lty = 3, col = 'red')

plot(sitar_pbmas_f_best, opt = c('d', 'v'), las = 1, apv = F, 
     legend = NULL, ylim = c(500,3500), xlim = c(5,40), vlim = c(0,400),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main="(e) SITAR BMC Trajectory and Velocity: PBMAS Females", 
     font.main = 1, adj  = 0, cex.main=1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams", side = 2, line = 3.25)
mtext("BMC velocity - grams / year", side = 4, line = 2.6)
abline(v = pbmas_f_mean_apv, lwd = 2, lty = 3, col = 'red')

plot(sitar_pbmas_m_best, opt = c('d', 'v'), las = 1, apv = F, 
     legend = NULL, ylim = c(500,3500), xlim = c(5,40), vlim = c(0,400),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main="(f) SITAR BMC Trajectory and Velocity: PBMAS Males", 
     font.main = 1, adj  = 0, cex.main=1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams", side = 2, line = 3.25)
mtext("BMC velocity - grams / year", side = 4, line = 2.6)
abline(v = pbmas_m_mean_apv, lwd = 2, lty = 3, col = 'red')

dev.off()