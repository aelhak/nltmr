#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al.
#
# MULTICOHORT SITAR MODELS
#

library(sitar)
library(splines)
library(tidyverse)

#### FEMALES ####

alsp_f2 <- alsp_f %>% mutate(bmcz = scale(tblh_bmc_))  
bmdcs_f2 <- bmdcs_f %>% mutate(bmcz = scale(tblh_bmc_))  
pbmas_f2 <- pbmas_f %>% mutate(bmcz = scale(tblh_bmc_))  
bmc_f2 <- rbind(alsp_f2, bmdcs_f2, pbmas_f2)
bmc_f2 <- bmc_f2 %>% filter(age >= 8.833333 & age <= 22.1)
bmc_f2$cohort <- as.factor(bmc_f2$cohort)

sitar_f <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "sitar(x = age, y = bmcz, id = id, data = bmc_f2, 
  a.form = ~cohort, b.form = ~cohort, c.form = ~cohort, df = ", .x, ")")))
}, otherwise = NA_real_))

(sitar_f <- Filter(Negate(anyNA), sitar_f))
(sitar_f_bic <- unlist(map(sitar_f, BICadj)))
(sitar_f_best <- sitar_f[[which.min(sitar_f_bic)]])  # k = 3
(sitar_f_mean_apv = plot(sitar_f_best, 'v', apv = TRUE)$apv[1])

#### MALES ####

alsp_m2 <- alsp_m %>% mutate(bmcz = scale(tblh_bmc_))  
bmdcs_m2 <- bmdcs_m %>% mutate(bmcz = scale(tblh_bmc_))  
pbmas_m2 <- pbmas_m %>% mutate(bmcz = scale(tblh_bmc_))  
bmc_m2 <- rbind(alsp_m2, bmdcs_m2, pbmas_m2)
bmc_m2 <- bmc_m2 %>% filter(age >= 8.75 & age <= 23.1) 
bmc_m2$cohort <- as.factor(bmc_m2$cohort)

sitar_m <- map(3:7, possibly(~ { eval(parse(text = paste0(
  "sitar(x = age, y = bmcz, id = id, data = bmc_m2, 
  a.form = ~cohort, b.form = ~cohort, c.form = ~cohort, df = ", .x, ")")))
}, otherwise = NA_real_))

(sitar_m <- Filter(Negate(anyNA), sitar_m))
(sitar_m_bic <- unlist(map(sitar_m, BICadj)))
(sitar_m_best <- sitar_m[[which.min(sitar_m_bic)]])  # k = 5
(sitar_m_mean_apv = plot(sitar_m_best, 'v', apv = TRUE)$apv[1])

#### FIGURE 9: SITAR PLOTS FROM BEST MODELS ####

plot(sitar_f_best, apv = T)
plot(sitar_m_best, apv = T)

graphics.off()
pdf("results/Fig9_mc_sitar.pdf",
    width = 7.4, height = 8.4)
par(
  mfrow = c(3,1),
  mar = c(5, 5, 3, 2)
  )

plot(sitar_f_best, opt = c('d', 'v'), las = 1, apv = F,
     legend = NULL, ylim = c(-1.5,2), xlim = c(5,25), vlim = c(0,0.6),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main = "(a) SITAR BMC Trajectory and Velocity: Females", 
     font.main = 1, adj  = 0, cex.main = 1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams in SD units", side = 2, line = 3.25)
mtext("Velocity - grams in SD units / year", side = 4, line = 2.6)
abline(v = sitar_f_mean_apv, lwd = 2, lty = 3, col = 'red')
text(13.2, -1.4, "12.2y", col = 'red', cex = 1.2)

plot(sitar_m_best, opt = c('d', 'v'), las = 1, apv = F, 
     legend = NULL, ylim = c(-1.5,2),  xlim = c(5,25), vlim = c(0,0.6),   
     vlab = "", ylab = "", xlab = "", lwd = 2, y2par = list(col = 'blue', lwd = 2),
     main = "(b) SITAR BMC Trajectory and Velocity: Males", 
     font.main = 1, adj  = 0, cex.main = 1.5)
mtext("Age - years", side = 1, line = 2)
mtext("BMC - grams in SD units", side = 2, line = 3.25)
mtext("Velocity - grams in SD units / year", side = 4, line = 2.6)
abline(v = sitar_m_mean_apv, lwd = 2, lty = 3, col = 'red')
text(14.89, -1.4, "13.9y", col = 'red', cex = 1.2)

dev.off()

rm(sitar_f, sitar_f_bic, alsp_f2, bmdcs_f2, pbmas_f2, bmc_f2, 
   sitar_m, sitar_m_bic, alsp_m2, bmdcs_m2, pbmas_m2, bmc_m2)
