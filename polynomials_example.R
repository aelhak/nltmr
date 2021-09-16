#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al: https://www.medrxiv.org/content/10.1101/2021.05.26.21257519v1

# POLYNOMIALS EXAMPLE

# The Berkeley Child Guidance Study dataset contains longitudinal 
# anthropometry data for 136 children from birth to 21 years.

library(lme4)
library(tidyverse)
library(RColorBrewer)

berkeley2 <- sitar::berkeley %>% filter(sex==2) %>% drop_na()

wt_line <- lmer(weight ~ age + (age | id), data = berkeley2)
wt_quad <- lmer(weight ~ poly(age, 2) + (age | id), data = berkeley2)
wt_cube <- lmer(weight ~ poly(age, 3) + (age | id), data = berkeley2)

(wt_preds <- data.frame(age = seq(min(berkeley2$age), max(berkeley2$age), length = 20)))
(wt_preds$predlm <- predict(wt_line, wt_preds, re.form = NA))
(wt_preds$predq2 <- predict(wt_quad, wt_preds, re.form = NA))
(wt_preds$predq3 <- predict(wt_cube, wt_preds, re.form = NA))

graphics.off()
pdf("results/Fig1_poly.pdf",
    width=6, height=4)
ggplot(
  data = berkeley2,
  mapping = aes(x = age, y = weight)) +
  geom_point(size = 0.1) +
  geom_line(data = wt_preds, aes(col = "(a) Linear Fit", y = predlm)) +
  geom_line(data = wt_preds, aes(col = "(b) Quadratic Poly (age\xB2)", y = predq2)) +
  geom_line(data = wt_preds, aes(col = "(c) Cubic Poly (age\xB3)", y = predq3)) +
  scale_color_brewer(palette = "Set1") + scale_x_continuous(breaks=c(8, 10, 12, 14, 16, 18, 20)) +
  labs(x = 'age in years', y = 'weight in kilograms') + theme_classic() +
  theme(legend.title = element_blank(), legend.position = "right", legend.direction = "vertical")
dev.off()

rm(berkeley2, wt_line, wt_quad, wt_cube, wt_preds)
