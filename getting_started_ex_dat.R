#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al
#
# GETTING STARTED WITH THE EXAMPLE DATASET
#

library(tidyverse)
library(sitar)
library(splines)

#### LOAD DATASET FROM GITHUB ####

dat <- read.csv(("https://raw.githubusercontent.com/aelhak/nltmr/main/synth_cohort.csv"))

#### SCATTERPLOT OF DATASETS ####

ggplot(
  data = dat, aes(x = age, y = bmc)) + theme_bw() +
  geom_point(size = 0.2) + facet_wrap(sex ~ ., strip.position="top") + 
  theme(legend.title = element_blank(), legend.position = "top", panel.grid.minor.x = element_blank()) +
  scale_y_continuous(breaks=seq(500,4500,1000)) + ylab('BMC - grams') + xlab("Age - years") +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40)) +
  guides(colour = guide_legend(override.aes = list(size=5)))

#### CREATE MALE AND FEMALE DATASETS ####

dat_m <- dat %>% filter(sex == "Males")
dat_f <- dat %>% filter(sex == "Females")
