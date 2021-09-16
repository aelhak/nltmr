#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al
#
# DESCRIPTIVE ANALYSES
#

library(tidyverse)
library(patchwork)

#### LOAD MAIN DATASET #####

bmc_all <- read.table("data/bmc_all_long.csv", header=TRUE, sep=",", na.strings=c("","NA"))
bmc_all$sex <- factor(bmc_all$sex, levels = c("Females", "Males"))

#### ANALYSIS DATASETS ####

bmc_m <- bmc_all %>% filter(sex == "Males") %>% drop_na()
bmc_f <- bmc_all %>% filter(sex == "Females") %>% drop_na()
alsp_m <- bmc_all %>% filter(cohort == "ALSPAC" & sex == "Males") %>% drop_na()
alsp_f <- bmc_all %>% filter(cohort == "ALSPAC" & sex == "Females") %>% drop_na() 
bmdcs_m <- bmc_all %>% filter(cohort == "BMDCS" & sex == "Males") %>% drop_na()
bmdcs_f <- bmc_all %>% filter(cohort == "BMDCS" & sex == "Females") %>% drop_na()
pbmas_m <- bmc_all %>% filter(cohort == "PBMAS" & sex == "Males") %>% drop_na()
pbmas_f <- bmc_all %>% filter(cohort == "PBMAS" & sex == "Females") %>% drop_na() 

#### FIGURE 2: SCATTERPLOT & LINEPLOTS OF OBSERVED BMC ####

bmc_scatterplot <- ggplot(
  data = subset(bmc_all, !is.na(tblh_bmc_)), aes(x = age, y = tblh_bmc_)) + theme_classic() +
  geom_point(aes(colour = cohort), size = 0.2) + facet_wrap(sex ~ ., strip.position="top") + 
  theme(legend.title = element_blank(), legend.position = "top", panel.grid.minor.x = element_blank()) +
  scale_y_continuous(breaks=seq(500,4500,1000)) + ylab('BMC - grams') + xlab("Age - years") +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggtitle("(a) Observed BMC Values")

bmc_lineplot <- ggplot(data = subset(bmc_all, !is.na(tblh_bmc_)), aes(x = age, y = tblh_bmc_)) + 
  geom_line(aes(colour = factor(id))) + scale_colour_identity() + facet_grid(cohort ~ sex) +
  scale_y_continuous(breaks=seq(500,4500,1000)) + ylab('BMC - grams') + xlab("Age - years") +
  scale_x_continuous(breaks=c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40)) + theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none", 
        panel.grid.minor.x = element_blank()) +
  ggtitle("(b) Observed BMC Individual Trajectories")  

graphics.off()
pdf("results/Fig2_descr.pdf",
    width = 5, height = 6.3)

(bmc_scatterplot / bmc_lineplot) + plot_layout(heights = c(1, 2.3))

dev.off() 

rm(bmc_scatterplot, bmc_lineplot)

#### eTABLE: N, MEAN, SD OF AGE AND BMC  BY VISIT, SEX AND COHORT ####

bmc_all %>% 
  group_by(visit, cohort, sex) %>%
  select(tblh_bmc_, age) %>% summarise_all(list(
    mean = ~round(mean(.,na.rm=T), 1), sd = ~round(sd(.,na.rm=T), 1), 
    min = ~round(min(.,na.rm=T), 1), max = ~round(max(.,na.rm=T), 1), 
    N = ~sum(!is.na(.)))) %>% drop_na() %>% rename(N = tblh_bmc__N) %>% select(
      cohort, sex, visit, N, age_mean, age_sd, age_min, age_max, tblh_bmc__mean, 
      tblh_bmc__sd,tblh_bmc__min, tblh_bmc__max) %>% arrange(cohort, sex) %>% 
  write_csv("results/descr_bmc.csv")

#### MEDIAN AND IQR NUMBER OF REPEAT BMC DATA ####

map(c("ALSPAC", "BMDCS", "PBMAS"), ~ {
  bmc_all %>% drop_na %>% filter(cohort == .x) %>%
    group_by(id) %>% summarise(n = n()) %>% summarise(across(.cols = n, .fns = list(
      Median = median, IQR = IQR), na.rm = TRUE, .names = "{col}_{fn}"))
})  
