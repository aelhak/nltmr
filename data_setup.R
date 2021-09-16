#
# Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise 
# nonlinear longitudinal growth trajectories in cohort studies
#
# Elhakeem et al: https://www.medrxiv.org/content/10.1101/2021.05.26.21257519v1
#
# DATA SETUP
#

library(tidyverse)
library(broom)

##############
### ALSPAC ###
##############

# LOAD DATA ----------------------
alsp_bmc_long <- read.table(
  "data/ALSPAC_data/ALSP_BMC_long.csv", 
                header=TRUE, sep=",", na.strings=c("","NA")
  )

# ORDER SEX -----
alsp_bmc_long$sex = factor(alsp_bmc_long$sex, levels=c("male","female"))
levels(alsp_bmc_long$sex) <- c("Males", "Females")

# CREATE VISIT ID  ---------
alsp_bmc_long <- alsp_bmc_long %>% mutate(
  cohort="ALSPAC",
  visit = ifelse(t==10, 1, ifelse(
    t==12, 2, ifelse(
      t==14, 3, ifelse(
        t==16, 4, ifelse(
          t==18, 5, 6))))))

# ORDER VISIT ----
alsp_bmc_long$visit = factor(alsp_bmc_long$visit, levels=c("1","2","3","4","5","6"))

#KEEP HARMONISED VARS ----
alsp_bmc_long <- alsp_bmc_long %>% select(-aln, -qlet, -t)

# SAVE HARMONISED DATASET AS CSV ----
write.csv(alsp_bmc_long,"data/alsp_bmc_long.csv", row.names = FALSE)

#############
### BMDCS ###
#############

# LOAD DATA ----------------------
bmdcs_bmc_long <- read.table(
  "data/BMDCS_data/BMDCS_BMC_long.csv", 
  header=TRUE, sep=",", na.strings=c("","NA")
)

# ORDER SEX -----
bmdcs_bmc_long <- bmdcs_bmc_long %>% mutate(
  sex=recode(sex,`1`="Males",`2`="Females"))

bmdcs_bmc_long$sex = factor(bmdcs_bmc_long$sex, levels=c("Males","Females"))
levels(bmdcs_bmc_long$sex) <- c("Males", "Females")

# VISIT ID VARIABLE ---------

bmdcs_bmc_long <- bmdcs_bmc_long %>% mutate(
  cohort="BMDCS",
  visit = ifelse(t==0, 1, ifelse(
    t==1, 2, ifelse(
      t==2, 3, ifelse(
        t==3, 4, ifelse(
          t==4, 5, ifelse(
            t==5, 6, 7)))))))

# ORDER VISIT ----
bmdcs_bmc_long$visit = factor(bmdcs_bmc_long$visit, levels=c("1","2","3","4","5","6","7"))

#KEEP HARMONISED VARS ----
bmdcs_bmc_long <- bmdcs_bmc_long %>% select(-randid, -t) %>% mutate(id = id+12075)

# SAVE HARMONISED DATASET AS CSV ----
write.csv(bmdcs_bmc_long,"data/bmdcs_bmc_long.csv", row.names = FALSE)

#############
### PBMAS ###
#############

# LOAD DATA ----------------------
pbmas_bmc_long <- read.table(
  "data/PBMAS_data/PBMAS_Trajectory_data_Nov_2019.csv", 
  header=TRUE, sep=",", na.strings=c("","NA")
)

# SUMMARY OF DATA ----
summary(pbmas_bmc_long)

# RESTRICT TO CUCASIAN ----
# RENAME VARIABLES ----
# REORDER SEX ----
pbmas_bmc_long <- pbmas_bmc_long %>% filter(
  Race!='Chinese' & Race!='East Indian' & Race!='Eritrean' & 
    Race!='Japanese(father)' & Race!='Native_Canadian') %>% rename(
  age = Age, id = ï..ID, sex = sex..1.male., tblh_bmc_ = TB.BMC...Head.BMC) %>% mutate(
    sex=recode(sex,`1`="Males",`0`="Females")) %>% select(
      id, sex, age, tblh_bmc_, Sequence)  %>% filter(Sequence!=32)

pbmas_bmc_long$sex = factor(pbmas_bmc_long$sex, levels=c("Males","Females"))
levels(pbmas_bmc_long$sex) <- c("Males", "Females")

# VISIT ID VARIABLE ---------
table(pbmas_bmc_long$Sequence)

pbmas_bmc_long <- pbmas_bmc_long %>% mutate(
  cohort="PBMAS",
  visit = ifelse(Sequence==2, 1, ifelse(
    Sequence==4, 2, ifelse(
      Sequence==6, 3, ifelse(
        Sequence==8, 4, ifelse(
          Sequence==10, 5, ifelse(
            Sequence==12, 6, ifelse(
              Sequence==14, 7, ifelse(
                Sequence==16, 8, ifelse(
                  Sequence==24, 9, ifelse(
                    Sequence==26, 10, ifelse(
                      Sequence==28, 11, ifelse(
                        Sequence==30, 12, ifelse(
                          Sequence==40, 13, ifelse(
                            Sequence==42, 14, ifelse(
                              Sequence==52, 15, 16))))))))))))))))
table(pbmas_bmc_long$visit)

# ORDER VISIT ----
pbmas_bmc_long$visit = factor(
  pbmas_bmc_long$visit, levels=c(
    "1", "2", "3", "4", "5", "6", "7", "8",
    "9", "10", "11", "12", "13", "14", "15", "16"))

pbmas_bmc_long <- pbmas_bmc_long %>% select(-Sequence) %>%
  group_by(id) %>% mutate(id = cur_group_id()+13028)

# SAVE HARMONISED DATASET AS CSV ----
write.csv(pbmas_bmc_long,"data/pbmas_bmc_long.csv", row.names = FALSE)

##################################
### COMBINED HARMONISED DATASET ##
##################################

bmc_all <- rbind(alsp_bmc_long, bmdcs_bmc_long, pbmas_bmc_long)

# ORDER SEX AND VISIT ----------

bmc_all$sex <- factor(
  bmc_all$sex, levels = c("Males", "Females")
)

bmc_all$visit = factor(
  bmc_all$visit, levels=c(
    "1", "2", "3", "4", "5", "6", "7", "8",
    "9", "10", "11", "12", "13", "14", "15", "16"))

bmc_all$sex = factor(bmc_all$sex, levels=c("Males","Females"))
levels(bmc_all$sex) <- c("Males", "Females")
bmc_all$visit <- as.factor(bmc_all$visit)
bmc_all$cohort <- as.factor(bmc_all$cohort)

# DUPLICATE VISIT DATA REMOVED [2 observations in PBMAS MALES]

bmc_all <- bmc_all %>% group_by(id, visit) %>% distinct()

# LIKELY ERRORS REMOVED [4 observations]

bmc_all %>%
  group_by(id) %>%
  arrange(id) %>%
  mutate(diff = tblh_bmc_ - lag(tblh_bmc_, default = first(tblh_bmc_))) %>% 
  filter(tblh_bmc_ <1000 & visit != 1 & age < 30 & diff <= -600)

bmc_all <- bmc_all[-c(2567, 32681, 59351, 79907), ]

summary(bmc_all)
write.csv(bmc_all,"data/bmc_all_long.csv", row.names = FALSE)
rm(alsp_bmc_long, bmdcs_bmc_long, pbmas_bmc_long)
