library(lme4)
library(lmerTest)
library(nlme)

library(readr)
rat_pup <- read_table("~/KULeuven/Bioinformatics_2024-2025/Statistical_Methods/Part I/Data and R-code/Lecture3&4_rat_pup.txt", 
                                 col_types = cols(pupid = col_integer(), 
                                                  litterid = col_integer(), litsize = col_integer()))
View(rat_pup)

rat_pup$treatment <- as.factor(rat_pup$treatment)
rat_pup$sex <- as.factor(rat_pup$sex)

levels(rat_pup$treatment)

model_pup <- lmer(weight ~ sex * treatment * litsize + (1 | litterid), data = rat_pup, REML = FALSE)

summary(model_pup)
