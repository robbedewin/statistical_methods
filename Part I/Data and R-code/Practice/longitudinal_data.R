library("lme4")

early.int1 <- read.table("C:\\Users\\robbe\\Documents\\KULeuven\\Bioinformatics_2024-2025\\Statistical_Methods\\Part I\\Data and R-code\\Lecture3&4_earlyint.txt", header=T, sep=",")

attach(early.int1)
lmer_int1 <- lmer(cog ~ age * program + (1 | id), REML = FALSE)

summary(lmer_int1)



## Mean:
early.mean=tapply(cog,list(age,program),mean)

early.sd=tapply(cog,list(age,program),sd)

   ## Variance:
   early.var=tapply(cog,list(age,program),var)

   ## Frequency:
   early.n=table(age,program)
early.mean



library(readr)
long_teaching <- read_csv("~/KULeuven/Bioinformatics_2024-2025/Statistical_Methods/Exam/long_teaching.txt")
View(long_teaching)
attach(long_teaching)

lmer1 <- lmer(y ~ time * prog + (1 |id), data = long_teaching, REML = FALSE)

summary(lmer1)
