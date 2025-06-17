################################################################################
#
#        STATISTICAL METHODS FOR BIOINFORMATICS - PART I (PROF. ALONSO)
#                     PRACTICAL EXAM STRATEGY SCRIPT
#
################################################################################

# This script is a playbook for the two most likely practical exam scenarios
# in Part I:
#   SCENARIO A: Longitudinal Data Analysis using a Multilevel Model
#   SCENARIO B: Missing Data Analysis using MI or IPW
#
# The exam will likely focus on one of these. Be prepared for both.

################################################################################
# S0: SETUP - LOAD LIBRARIES
################################################################################

# Load all the packages you might need at the beginning.
# The exam environment should have these pre-installed, but you should be
# prepared to install them yourself with install.packages("package_name").

# Install packages if needed (uncomment the next line):
# install.packages(c("lme4", "lmerTest", "mice", "VIM"))

# For Longitudinal Models
library(lme4)       # The primary package for linear mixed-effects models (lmer)
library(lmerTest)   # Provides p-values for lmer models

# For Missing Data
library(mice)       # The primary package for Multiple Imputation
library(VIM)        # For visualizing missing data patterns (e.g., aggr plot)

################################################################################
#
# SCENARIO A: LONGITUDINAL DATA ANALYSIS WITH MULTILEVEL MODELS
#
# This is the most common and critical scenario for Part I. The 2024 exam was
# entirely based on this.
#
################################################################################

# --- Step 1: Understand the Research Question & Data Structure ---

# The question will involve tracking an outcome over time for multiple subjects.
# Load the data (e.g., from an .xlsx or .csv file).
# library(readxl)
# data <- read_excel("Data Alonso Practical Part 2024.xlsx")

# Inspect the data
# head(data)
# str(data)

# Identify the key variables:
# - The response variable (e.g., 'y' - the SAT score)
# - The subject identifier (e.g., 'id')
# - The time variable (e.g., 'time')
# - The between-subject predictor (e.g., 'prog' - the program type)

# Ensure the data is in "long" format (multiple rows per subject), which is
# what lmer() requires.

# --- Step 2: Formulate the Hierarchical Model (CRITICAL STEP) ---

# The exam question will likely provide the mathematical model specification.
# Write it down in your script as comments to guide your coding and interpretation.
# Example from the 2024 exam:

# Level 1 (Within-subject change):
# Y_ij = pi_0i + pi_1i * time_j + epsilon_ij

# Level 2 (Between-subject differences):
# pi_0i = gamma_00 + gamma_01 * PROG_i + b_0i  (Model for the intercept)
# pi_1i = gamma_10 + gamma_11 * PROG_i          (Model for the slope - NO random slope in this example)

# --- Step 3: Fit the Model in R using lmer() ---

# Translate the mathematical model into the lmer() formula.

# The fixed effects part (the gammas) corresponds to the main formula:
# response ~ time * predictor
# This expands to: 1 (for gamma_00) + time (for gamma_10) + predictor (for gamma_01) + time:predictor (for gamma_11)

# The random effects part (the b's) corresponds to the `( ... | id)` part.
# A random intercept is (1 | id).
# A random intercept and random slope for time is (1 + time | id) or (time | id).

# Example from 2024 Exam: The model had a random intercept (b_0i) but no random slope.
# So the random part is `(1 | id)`.

# set.seed(123) # For reproducibility
# lmm_fit <- lmer(y ~ time * prog + (1 | id), data = data, REML = FALSE)

# NOTE: If the model included a random slope for time (b_1i), the formula would be:
# lmm_fit <- lmer(y ~ time * prog + (1 + time | id), data = data, REML = FALSE)

# --- Step 4: Interpret the Output (Key Skill for MCQs) ---

# Use summary() to get the detailed output.
# model_summary <- summary(lmm_fit)
# print(model_summary)

# HOW TO READ THE OUTPUT:

# 4.1 -- Random Effects Section:
#       Groups   Name        Variance  Std.Dev.
#       id       (Intercept) 84.02     9.166    <-- This is sigma_0_sq (Var of b_0i)
#       Residual             60.31     7.766    <-- This is sigma_epsilon_sq (Var of epsilon_ij)
#
#       * To answer "compare variability between individuals vs. within individuals":
#         Compare the `(Intercept)` Variance (between-subject) to the `Residual`
#         Variance (within-subject). In the example above, 84.02 > 60.31.

# 4.2 -- Fixed Effects Section:
#                  Estimate Std. Error t value Pr(>|t|)
# (Intercept)      104.3007   1.7274   60.38   <2e-16   <-- This is gamma_00
# time             -16.2555   1.8860   -8.62   <2e-16   <-- This is gamma_10
# prog              -0.9646   2.3020   -0.42   0.6761   <-- This is gamma_01
# time:prog          6.3187   2.5133    2.51   0.0135   <-- This is gamma_11
#
#       * INTERPRETATION:
#         - gamma_00: Average score for the reference group (prog=0) at time=0.
#         - gamma_10: Average slope (rate of change) for the reference group.
#         - gamma_01: Average difference in initial score between prog=1 and prog=0.
#                     A non-significant p-value (like 0.6761) means no significant
#                     difference at the start (as expected with randomization).
#         - gamma_11: The key interaction term. This is the average *difference in slopes*
#                     between prog=1 and prog=0. A significant positive value (like 6.3187)
#                     means the intervention group's score decreased *less steeply*
#                     (or increased more steeply) than the control group's.

# 4.3 -- Calculating Predicted Values:
# Q: "What is the average SAT score for children in the standard program (prog=0)
#     one month after starting (time=1)?"
#
# Use the estimated coefficients:
# Y_hat = gamma_00 + gamma_01*prog + gamma_10*time + gamma_11*prog*time
#
# predicted_score <- 104.30 + (-0.9646 * 0) + (-16.2555 * 1) + (6.3187 * 0 * 1)
# print(predicted_score)


################################################################################
#
# SCENARIO B: MISSING DATA ANALYSIS
#
# This is a possible alternative practical task, based on older exams and the
# fact that it's the final topic of the lecture series.
#
################################################################################

# --- Step 1: Diagnose the Missingness ---

# Load the data (e.g., from the titanic.txt file)
# titanic_data <- read.csv("titanic.txt")

# Use the VIM package to visualize the pattern of missingness.
# aggr_plot <- aggr(titanic_data, col=c('navyblue','red'), numbers=TRUE,
#                   sortVars=TRUE, labels=names(titanic_data), cex.axis=.7,
#                   gap=3, ylab=c("Histogram of missing data","Pattern"))

# Think about the mechanism. Is it likely MCAR, MAR, or MNAR?
# For the Titanic data, missing 'age' is likely related to 'survived' status,
# which is observed. This suggests the data are MAR, justifying MI or IPW.

# --- Step 2A: Multiple Imputation (MI) using the `mice` package ---

# 1. IMPUTE: Create m complete datasets.
# m <- 100 # Number of imputations
# imp_data <- mice(titanic_data, m = m, seed = 123)

# 2. ANALYZE: Run your desired analysis on each imputed dataset.
#    For example, a logistic regression to predict survival.
# analysis_fit <- with(imp_data, glm(survived ~ pclass + sex + age, family = binomial))

# 3. POOL: Combine the results from the `m` analyses using Rubin's Rules.
# pooled_results <- pool(analysis_fit)

# Get the final, pooled summary of the model.
# summary(pooled_results)


# --- Step 2B: Inverse Probability Weighting (IPW) ---

# 1. MODEL THE MISSINGNESS: Create an indicator for being a "complete case".
# titanic_data$is_complete <- as.numeric(complete.cases(titanic_data))

# Fit a model to predict the probability of being a complete case.
# Use only variables that are fully observed as predictors.
# prob_model <- glm(is_complete ~ survived + pclass, data = titanic_data, family = "binomial")

# 2. CALCULATE WEIGHTS: The weight is the inverse of the predicted probability.
# titanic_data$ipw_weights <- 1 / fitted(prob_model)

# 3. RUN WEIGHTED ANALYSIS: Run the analysis only on the complete cases,
#    using the calculated weights.
# ipw_fit <- glm(survived ~ pclass + sex + age,
#                data = titanic_data[titanic_data$is_complete == 1, ],
#                weights = ipw_weights[is_complete == 1],
#                family = "binomial")

# summary(ipw_fit)