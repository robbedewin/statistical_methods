# =========================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE LOGISTIC REGRESSION PROBLEMS
# =========================================================================
# This script provides a complete workflow for analyzing a logistic regression
# model. It covers model fitting, converting coefficients to odds ratios,
# performing model selection with AIC, and predicting probabilities.
# We will use the "Donner Party" example from Lecture 2.

# --- 1. Load Necessary Libraries ---
# -----------------------------------
# We need 'AICcmodavg' for easy model selection tables.
# You may need to install it first: install.packages("AICcmodavg")
library(AICcmodavg)

# --- 2. Load and Prepare Your Data ---
# -------------------------------------
# For this example, we assume you have a dataframe 'donner_data' loaded.
# It contains 'Outcome' (1=Survived, 0=Died), 'Age', and 'fem' (1=Female, 0=Male).
#
# donner_data <- read.csv("your_donner_data.csv") # Example

# =========================================================================
# STEP A: MODEL BUILDING (DEFINING CANDIDATE MODELS)
# =========================================================================
# In logistic regression, we often have several plausible theories (models)
# to explain our binary outcome. Based on the lecture, we define four
# candidate models for the Donner Party survival.

# We use glm() with family = "binomial" for logistic regression.
model1 <- glm(Outcome ~ Age, data = donner_data, family = "binomial")
model2 <- glm(Outcome ~ fem, data = donner_data, family = "binomial")
model3 <- glm(Outcome ~ Age + fem, data = donner_data, family = "binomial")
model4 <- glm(Outcome ~ Age * fem, data = donner_data, family = "binomial") # Includes interaction

# =========================================================================
# STEP B: MODEL SELECTION USING AIC
# =========================================================================
# We use AIC to find the "best" model among our candidates. AIC balances
# model fit with complexity. The model with the lowest AIC is preferred.

# Create a list of the candidate models.
candidate_models <- list(model1, model2, model3, model4)
model_names <- c("Age only", "Sex only", "Age + Sex", "Interaction")

# The aictab() function creates a summary table with AIC values and weights.
aic_table <- aictab(cand.set = candidate_models, modnames = model_names)
print(aic_table)

# --- How to Interpret the AIC Table ---
# - The table is ranked from best (lowest AICc) to worst model.
# - Delta_AICc: The difference in AICc from the best model.
# - AICcWt: The Akaike Weight. This is the probability that the model is the
#   best among the set.
#
# For the Donner data, 'Age + Sex' (model3) has the highest weight (~0.56),
# making it the most plausible model. The interaction model is second (~0.25).
# This provides a quantitative way to compare our theories.

# =========================================================================
# STEP C: INTERPRET THE CHOSEN MODEL (ODDS RATIOS)
# =========================================================================
# Let's focus on our best model, model3.
final_model <- model3
summary(final_model)

# The coefficients from summary() are in LOG-ODDS, which are hard to interpret.
# We must convert them to ODDS RATIOS by exponentiating.
# OR = exp(coefficient)

# Get the odds ratios for the coefficients.
exp(coef(final_model))

# Get the 95% confidence intervals for the odds ratios.
exp(confint(final_model))

# --- How to Interpret Odds Ratios ---
# - OR > 1: The predictor increases the odds of the outcome.
# - OR < 1: The predictor decreases the odds of the outcome.
# - OR = 1: The predictor has no effect on the odds of the outcome.
#
# For 'fem', the OR is ~3.0. Interpretation: The odds of survival for a
# woman are 3 times the odds of survival for a man of the same age.
# For a 10-year increase in 'Age', the OR is exp(coef*10) ≈ 0.70. Interpretation:
# A 10-year increase in age is associated with a 30% decrease in the odds of survival.

# =========================================================================
# STEP D: PREDICTING PROBABILITIES
# =========================================================================
# We can use our final model to predict the probability of survival for new individuals.

# Create a dataframe with the new data.
newdata <- data.frame(Age = 25, fem = 1) # A 25-year-old female

# Use predict() with type = "response" to get the probability directly.
predicted_prob <- predict(final_model, newdata, type = "response")
print(predicted_prob)
# The result will be a probability, e.g., ~0.70 or 70% chance of survival.
