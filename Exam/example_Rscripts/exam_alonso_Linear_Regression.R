# =========================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE LINEAR REGRESSION PROBLEMS
# =========================================================================
# This script provides a complete workflow for analyzing a multiple linear
# regression model. It covers data exploration, model fitting, interpretation
# of both summary() and anova() outputs, model comparison, and prediction.
# We will use the "Patient Satisfaction" example from Lecture 1.

# --- 1. Load and Prepare Your Data ---
# -------------------------------------
# For this example, we assume you have a dataframe 'satisfaction' loaded.
# It contains 'satis', 'age', 'severity', and 'anxiety'.
#
# satisfaction <- read.csv("your_satisfaction_data.csv") # Example

# --- 2. Initial Data Exploration ---
# -----------------------------------
# Before modeling, always explore the relationships between your variables.
# The correlation matrix gives a first look at the strength and direction
# of linear relationships.
cor(satisfaction)

# A pairs plot (scatterplot matrix) is even better for visualizing these.
plot(satisfaction)
# From this, we'd see that 'age', 'severity', and 'anxiety' are all
# negatively correlated with patient 'satis'faction.

# =========================================================================
# STEP A: FIT AND INTERPRET THE FULL MODEL
# =========================================================================
# We fit a multiple linear regression model using the lm() function.
full_model <- lm(satis ~ age + severity + anxiety, data = satisfaction)

# The summary() function is your primary tool for interpretation.
summary(full_model)

# --- How to Interpret the summary() Output ---
# Coefficients:
#   - Estimate: The value for each β. For 'age', it's -1.142, meaning for
#     a one-year increase in age, satisfaction is predicted to decrease by
#     1.142 points, holding other variables constant.
#   - Pr(>|t|): The p-value for each coefficient. It tests if the coefficient
#     is significantly different from zero, GIVEN the other variables are in
#     the model. Here, 'severity' (p=0.374) is not significant.
#
# R-squared:
#   - Multiple R-squared: 0.682. This means the model explains 68.2% of the
#     variance in patient satisfaction.
#
# F-statistic:
#   - p-value: < 2.2e-16. This tests the overall model significance (H₀: β₁=β₂=β₃=0).
#     A tiny p-value means the model as a whole is statistically significant.

# =========================================================================
# STEP B: UNDERSTAND THE anova() FUNCTION (SEQUENTIAL TESTS)
# =========================================================================
# The anova() function provides a different perspective. It tests the
# variables sequentially, based on the order in the formula.
anova(full_model)

# --- How to Interpret the anova() Output ---
# - Row 'age': Tests if 'age' adds value compared to a model with only an intercept.
# - Row 'severity': Tests if 'severity' adds value to a model that ALREADY has 'age'.
# - Row 'anxiety': Tests if 'anxiety' adds value to a model that ALREADY has 'age' AND 'severity'.
# The p-values can change if you change the order of variables in the formula.

# =========================================================================
# STEP C: MODEL COMPARISON (REFINING THE MODEL)
# =========================================================================
# Since 'severity' was not significant in the full model summary, we can
# create a more parsimonious (simpler) model without it.
reduced_model <- lm(satis ~ age + anxiety, data = satisfaction)

# We can formally compare this 'nested' model to the full model using an LRT.
# A non-significant p-value means the simpler model is not significantly worse,
# so we prefer it based on parsimony.
anova(reduced_model, full_model)
# The result will show a high p-value, justifying the removal of 'severity'.

# Our final, chosen model is the reduced_model.
summary(reduced_model)

# =========================================================================
# STEP D: PREDICTION WITH THE FINAL MODEL
# =========================================================================
# Now we can use our final model to make predictions.
# Let's predict for a new patient who is 43 years old with an anxiety of 2.7.
newdata <- data.frame(age = 43, anxiety = 2.7)

# To get the 95% CONFIDENCE INTERVAL for the AVERAGE satisfaction of this group:
predict(reduced_model, newdata, interval = "confidence")
# This gives a narrower interval, e.g., (44, 54).

# To get the 95% PREDICTION INTERVAL for a SINGLE new patient's satisfaction:
predict(reduced_model, newdata, interval = "prediction")
# This gives a much wider interval, e.g., (28, 70), because it also accounts
# for individual variability (the error term ε).
