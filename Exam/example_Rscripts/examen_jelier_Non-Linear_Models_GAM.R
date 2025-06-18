# =========================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE NON-LINEAR MODELS (GAM)
# =========================================================================
# This script provides a workflow for fitting and evaluating Generalized
# Additive Models (GAMs), a key topic from Lecture 4. It covers identifying
# non-linear effects, comparing models with ANOVA, and evaluating performance.

# --- 1. Load Necessary Libraries ---
# -----------------------------------
# 'ISLR' for the 'College' dataset.
# 'gam' for fitting Generalized Additive Models.
library(ISLR)
library(gam)

# --- 2. Load and Prepare the Data ---
# ------------------------------------
data <- College

# Create training and test sets
set.seed(1)
train_indices <- sample(1:nrow(data), nrow(data) / 2)
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]


# =========================================================================
# --- Question 1: Fit a GAM and Evaluate Non-Linearity ---
# =========================================================================
# Your task: Fit a GAM. Are there evidence for non-linear effects?

# a) Fit the GAM
# We use the gam() function. For potentially non-linear predictors, we wrap
# them in `s()` for a smoothing spline. The `df` argument controls the flexibility.
gam_fit <- gam(Outstate ~ Private + s(Room.Board, 4) + s(PhD, 4) +
                 s(perc.alumni, 4) + s(Expend, 4) + s(Grad.Rate, 4),
               data = train_data)

# b) Visualize the effects
# Plotting the GAM fit shows the relationship for each predictor.
# If a curve is not a straight line, it indicates a non-linear effect.
# The dashed lines represent the confidence interval for the fit.
par(mfrow = c(2, 3))
plot(gam_fit, se = TRUE, col = "blue")
# Observation: The plot for `Expend` shows a clear non-linear relationship.
# Other variables like `Room.Board` appear mostly linear.

# c) Formal Test for Non-Linearity
# The summary of the GAM provides an "Anova for Parametric Effects" and an
# "Anova for Non-parametric Effects". The p-values (Pr(F)) in the second
# table test if the smooth term is significantly non-linear.
summary(gam_fit)
# For `s(Expend, 4)`, the p-value is very small, confirming that its
# non-linear effect is statistically significant.


# =========================================================================
# --- Question 2: Compare GAM to a Linear Model ---
# =========================================================================
# Your task: Can you build a better model with GAMs? Compare them.

# a) Fit a standard Linear Model (which is just a GAM with all linear terms)
lm_fit <- lm(Outstate ~ Private + Room.Board + PhD + perc.alumni + Expend + Grad.Rate,
             data = train_data)

# b) Compare models using ANOVA
# We can use an ANOVA test to formally check if the more complex GAM provides
# a significantly better fit than the simpler linear model.
# NOTE: To do this, you would create a nested GAM with only specific terms
# being smooth, like in your exercise file.
# For example, let's test if making `Expend` non-linear is justified.
gam_simplified <- gam(Outstate ~ Private + Room.Board + PhD + perc.alumni + s(Expend, 4) + Grad.Rate,
                      data = train_data)
anova(lm_fit, gam_simplified, test = "F")

# A significant p-value from this ANOVA would indicate that the GAM with the
# non-linear term for 'Expend' is significantly better than the purely linear model.

# c) Compare models based on Test Set Performance
# The ultimate test is how well the models predict on unseen data.
pred_lm <- predict(lm_fit, newdata = test_data)
pred_gam <- predict(gam_fit, newdata = test_data)

# Calculate Test Mean Squared Error (MSE)
mse_lm <- mean((pred_lm - test_data$Outstate)^2)
mse_gam <- mean((pred_gam - test_data$Outstate)^2)

print(paste("Linear Model Test MSE:", mse_lm))
print(paste("GAM Test MSE:", mse_gam))

# Discussion: If the GAM's test MSE is lower than the linear model's, it confirms
# that the flexible GAM has captured true non-linear patterns in the data,
# leading to better predictive performance. This demonstrates the value of
# moving "beyond linearity".

