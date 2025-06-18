# =========================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE LINEAR REGRESSION PROBLEMS
# =========================================================================
# This script provides a complete workflow for tackling a linear regression
# problem on a high-dimensional dataset, as seen in Lecture 1 & 2.
# It covers data exploration, fitting a Lasso model, using cross-validation
# to tune the penalty parameter (lambda), and interpreting the results.

# --- 1. Load Necessary Libraries ---
# -----------------------------------
# 'glmnet' is essential for Ridge and Lasso regression.
# 'leaps' is useful for best subset and stepwise selection if needed.
# You may need to install these first: install.packages(c("glmnet", "leaps"))
library(glmnet)
library(leaps)

# --- 2. Load and Explore the Data ---
# ------------------------------------
# Load your data (e.g., from an .Rdata file or a CSV)
# load("Cortex2.rdata")
# For demonstration, we'll create simulated data similar to your exercises.
set.seed(1)
p <- 50   # Number of predictors (features)
n <- 100  # Number of observations (samples)
x <- matrix(rnorm(n * p), n, p)
B <- rnorm(p)
B[c(3, 10, 22, 31, 45)] <- 0 # Some coefficients are truly zero (sparse model)
y <- x %*% B + rnorm(n)

# Convert matrix to a data frame for easier handling
data <- as.data.frame(cbind(y, x))
colnames(data)[1] <- "response"


# --- Question 1: Data Exploration and Potential Issues ---
# ---------------------------------------------------------
# Your task: Study and describe the data. Identify potential issues.

# a) Check dimensions: Is the number of predictors (p) large relative to
#    the number of samples (n)? This is the "high-dimensionality" problem.
dim(data)
# If p >= n, standard OLS regression is not appropriate. This is a key issue
# to mention as it motivates the use of regularization methods like Lasso.

# b) Check for collinearity (correlations between predictors)
# High collinearity can make coefficient estimates unstable in standard
# regression. While Lasso handles this, it's an important issue to note.
cor_matrix <- cor(data[, -1]) # Exclude the response variable
# Look for high absolute values in the matrix. A heatmap is a great visual.
# For example, find the two most correlated predictors:
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
which(abs(cor_matrix) == max(abs(cor_matrix), na.rm = TRUE), arr.ind = TRUE)

# c) Check the response variable distribution
hist(data$response, main = "Distribution of the Response Variable")
# Is it roughly normal? Significant skew could be an issue.


# =========================================================================
# --- Question 2: Train a Lasso Model and Interpret ---
# =========================================================================
# Your task: Build a Lasso model, compare coefficients, and discuss.

# a) Prepare data for glmnet
# glmnet requires a matrix of predictors (x) and a vector for the response (y).
x_matrix <- as.matrix(data[, -1])
y_vector <- data$response

# b) Find the optimal lambda (penalty parameter) using cross-validation
# This is the most important step for tuning the model.
# CV helps prevent over-learning (overfitting) by selecting a lambda that
# performs well on unseen data, not just the training data.
set.seed(1) # for reproducibility
cv_lasso <- cv.glmnet(x_matrix, y_vector, alpha = 1) # alpha=1 for Lasso

# Plot the cross-validation error to visualize the tuning
plot(cv_lasso)
title("Cross-Validation Error vs. Log(Lambda) for Lasso", line = 3)
# The plot shows the MSE for different lambda values. The dotted lines indicate
# the lambda that gives the minimum MSE (`lambda.min`) and the lambda that is
# within one standard error of the minimum (`lambda.1se`). `lambda.1se`
# gives a more parsimonious (simpler) model and is often preferred.

# Get the optimal lambda values
lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se

# c) Fit the final Lasso model and extract coefficients
# We can look at the coefficients for both optimal lambdas.
coef_min <- coef(cv_lasso, s = lambda_min)
coef_1se <- coef(cv_lasso, s = lambda_1se)

print("Coefficients for lambda.min:")
print(coef_min)
print("Coefficients for lambda.1se (more parsimonious):")
print(coef_1se)

# Discussion Point: Can a subset of variables have good predictive power?
# Yes. Lasso performs variable selection by shrinking some coefficients
# to exactly zero. The non-zero coefficients in `coef_1se` represent the
# subset of variables selected by the model as most important.


# --- Question 3: Compare Lasso Coefficient to Correlation ---
# ------------------------------------------------------------
# Your task: Compare a specific coefficient with its correlation to the response.

# Let's pick a variable that was shrunk to zero, e.g., V3 (since B[3] was 0).
lasso_coef_V3 <- coef_1se["V3", 1]
correlation_V3 <- cor(data$response, data$V3)

print(paste("Lasso coefficient for V3:", round(lasso_coef_V3, 4)))
print(paste("Correlation of V3 with response:", round(correlation_V3, 4)))

# Discussion Point: Does this make sense?
# Yes, it can make sense. A variable might have a non-zero simple correlation
# with the response. However, in a multivariate context, its predictive
# contribution might be redundant if other, more correlated predictors are
# already in the model. Lasso considers all predictors simultaneously and will
# shrink the coefficient of a redundant predictor, even to zero, if it
# doesn't add unique explanatory power. This is a key strength of Lasso.

# =========================================================================
# --- Bonus: Compare Lasso to Best Subset Selection ---
# =========================================================================
# Best subset can be used if p is not too large. It helps confirm findings.
# We'll use the 'regsubsets' function from the 'leaps' package.

# Note: regsubsets has a limit on the number of variables (nvmax).
# This is computationally intensive and not suitable for very large p.
# if (p <= 50) {
#   best_subset_model <- regsubsets(response ~ ., data = data, nvmax = p)
#   bs_summary <- summary(best_subset_model)
#
#   # Plot BIC, Cp, and Adjusted R2 to find the best model size
#   par(mfrow = c(1, 3))
#   plot(bs_summary$bic, type = 'b', xlab = "Number of Variables", ylab = "BIC")
#   plot(bs_summary$cp, type = 'b', xlab = "Number of Variables", ylab = "Cp")
#   plot(bs_summary$adjr2, type = 'b', xlab = "Number of Variables", ylab = "Adjusted R2")
#
#   # The model with the lowest BIC/Cp or highest AdjR2 is often chosen.
#   # This provides an alternative set of "important" variables to compare with Lasso's.
# }

