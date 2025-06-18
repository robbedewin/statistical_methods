# =========================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE HIGH-DIMENSIONAL REGRESSION
# =========================================================================
# This script covers the standard workflow for a high-dimensional dataset
# (where predictors p > samples n), a key topic from Lectures 3 & 4.
# It includes data exploration, PCR, Ridge, and Lasso regression.

# --- 1. Load Necessary Libraries ---
# -----------------------------------
# 'glmnet' for Ridge and Lasso regression.
# 'pls' for Principal Component Regression (PCR) and Partial Least Squares (PLS).
# You may need to install these first: install.packages(c("glmnet", "pls"))
library(glmnet)
library(pls)

# --- 2. Load and Prepare the Data ---
# ------------------------------------
# Load your data. For demonstration, we'll use the 'Hitters' dataset
# from your lab, as it has a moderate number of predictors (p=19).
# For a true p > n problem like 'Vijver', the principles are the same.
library(ISLR)
data <- na.omit(Hitters) # Remove rows with missing values

# --- Question 1: Data Exploration and Potential Issues ---
# ---------------------------------------------------------
# Your task: Explore the dataset. What challenges does it present?

# a) High-Dimensionality: Check if p is large relative to n.
dim(data)
# With many predictors, there's a high risk of overfitting and model instability.
# This motivates using dimension reduction or regularization.

# b) Collinearity: Check for correlations between predictors.
# Highly correlated predictors make standard regression coefficients unreliable.
cor_matrix <- cor(data[, sapply(data, is.numeric)]) # Correlation of numeric columns
# A heatmap is a good way to visualize this.
# High collinearity is a major challenge that Ridge and PCR are good at handling.

# =========================================================================
# --- Question 2: Fit and Compare Models ---
# =========================================================================
# Your task: Use Lasso, Ridge, and PCR. Evaluate performance and compare.

# a) Prepare data for modeling
# Create a matrix of predictors (x) and a vector of the response (y).
x <- model.matrix(Salary ~ ., data)[, -1] # Predictor matrix
y <- data$Salary                               # Response vector

# Create training and test sets for evaluation
set.seed(1)
train_indices <- sample(1:nrow(x), nrow(x) / 2)
test_indices <- (-train_indices)
y_test <- y[test_indices]


# --- Model A: Principal Component Regression (PCR) ---
# PCR reduces dimensionality by creating uncorrelated principal components.
set.seed(1)
pcr_fit <- pcr(Salary ~ ., data = data, subset = train_indices,
               scale = TRUE, validation = "CV")

# Use cross-validation plot to select the optimal number of components
validationplot(pcr_fit, val.type = "MSEP")
# The plot helps identify the number of components that minimizes test error.
# Let's assume the optimal number is 7 based on the plot.

# Calculate the test MSE for the optimal PCR model
pcr_pred <- predict(pcr_fit, x[test_indices, ], ncomp = 7)
pcr_mse <- mean((pcr_pred - y_test)^2)
print(paste("PCR Test MSE:", pcr_mse))


# --- Model B: Ridge Regression (alpha = 0) ---
# Ridge shrinks coefficients but rarely sets them to exactly zero.
# It's good for handling collinearity.
set.seed(1)
cv_ridge <- cv.glmnet(x[train_indices, ], y[train_indices], alpha = 0)
plot(cv_ridge) # Plot CV error vs. log(lambda)
bestlam_ridge <- cv_ridge$lambda.min

# Calculate the test MSE for the optimal Ridge model
ridge_pred <- predict(cv_ridge, s = bestlam_ridge, newx = x[test_indices, ])
ridge_mse <- mean((ridge_pred - y_test)^2)
print(paste("Ridge Test MSE:", ridge_mse))


# --- Model C: Lasso Regression (alpha = 1) ---
# Lasso performs variable selection by shrinking some coefficients to zero.
set.seed(1)
cv_lasso <- cv.glmnet(x[train_indices, ], y[train_indices], alpha = 1)
plot(cv_lasso)
bestlam_lasso <- cv_lasso$lambda.min

# Calculate the test MSE for the optimal Lasso model
lasso_pred <- predict(cv_lasso, s = bestlam_lasso, newx = x[test_indices, ])
lasso_mse <- mean((lasso_pred - y_test)^2)
print(paste("Lasso Test MSE:", lasso_mse))


# --- Model Interpretation and Comparison ---
# -------------------------------------------
# Your task: How many genes/predictors are used? Evaluate performance.

# 1. Performance: Compare the MSE values. Which model performs best?
#    The best model is the one with the lowest test MSE.

# 2. Variable Selection (Lasso):
#    Look at the coefficients from the final Lasso model.
final_lasso_coef <- predict(cv_lasso, type = "coefficients", s = bestlam_lasso)
print(final_lasso_coef)
# The non-zero coefficients are the variables selected by Lasso. This provides
# a sparse, interpretable model.

# 3. Discussion on Overfitting:
#    The use of cross-validation (CV) in all three methods is the key strategy
#    to combat overfitting. CV estimates the test error, which allows us to tune
#    the model's complexity (number of components for PCR, lambda for Ridge/Lasso)
#    to a level that generalizes well to new data, rather than just memorizing
#    the training data.

