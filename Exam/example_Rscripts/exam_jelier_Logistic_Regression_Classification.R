# ============================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE LOGISTIC REGRESSION (CLASSIFICATION)
# ============================================================================
# This script provides a workflow for a binary classification problem,
# focusing on the types of questions seen in Lecture 1 & 2 solutions.
# It covers fitting and comparing a standard Logistic Model, a non-linear
# Generalized Additive Model (GAM), and a Gradient Boosting Machine (GBM).

# --- 1. Load Necessary Libraries ---
# -----------------------------------
# 'gam' for Generalized Additive Models.
# 'gbm' for Gradient Boosting Machines.
# 'pROC' for creating and plotting ROC curves to compare models.
# You may need to install these first: install.packages(c("gam", "gbm", "pROC"))
library(gam)
library(gbm)
library(pROC)

# --- 2. Load and Prepare the Data ---
# ------------------------------------
# Load your data. For this example, we'll simulate data.
set.seed(1)
n <- 200
p <- 10
x <- matrix(rnorm(n * p), n, p)
# Create a true underlying model with non-linear and interaction effects
y_log_odds <- -1 + 2 * x[, 1] - 3 * x[, 2]^2 + 1.5 * x[, 3] * x[, 4]
y_prob <- exp(y_log_odds) / (1 + exp(y_log_odds))
y_binary <- rbinom(n, 1, y_prob)

# Create a data frame and ensure the response is a factor
data <- as.data.frame(x)
data$Treatment <- as.factor(ifelse(y_binary == 1, "Memantine", "Saline"))


# --- 3. Split Data into Training and Test Sets ---
# -------------------------------------------------
set.seed(1)
train_indices <- sample(1:n, 0.7 * n)
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]


# =========================================================================
# --- Question 1: Train and Compare Models ---
# =========================================================================

# --- Model A: Standard Logistic Regression ---
log_model <- glm(Treatment ~ ., data = train_data, family = "binomial")

# --- Model B: Generalized Additive Model (GAM) ---
# GAMs can capture non-linear relationships. We use smoothing splines `s()`
# on predictors we suspect might have a non-linear effect.
# The formula is built programmatically for all predictors.
gam_formula <- as.formula(paste("Treatment ~", paste(paste0("s(", names(train_data)[1:p], ", 4)"), collapse = " + ")))
gam_model <- gam(gam_formula, data = train_data, family = "binomial")

# Discussion Point: Is there evidence for non-linear effects?
# Yes, if the smooth terms in the GAM are significant. You can check this:
# summary(gam_model)
# And plot them to visualize the non-linear relationships.
# par(mfrow = c(2, 5)); plot(gam_model, se = TRUE, col = "blue")

# --- Model C: Gradient Boosting Machine (GBM) ---
# GBM is a powerful tree-based ensemble method that can capture
# complex interactions automatically.

# GBM requires numeric response (0/1)
train_data_gbm <- train_data
train_data_gbm$Treatment <- as.numeric(train_data_gbm$Treatment) - 1 # Saline=0, Memantine=1

# Note: Tuning GBM can be complex. For an exam, using reasonable defaults
# and explaining the key parameters is usually sufficient.
gbm_model <- gbm(
  Treatment ~ .,
  data = train_data_gbm,
  distribution = "bernoulli", # For logistic regression
  n.trees = 1000,              # Number of trees
  interaction.depth = 3,       # Allows for up to 3-way interactions
  shrinkage = 0.01,            # Learning rate
  cv.folds = 5                 # Use 5-fold CV to find the best number of trees
)

# Discussion Point: Is over-learning (overfitting) a problem?
# Yes, especially for flexible models like GBM. Using too many trees will
# overfit. Cross-validation is key to find the optimal number of trees that
# minimizes the CV error.
best_iter <- gbm.perf(gbm_model, method = "cv") # Find best number of trees from CV
# summary(gbm_model, n.trees = best_iter) # Shows variable importance

# Discussion Point: Is there evidence for important interactions?
# Yes, if the interaction.depth > 1 improves performance. GBM automatically
# models interactions between the variables deemed most important.


# =========================================================================
# --- Question 2: Evaluate Model Performance on the Test Set ---
# =========================================================================

# Use the models to predict probabilities on the test set
pred_log <- predict(log_model, newdata = test_data, type = "response")
pred_gam <- predict(gam_model, newdata = test_data, type = "response")
pred_gbm <- predict(gbm_model, newdata = test_data, n.trees = best_iter, type = "response")

# Calculate ROC curves and AUC values
roc_log <- roc(test_data$Treatment, pred_log, quiet = TRUE)
roc_gam <- roc(test_data$Treatment, pred_gam, quiet = TRUE)
roc_gbm <- roc(test_data$Treatment, pred_gbm, quiet = TRUE)

# Print AUC values (higher is better)
cat("AUC for Logistic Regression:", auc(roc_log), "\n")
cat("AUC for GAM:", auc(roc_gam), "\n")
cat("AUC for GBM:", auc(roc_gbm), "\n")

# Plot ROC curves to visually compare performance
plot(roc_log, col = "blue", main = "ROC Curve Comparison")
lines(roc_gam, col = "red")
lines(roc_gbm, col = "green")
legend("bottomright", legend = c("Logistic", "GAM", "GBM"), col = c("blue", "red", "green"), lwd = 2)

# Discussion: Based on AUC, which model performs best at separating the two
# treatment groups? The GBM is often the best performer due to its ability to
# model non-linearities and interactions, but this is data-dependent.


# =========================================================================
# --- Question 3: Evaluate Specific Changes (e.g., protein levels) ---
# =========================================================================
# Your task: Does Memantine treatment induce specific changes in protein
# levels of a certain group?

# This is often an "interaction" question. For example, you would fit a model
# that predicts a specific protein level (let's say V1) and include an
# interaction term between Treatment and another factor (e.g., Genotype).

# Hypothetical linear model:
# protein_model <- lm(V1 ~ Treatment * Genotype, data = your_subset_data)
# summary(protein_model)

# The answer lies in the significance of the 'Treatment:Genotype' interaction
# term. If this term is significant, it means the effect of the Treatment
# on protein V1 is DIFFERENT for different genotypes, allowing you to answer
# the question directly.

