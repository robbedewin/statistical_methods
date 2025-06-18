# =========================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE TREE-BASED MODELS
# =========================================================================
# This script provides a complete workflow for classification problems using
# decision trees, bagging, Random Forests, and boosting, as covered in
# Lecture 5. It is structured to answer typical exam questions on model
# fitting, tuning, performance evaluation, and interpretation.

# --- 1. Load Necessary Libraries ---
# -----------------------------------
# 'tree' for fitting single classification and regression trees.
# 'randomForest' for bagging and Random Forest models.
# 'gbm' for Gradient Boosting Machines.
# You may need to install these first: install.packages(c("tree", "randomForest", "gbm"))
library(tree)
library(randomForest)
library(gbm)

# --- 2. Load and Prepare the Data ---
# ------------------------------------
# For this example, we'll use the 'SAheart.data' dataset from your exercise.
# We assume it is in your working directory.
heart <- read.csv('SAheart.data', header = TRUE, row.names = 1)

# For classification, the response variable 'chd' (coronary heart disease)
# must be a factor. We'll also save the original numeric version for gbm.
chd_original <- heart$chd
heart$chd <- as.factor(heart$chd)

# Create training and test sets
set.seed(1)
train_indices <- sample(1:nrow(heart), nrow(heart) * 2/3)
test_indices <- (-train_indices)
test_data <- heart[test_indices, ]


# =========================================================================
# --- Question 1: Fit and Prune a Single Classification Tree ---
# =========================================================================
# Your task: Build a single tree, prune it using cross-validation, and
# evaluate its performance.

# a) Fit a large, unpruned tree
# The `tree()` function fits a classification tree because the response 'chd' is a factor.
tree_full <- tree(chd ~ ., data = heart, subset = train_indices)

# Inspect and plot the tree
summary(tree_full)
plot(tree_full)
text(tree_full, pretty = 0, cex = 0.7)
# Interpretation: The summary shows the variables used for splitting and the
# misclassification error on the training set. The plot visualizes the decisions.

# b) Prune the tree using cross-validation
# `cv.tree` performs CV to see how the error rate changes with tree size.
# `FUN = prune.misclass` specifies that classification error should be the metric.
set.seed(1)
cv_tree <- cv.tree(tree_full, FUN = prune.misclass)
plot(cv_tree$size, cv_tree$dev, type = "b",
     xlab = "Tree Size", ylab = "CV Misclassification Error")
# The plot shows that the CV error is lowest for a certain tree size.
# Let's choose the best size, for instance, a tree with 6 terminal nodes.
best_size <- 6

# c) Create the pruned tree
tree_pruned <- prune.misclass(tree_full, best = best_size)
plot(tree_pruned)
text(tree_pruned, pretty = 0, cex = 0.8)

# d) Evaluate performance on the test set
# Predict the class for the test data using both trees.
pred_full_tree <- predict(tree_full, newdata = test_data, type = "class")
pred_pruned_tree <- predict(tree_pruned, newdata = test_data, type = "class")

# Calculate accuracy
accuracy_full <- mean(pred_full_tree == test_data$chd)
accuracy_pruned <- mean(pred_pruned_tree == test_data$chd)

print(paste("Accuracy of Full Tree:", accuracy_full))
print(paste("Accuracy of Pruned Tree:", accuracy_pruned))
# Discussion: Often, the pruned tree performs better on test data because
# it is less overfit to the training data. A good score is one that is
# significantly better than the null error rate (predicting the majority class).

# =========================================================================
# --- Question 2: Fit Bagging and Random Forest Models ---
# =========================================================================
# Your task: Compare Bagging and Random Forest. Interpret variable importance.

# a) Fit a Bagged Trees model
# Bagging is a special case of Random Forest where all predictors are considered
# at each split (mtry = number of predictors).
set.seed(1)
bag_fit <- randomForest(chd ~ ., data = heart, subset = train_indices,
                        mtry = ncol(heart) - 1, # Use all predictors
                        importance = TRUE, ntree = 1000)

# b) Fit a Random Forest model
# Here, mtry is typically set to sqrt(p) for classification.
set.seed(1)
rf_fit <- randomForest(chd ~ ., data = heart, subset = train_indices,
                       mtry = floor(sqrt(ncol(heart) - 1)), # Use a subset of predictors
                       importance = TRUE, ntree = 1000)

# c) Evaluate performance
pred_bag <- predict(bag_fit, newdata = test_data)
pred_rf <- predict(rf_fit, newdata = test_data)
accuracy_bag <- mean(pred_bag == test_data$chd)
accuracy_rf <- mean(pred_rf == test_data$chd)

print(paste("Bagging Accuracy:", accuracy_bag))
print(paste("Random Forest Accuracy:", accuracy_rf))
# Discussion: Random Forest often outperforms Bagging because de-correlating
# the trees by using a random subset of predictors at each split can reduce variance.

# d) Interpret Variable Importance
# The `importance()` function shows how much each variable contributes to
# reducing impurity (e.g., Gini index). `varImpPlot` visualizes this.
par(mfrow = c(1, 2))
varImpPlot(bag_fit, main = "Variable Importance: Bagging")
varImpPlot(rf_fit, main = "Variable Importance: Random Forest")
# Discussion: Compare the plots. Does Random Forest give more balanced importance
# scores? This often happens because it prevents one very strong predictor from
# dominating all the trees.

# =========================================================================
# --- Question 3: Fit a Boosted Tree Model ---
# =========================================================================
# Your task: Fit a GBM, evaluate its performance, and interpret its results.

# a) Prepare data (GBM for classification needs a 0/1 numeric response)
heart_gbm <- heart
heart_gbm$chd <- chd_original

# b) Fit the GBM model
# Key parameters:
# - distribution = "bernoulli": For logistic regression (0/1 outcome).
# - n.trees: The number of trees to grow. CV will find the best number.
# - shrinkage: The learning rate. Small values (e.g., 0.01) are typical.
# - interaction.depth: Controls the complexity of interactions.
set.seed(1)
gbm_fit <- gbm(chd ~ ., data = heart_gbm[train_indices, ],
               distribution = "bernoulli",
               n.trees = 5000,
               interaction.depth = 4,
               shrinkage = 0.01,
               cv.folds = 5,
               verbose = FALSE)

# c) Find the optimal number of trees using CV
best_iter_gbm <- gbm.perf(gbm_fit, method = "cv")
print(paste("Optimal number of trees for GBM:", best_iter_gbm))

# d) Evaluate performance on the test set
pred_gbm_prob <- predict(gbm_fit, newdata = test_data,
                         n.trees = best_iter_gbm, type = "response")
pred_gbm_class <- ifelse(pred_gbm_prob > 0.5, 1, 0)
accuracy_gbm <- mean(pred_gbm_class == chd_original[test_indices])

print(paste("Boosting Accuracy:", accuracy_gbm))

# e) Interpret the model
# Variable Importance
summary(gbm_fit, n.trees = best_iter_gbm)

# Partial Dependence Plots: Show the marginal effect of a variable on the
# response after accounting for other predictors.
par(mfrow = c(1, 2))
plot(gbm_fit, i.var = "tobacco", n.trees = best_iter_gbm)
plot(gbm_fit, i.var = "age", n.trees = best_iter_gbm)
# Discussion: These plots can reveal non-linear relationships. For example,
# the plot for 'age' might show that the risk of heart disease increases
# non-linearly after a certain age. This is a key advantage of boosting.

