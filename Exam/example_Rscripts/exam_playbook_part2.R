################################################################################
#
#        STATISTICAL METHODS FOR BIOINFORMATICS - PART II (PROF. JELIER)
#                     PRACTICAL EXAM STRATEGY SCRIPT
#
################################################################################

# This script is a general-purpose playbook for tackling the practical exam.
# The goal is not just to run code, but to demonstrate a systematic and
# thoughtful approach to analyzing a new dataset.

################################################################################
# S0: SETUP - LOAD LIBRARIES
################################################################################

# Load all the packages you might need at the beginning.
# The exam environment should have these pre-installed.

# Install the libraries if they are not already installed.
install.packages(c("glmnet", "mgcv", "randomForest", "gbm", "corrplot", "pROC"))

library(glmnet)       # For Ridge and LASSO
library(mgcv)         # For GAMs
library(randomForest) # For Random Forests
library(gbm)          # For Boosting
library(corrplot)     # For visualizing correlations
library(pROC)         # For AUC curves in classification

################################################################################
# S1: DATA TRIAGE & INITIAL EXPLORATION
#
# GOAL: Understand the data, frame the problem, and identify potential issues.
#       Your notes from this section are the material for the first part of
#       your answer (e.g., "Study and describe the data...").
################################################################################

# --- 1.1 Load and Frame the Problem ---

# Load the dataset (the exam will provide the file name)
# load("exam_data.Rdata")
# data <- read.csv("exam_data.csv")

# Get a first look at the data object
# data <- Vijver # Example from a lab

# Check dimensions: How many observations (n) and predictors (p)?
# Is it a high-dimensional problem (p > n)? This is a critical first question.
# If p > n, you MUST use regularization (LASSO/Ridge) or dimension reduction.
print("Dimensions of the data:")
# dim(data)

# Inspect variable names and structure
print("Variable names:")
# names(data)
print("Structure of the data:")
# str(data)

# Identify the response variable and the predictors based on the exam question.
# response_var <- "your_response_variable_name"
# predictor_vars <- names(data)[!names(data) %in% response_var]

# Determine the problem type: Regression or Classification?
# This dictates the models and evaluation metrics you will use.
# is.numeric(data[[response_var]]) # --> Regression
# is.factor(data[[response_var]]) || is.character(data[[response_var]]) # --> Classification


# --- 1.2 Checklist for Potential Issues ---

# Use summary() for a quick overview of ranges, missings (NA), etc.
print("Summary of the data:")
# summary(data)

# Check for outliers in key predictors
# boxplot(data$predictor1, main="Boxplot for Predictor 1")

# Check for correlation between predictors (crucial for LASSO vs. Ridge discussion)
# Note: model.matrix handles categorical variables for the correlation matrix
# x_matrix <- model.matrix(~ . -1, data = data[, predictor_vars])
# cor_matrix <- cor(x_matrix)
# corrplot(cor_matrix, method = "circle", type = "upper", tl.cex = 0.7)
# Finding: "I observe strong positive/negative correlations between predictors A and B..."

# Check for class imbalance (for CLASSIFICATION problems only)
# print("Class distribution of the response:")
# table(data[[response_var]])
# Finding: "The response variable shows significant class imbalance, with X% in class 0
#           and Y% in class 1. This means accuracy alone may be a misleading metric."


################################################################################
# S2: MODEL BUILDING TOOLKIT
#
# GOAL: Build and tune a few relevant models. LASSO and GAM are your primary
#       tools as they directly answer most of the likely exam questions.
################################################################################

# --- 2.1 Penalized Regression: LASSO (alpha=1) and Ridge (alpha=0) ---
# Use this to handle high-dimensionality and perform automatic variable selection.

# Prepare data for glmnet
# x <- model.matrix(response ~ ., data)
# y <- data$response

# Use cross-validation to find the best lambda (tuning parameter)
# Set alpha=1 for LASSO, alpha=0 for Ridge.
# Set family="binomial" for classification, "gaussian" for regression.

# print("Running cross-validation for LASSO...")
# set.seed(123) # for reproducibility
# cv_lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")

# Plot the CV error curve
# plot(cv_lasso)
# title("LASSO Cross-Validation", line = 2.5)

# Get the optimal lambda values
# lambda_min <- cv_lasso$lambda.min   # Lambda that gives minimum CV error
# lambda_1se <- cv_lasso$lambda.1se   # Simplest model within 1 SE of the minimum

# Inspect the coefficients of the final model
# print("LASSO Coefficients at lambda.1se:")
# coef(cv_lasso, s = lambda_1se)
# Finding: "The LASSO model with lambda chosen by the 1-standard-error rule selects
#           X non-zero coefficients, suggesting a sparse solution is possible."


# --- 2.2 Generalized Additive Models (GAMs) ---
# Use this to explicitly check for non-linear relationships.

# You can't fit a GAM with 100s of predictors.
# Choose a few important ones (e.g., from LASSO results or based on the problem description).

# print("Fitting a GAM model...")
# gam_fit <- gam(response ~ s(predictor1) + s(predictor2) + predictor3,
#                data = data, family = binomial) # or gaussian

# Check the summary for the Effective Degrees of Freedom (edf)
# summary(gam_fit)
# Finding: "The summary shows that the smooth term for `predictor1` has an edf of 3.8
#           and is highly significant (p < 0.001), providing strong evidence of a
#           non-linear relationship."

# Plot the smooth terms to visualize the relationships
# plot(gam_fit, pages = 1, all.terms = TRUE)
# Finding: "The plot for `predictor1` shows a U-shaped effect on the response."


# --- 2.3 Ensemble Methods (e.g., Boosting) ---
# Use this if you need a high-performance model for comparison or if interactions are suspected.

# print("Fitting a Gradient Boosting Machine (GBM)...")
# gbm_fit <- gbm(response ~ .,
#                data = data,
#                distribution = "bernoulli", # "gaussian" for regression
#                n.trees = 500,
#                interaction.depth = 4,
#                shrinkage = 0.01)

# summary(gbm_fit) # Shows variable importance


################################################################################
# S3: COMPARE MODELS AND ANSWER EXAM QUESTIONS
#
# GOAL: Synthesize your findings from S1 and S2 into coherent answers.
#       Structure your report around the questions provided.
################################################################################

# --- Question Template 1: Data Exploration & Potential Issues ---
# Use your findings from Step 1.
# "The dataset contains N observations and p predictors. This is a high-dimensional
# problem, which presents a risk of overfitting. The response is binary, indicating
# a classification task. A correlation plot revealed a strong positive correlation
# between predictors A and B. Boxplots showed potential outliers in predictor C.
# The response variable is imbalanced..."

# --- Question Template 2: Model Comparison & Interpretation ---

# "To model the data, a LASSO logistic regression and a GAM were trained..."

# Q: Is over-learning a problem?
# A: "Yes, over-learning is a major concern. To combat this, 10-fold cross-validation
#     was used to select the optimal tuning parameter (lambda) for the LASSO model.
#     This ensures the model is optimized for performance on unseen data, not just
#     the training set."

# Q: Do correlations between variables influence the results?
# A: "Yes. The strong correlation between predictors A and B was noted. The LASSO
#     model selected predictor A while setting the coefficient for B to zero. This is
#     a known characteristic of LASSO. A Ridge model, by contrast, would likely have
#     kept both predictors, shrinking their coefficients together."

# Q: Is there evidence for non-linear effects?
# A: "Yes. The GAM model provided strong evidence of non-linearity for predictor A.
#     The summary of the GAM fit showed its smooth term was highly significant with
#     an effective degrees of freedom of 4.1. The plot of this term revealed a
#     clear U-shaped relationship with the response."

# Q: Which model performs best?
# A: (This requires a formal comparison, often using cross-validation on all models)
#    "Comparing a 10-fold CV error, the LASSO model achieved an AUC of 0.85, while
#     the GAM achieved an AUC of 0.88, suggesting the non-linear relationships
#     captured by the GAM provide a slight performance edge. The boosted tree model
#     achieved the highest AUC of 0.91, but at the cost of interpretability."

# --- Question Template 3: Specific biological/domain question ---
# Use the results of your best/most interpretable model to answer.

# Q: Does Treatment X induce specific changes in protein levels?
# A: "To answer this, we used the LASSO model trained to predict Treatment. The
#     model selected a subset of 5 proteins with non-zero coefficients. The protein
#     with the largest coefficient was P53_N, suggesting it is the most significant
#     biomarker for distinguishing the treatments."