################################################################################
#
#        EXAM SCRIPT: SOLVING THE EDUCATIONAL STUDY PROBLEM (ALONSO)
#
################################################################################

# This script provides a complete walkthrough for the 2024 exam question.
# It covers:
# 1. Loading and preparing the data.
# 2. Translating the mathematical model from the question into R code.
# 3. Fitting the specified hierarchical model.
# 4. A detailed guide on how to interpret the output to answer each
#    multiple-choice question systematically.

################################################################################
# STEP 0: SETUP - LOAD REQUIRED PACKAGES
################################################################################

# The lme4 package is essential for fitting linear mixed-effects models.
# The lmerTest package is useful as it adds p-values to the summary output.
# You may need to install the package if it's not available in the exam environment:
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("readxl") # If the data is in an Excel file

library(lme4)
library(lmerTest)
library(readxl)


################################################################################
# STEP 1: LOAD AND INSPECT THE DATA
################################################################################

# The exam will provide the data file. Let's assume it's an Excel file.
# Replace "path/to/your/file.xlsx" with the actual file path.
# educational_data <- read_excel("path/to/your/file.xlsx")

# For demonstration, let's create a simulated dataset that matches the description
set.seed(42) # for reproducibility
n_students <- 200
time_points <- 0:5
student_ids <- 1:n_students

# Create the base data frame
educational_data <- expand.grid(time = time_points, id = student_ids)

# Assign students to programs (randomized)
programs <- data.frame(id = student_ids, prog = rbinom(n_students, 1, 0.5))
educational_data <- merge(educational_data, programs, by = "id")

# --- Define the true parameters for simulation ---
# These are the "real" values we are trying to estimate.
# Fixed effects
gamma_00 <- 5.0  # Avg. initial score for standard program (prog=0)
gamma_10 <- 0.5  # Avg. slope for standard program
gamma_01 <- -0.2 # Difference in initial score (should be near 0 due to randomization)
gamma_11 <- 1.5  # EXTRA increase in slope for new program (prog=1)

# Random effects and residual variance
sigma_b0 <- 2.0  # SD of random intercepts
sigma_eps <- 3.0 # SD of residuals

# Generate the data based on the model
b_0i <- rnorm(n_students, mean = 0, sd = sigma_b0)
random_effects <- data.frame(id = student_ids, b_0i = b_0i)
educational_data <- merge(educational_data, random_effects, by = "id")

# Calculate the response variable 'y'
educational_data$y <- with(educational_data,
                           (gamma_00 + gamma_01*prog + b_0i) +        # Intercept part
                             (gamma_10 + gamma_11*prog) * time +       # Slope part
                             rnorm(nrow(educational_data), 0, sigma_eps) # Residual error
)
# Clean up extra columns for clarity
educational_data <- educational_data[, c("id", "y", "time", "prog")]


# --- Initial Inspection ---
# Always a good first step.
head(educational_data)
str(educational_data)

View(educational_data)
################################################################################
# STEP 2: FORMULATE AND FIT THE HIERARCHICAL MODEL
################################################################################

# The question explicitly provides the model. Our task is to translate it to R code.

# Level 1 Model: Y_ij = pi_0i + pi_1i * time_j + epsilon_ij
# Level 2 Model:
#   pi_0i = gamma_00 + gamma_01*PROG_i + b_0i
#   pi_1i = gamma_10 + gamma_11*PROG_i
#
# This model specifies:
# - FIXED effects for the intercept, time, program, and the time:program interaction.
# - A RANDOM effect for the intercept (b_0i). Each student has their own baseline.
# - NO random effect for the slope. The slope is assumed to be the same for all
#   students within the same program group.

# Translating to lmer() formula:
# Fixed effects: y ~ time * prog  (This expands to intercept + time + prog + time:prog)
# Random effects: (1 | id)        (This specifies a random intercept for each student 'id')

# The question requires fitting using Maximum Likelihood (ML), not the default REML.
# Therefore, we must add the argument `REML = FALSE`.

# Fit the model
model_fit <- lmer(y ~ time * prog + (1 | id), data = educational_data, REML = FALSE)

#Fit model from exam
long_teaching <- read_csv("~/KULeuven/Bioinformatics_2024-2025/Statistical_Methods/Exam/long_teaching.txt")
View(long_teaching)
attach(long_teaching)

model_fit <- lmer(y ~ time * prog + (1 |id), data = long_teaching, REML = FALSE)

################################################################################
# STEP 3: INTERPRET THE OUTPUT TO ANSWER THE QUESTIONS
################################################################################

# Get the detailed summary of our fitted model
model_summary <- summary(model_fit)
print(model_summary)


# --- Question 1: Compare variability between vs. within individuals ---
# "Based on the point estimates obtained from the previous model which of the
#  following statements is correct?"
#  - The variability between individuals is larger/smaller/equal to the variability within individuals.

# HOW TO ANSWER:
# Look at the "Random effects" section of the summary.
# - "Variability between individuals" corresponds to the variance of the random intercept.
#   This is the `Variance` for `id (Intercept)`.
# - "Variability within individuals" corresponds to the residual variance.
#   This is the `Variance` for `Residual`.

var_between <- VarCorr(model_fit)$id[1,1]
var_within <- sigma(model_fit)^2

print(paste("Variance BETWEEN individuals (sigma_0^2):", round(var_between, 2)))
print(paste("Variance WITHIN individuals (sigma_epsilon^2):", round(var_within, 2)))

if (var_between > var_within) {
  print("Conclusion 1: The variability between individuals is larger than the variability within individuals.")
} else {
  print("Conclusion 1: The variability between individuals is smaller than the variability within individuals.")
}


# --- Question 2: Is the new program significantly better? ---
# This question is about the RATE OF INCREASE over time. Is the slope for the new
# program significantly steeper than the slope for the standard program?

# HOW TO ANSWER:
# Look at the "Fixed effects" section for the interaction term `time:prog`.
# This coefficient represents the *difference in slopes* between the new and standard programs (gamma_11).
# A positive and significant coefficient means the new program produces a faster increase.

# Extract the interaction term's results
interaction_term <- coef(model_summary)["time:prog", ]
estimate_g11 <- interaction_term["Estimate"]
p_value_g11 <- interaction_term["Pr(>|t|)"]

print(paste("Estimate for gamma_11 (time:prog):", round(estimate_g11, 2)))
print(paste("P-value for gamma_11 (time:prog):", round(p_value_g11, 4)))


if (estimate_g11 > 0 & p_value_g11 < 0.05) {
  print("Conclusion 2: As an average the new program is significantly better than the standard program.")
} else if (estimate_g11 > 0 & p_value_g11 >= 0.05) {
  print("Conclusion 2: As an average based on the point estimates the new program seems to be better than the standard program but the result is not significant.")
} else {
  # Handle other cases if necessary
}


# --- Question 3: What is the main effect of the program? ---
# This question is about the `prog` coefficient (gamma_01), which represents the
# difference in the starting point (at time=0) between the two programs.

# HOW TO ANSWER:
# Look at the "Fixed effects" section for the `prog` term. Check its sign and significance.
# Because of randomization, we expect this effect to be non-significant.

main_effect_prog <- coef(model_summary)["prog", ]
estimate_g01 <- main_effect_prog["Estimate"]
p_value_g01 <- main_effect_prog["Pr(>|t|)"]

print(paste("Estimate for gamma_01 (prog):", round(estimate_g01, 2)))
print(paste("P-value for gamma_01 (prog):", round(p_value_g01, 4)))

if (p_value_g01 >= 0.05) {
  print(paste("Conclusion 3: The point estimate for the main effect of program is",
              ifelse(estimate_g01 > 0, "positive", "negative"), "and not significant."))
} else {
  # Handle other cases
}


# --- Question 4: Calculate the average SAT score for a specific case ---
# "Based on the previous model what is the average SAT score for children
#  in the standard program one month after starting the experiment?"

# HOW TO ANSWER:
# Use the fixed effects estimates to calculate the predicted value.
# The model for the average response is:
# E(Y) = gamma_00 + gamma_01*prog + gamma_10*time + gamma_11*prog*time

# For this question: prog = 0 (standard program) and time = 1.

coeffs <- fixef(model_fit)
gamma_00_hat <- coeffs["(Intercept)"]
gamma_01_hat <- coeffs["prog"]
gamma_10_hat <- coeffs["time"]
gamma_11_hat <- coeffs["time:prog"]

predicted_score <- gamma_00_hat + (gamma_01_hat * 0) + (gamma_10_hat * 1) + (gamma_11_hat * 0 * 1)

print(paste("Conclusion 4: The predicted average score is:", round(predicted_score, 2)))


# --- Question 5: Assess the p-value for gamma_01 ---
# "Based on the previous model which of the following statements is correct?"
# This asks for the range of the p-value for gamma_01 (the main effect of `prog`).

# HOW TO ANSWER:
# We already extracted this value for Question 3. We just need to check its range.

print(paste("P-value for gamma_01 (prog) is:", round(p_value_g01, 4)))

if (p_value_g01 < 0.05) {
  print("Conclusion 5: The p-value associated with gamma_01 is smaller than 5%")
} else if (p_value_g01 >= 0.05 & p_value_g01 <= 0.80) {
  print("Conclusion 5: The p-value associated with gamma_01 is between 5% and 80%")
} else if (p_value_g01 > 0.90) {
  print("Conclusion 5: The p-value associated with gamma_01 is larger than 90%")
} else {
  print("Conclusion 5: The p-value is in another range.")
}

