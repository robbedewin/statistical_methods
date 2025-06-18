# --- A General Workflow for Multilevel Model Analysis ---
# This script guides you through the key steps for analyzing a multilevel
# model, from model building to testing the variance structure and fixed effects.

# 1. LOAD NECESSARY PACKAGES
# ---------------------------
# nlme: For the lme() function, which is necessary for modeling unequal variances.
# lmerTest: For getting p-values for fixed effects in lmer() models.
# RLRsim: For correctly testing if a random effect variance is zero.
library(nlme)
library(lmerTest) # For lmer(), though we focus on lme() here
library(RLRsim)   # For the exact Restricted Likelihood Ratio Test

# It's good practice to simulate data to understand the structure
# For this example, we assume you have the 'ratpup' dataframe loaded.
# Let's make sure the factor levels are set correctly for interpretation.
# ratpup$treatment <- relevel(ratpup$treatment, ref = "Control")
# ratpup$sex <- as.factor(rat_pup$sex) # if not already a factor

# 2. MODEL BUILDING & STRATEGY
# ----------------------------
# The general strategy is to find the best model for the VARIANCE structure first,
# using REML estimation. Then, switch to ML estimation to test the FIXED EFFECTS.

# --- STEP A: TEST THE VARIANCE OF THE RESIDUALS (HETEROSCEDASTICITY) ---
# Question: Is the within-group variance the same across all treatment groups?

# Model 1: Homoscedastic (assumes one, constant residual variance)
# We use lme() and method = "REML" because we are testing variance components.
model.hom <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1,
                 random = ~1 | litterid,
                 data = ratpup,
                 method = "REML")

# Model 2: Heteroscedastic (allows different residual variance per treatment)
# We add the 'weights' argument to specify this.
model.het <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1,
                 random = ~1 | litterid,
                 weights = varIdent(form = ~ 1 | treatment), # Key difference
                 data = ratpup,
                 method = "REML")

# Compare the two models using a Likelihood Ratio Test (LRT)
# This is valid with REML because the fixed effects are identical.
anova(model.hom, model.het)

# INTERPRETATION:
# If p-value is small (< 0.05), the heteroscedastic model (model.het) is
# a significantly better fit. You should use it for the next steps.
# If p-value is large, the simpler homoscedastic model is sufficient.
# (For the ratpup data, the p-value is very small, so we proceed with model.het)


# --- STEP B: TEST THE SIGNIFICANCE OF THE RANDOM EFFECT (CLUSTER EFFECT) ---
# Question: Is there a significant litter effect? (i.e., is the variance of the
# random intercepts, σ²_b, significantly greater than zero?)

# Model A: The best model so far (model.het, which includes the random effect).
# Model B: A model with NO random effect. We use gls() for this.
model.no_random <- gls(weight ~ treatment + sex1 + litsize + treatment:sex1,
                       weights = varIdent(form = ~ 1 | treatment), # Keep the weights
                       data = ratpup,
                       method = "REML")

# Compare them with an LRT.
# NOTE: This p-value is "conservative" due to the boundary problem.
anova(model.no_random, model.het)

# For a more accurate test that accounts for the boundary problem:
exactRLRT(model.het)

# INTERPRETATION:
# If p-value is small (< 0.05), the random effect is significant. You need it
# in your model. The litter-to-litter variation is real.
# If p-value is large, you could simplify the model to the gls() version.
# (For the ratpup data, the p-value is tiny, so we confirm we need the random effect).


# --- STEP C: TEST THE FIXED EFFECTS ---
# Now that we've confirmed our variance structure (heteroscedastic residuals
# AND a random litter effect), we can test our main research questions.
# To do this fairly, we MUST refit our final model using method = "ML".

# Final Model using Maximum Likelihood (ML)
final.model.ml <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1,
                      random = ~1 | litterid,
                      weights = varIdent(form = ~ 1 | treatment),
                      data = ratpup,
                      method = "ML") # Switched to ML

# Now we can test the fixed effects using anova() on this single model.
# This performs a sequential (Type I) test.
anova(final.model.ml)

# INTERPRETATION of anova(final.model.ml):
# - treatment p-value: Is treatment significant on its own?
# - sex1 p-value: Is sex significant after accounting for treatment?
# - litsize p-value: Is litter size significant after accounting for treatment and sex?
# - treatment:sex1 p-value: Is the interaction significant after accounting for all main effects?
# (From the lecture, this interaction was not significant, so it could be dropped
# for a final, more parsimonious model).

# Finally, look at the summary of the final ML model to get the coefficients.
summary(final.model.ml)

# INTERPRETATION of summary():
# This is where you get all the numbers (point estimates) to calculate predicted
# values and answer specific questions about the magnitude and direction of effects,
# just like we did for your exam questions.