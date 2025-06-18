# =========================================================================
# GENERALIZED SCRIPT FOR EXAM-STYLE MULTILEVEL MODELING PROBLEMS
# =========================================================================
# This script provides a complete workflow for analyzing a multilevel model,
# focusing on the logical steps needed to build the right model and answer
# questions about it. We will use the clustered data example (rat pups),
# but the principles apply to longitudinal data as well.

# --- 1. Load Necessary Libraries ---
# -----------------------------------
# We need 'nlme' for the lme() function, which is essential for testing
# different variance structures (which the more modern lmer() cannot do).
library(nlme)
# We need 'RLRsim' for the most accurate way to test if a random effect is needed.
library(RLRsim)

# --- 2. Understand Your Data and Model ---
# ------------------------------------------
# First, identify the structure of the problem.
# - Outcome Variable (Level 1): e.g., 'weight' of a pup
# - Cluster Variable (Level 2): e.g., 'litterid' (the mother)
# - Predictors (can be Level 1 or 2): e.g., 'sex', 'treatment', 'litsize'

# Your initial, most complex plausible model might be:
# Y ~ predictors + (random_intercept | cluster)
# e.g., weight ~ treatment*sex1 + litsize + (1 | litterid)

# --- 3. The Analysis Strategy: Variance First, then Fixed Effects ---
# ---------------------------------------------------------------------
# The correct strategy is to first determine the best structure for the
# RANDOM part of your model using REML estimation. Once that is settled,
# you switch to ML estimation to test the FIXED effects.

# =========================================================================
# STEP A: FIND THE BEST RESIDUAL VARIANCE STRUCTURE (using REML)
# =========================================================================
# Your question: Is the within-cluster variance the same for all groups?
# (e.g., is the pup weight variability the same in all treatment groups?)

# Model 1: Homoscedastic Model (assumes one constant residual variance)
# We must use lme() and method="REML" for this part of the analysis.
# Note: In lme(), random effects are specified as 'random = ~1 | litterid'
model_hom <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1,
                 random = ~1 | litterid,
                 data = ratpup,
                 method = "REML")

# Model 2: Heteroscedastic Model (allows different variance per group)
# We add the 'weights' argument. This tells lme() to estimate a separate
# residual variance for each level of the 'treatment' variable.
model_het <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1,
                 random = ~1 | litterid,
                 weights = varIdent(form = ~ 1 | treatment), # This is the key change
                 data = ratpup,
                 method = "REML")

# Compare the two models using a Likelihood Ratio Test (LRT).
# This is valid with REML because the fixed effects part is identical.
anova(model_hom, model_het)

# --- How to Interpret the anova() Output ---
# It will produce a table with a p-value.
# - If p < 0.05: The more complex model (model_het) is a significantly better
#   fit. You should proceed using a heteroscedastic model.
# - If p > 0.05: The simpler model (model_hom) is sufficient. The variances
#   can be considered equal across groups.
#
# FOR THE RAT PUP DATA, THIS TEST IS HIGHLY SIGNIFICANT, SO WE CHOOSE 'model_het'.

# =========================================================================
# STEP B: TEST THE SIGNIFICANCE OF THE RANDOM EFFECT (using REML)
# =========================================================================
# Your question: Is the cluster effect necessary? (e.g., Is there a significant
# difference between litters after accounting for everything else?)
# This tests the hypothesis H₀: σ²_b = 0 (the variance of the random intercepts is zero).

# Model A: Our best model so far, which includes the random effect (model_het).
# Model B: The same model but with NO random effect. We use gls() for this.
model_no_random <- gls(weight ~ treatment + sex1 + litsize + treatment:sex1,
                       weights = varIdent(form = ~ 1 | treatment),
                       data = ratpup,
                       method = "REML")

# Compare them with an LRT. This p-value is "conservative" due to the boundary problem.
anova(model_no_random, model_het)

# To get a more accurate p-value that handles the boundary problem correctly:
exactRLRT(model_het) # Use your best LME model here.

# --- How to Interpret the Test Output ---
# - If p < 0.05: The random effect is significant and necessary. You must keep
#   it in your model.
# - If p > 0.05: The cluster effect is not significant. You could simplify your
#   analysis to the 'gls' model.
#
# FOR THE RAT PUP DATA, THIS TEST IS HIGHLY SIGNIFICANT, CONFIRMING THE LITTER EFFECT.

# =========================================================================
# STEP C: TEST THE FIXED EFFECTS (using ML)
# =========================================================================
# Now that we know our final variance structure is HETEROSCEDASTIC and includes a
# RANDOM LITTER EFFECT, we can test our main research questions.
# We MUST refit our chosen model with method = "ML" to fairly compare fixed effects.

final_model_ml <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1,
                      random = ~1 | litterid,
                      weights = varIdent(form = ~ 1 | treatment),
                      data = ratpup,
                      method = "ML") # Switched to ML for fair fixed effect tests

# To test the overall significance of each factor sequentially:
anova(final_model_ml)

# To see the specific coefficients, standard errors, and p-values for each term:
summary(final_model_ml)

# --- How to Interpret the Final summary() Output ---
# This is the table you use to answer most exam questions.

# Random effects:
#  Formula: ~1 | litterid
#          (Intercept) Residual
# StdDev:     0.2882     0.5123   <- Point estimates for σ_b and σ_ε

# Variance function: ... Parameter estimates:
#  Control       Low      High 
# 1.000000  0.589771  0.639438   <- Multipliers for residual variance

# Fixed effects:
#                             Value Std.Error  DF  t-value p-value
# (Intercept)               8.3506  0.26150   ...   31.93   0.0000  <- This is γ₀₀
# treatmentHigh            -0.9047  0.18092   ...   -5.00   0.0000  <- Effect of High vs Control
# treatmentLow             -0.4668  0.15105   ...   -3.09   0.0052  <- Effect of Low vs Control
# sex1                     -0.4065  0.09357   ...   -4.34   0.0000  <- Effect of Female vs Male
# litsize                  -0.1304  0.01755   ...   -7.42   0.0000  <- Effect of Litter Size
# treatmentHigh:sex1        0.0930  0.12521   ...    0.74   0.4581  <- Interaction Term
# treatmentLow:sex1         0.0756  0.10998   ...    0.68   0.4924  <- Interaction Term
#
# Use this table to calculate predicted values and determine the significance,
# magnitude, and direction of the effects of your predictors.
