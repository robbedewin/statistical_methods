# Ensure necessary libraries are installed
# %pip install seaborn matplotlib pandas scikit-learn pygam statsmodels nbformat

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
# from sklearn.preprocessing import StandardScaler, LabelEncoder
# from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
# from sklearn.linear_model import LogisticRegressionCV
# from pygam import LogisticGAM, s, l, te
# from sklearn.ensemble import GradientBoostingClassifier
# from sklearn.metrics import roc_auc_score, confusion_matrix, roc_curve
# from sklearn.inspection import PartialDependenceDisplay
# from scipy.stats import ttest_ind
# import statsmodels.stats.multitest as smm
import nbformat

# Create a new notebook object
nb = nbformat.v4.new_notebook()

# Helper function to add code cells to the notebook
def add_code_cell(code_string):
    nb.cells.append(nbformat.v4.new_code_cell(code_string))

# Helper function to add markdown cells to the notebook
def add_markdown_cell(markdown_string):
    nb.cells.append(nbformat.v4.new_markdown_cell(markdown_string))

# --- Notebook Title and Introduction ---
add_markdown_cell("""\
# Analysis of Protein Expression in Control and Down Syndrome Mice

This notebook analyzes a dataset containing protein expression levels in the cerebral cortex of mice. The goal is to understand the effects of genotype, behavior, and drug treatment (Memantine vs. Saline) on these protein levels.

The analysis addresses the following questions:
1.  Study and describe the data. Do you see indications of potential issues when statistically modeling the data?
2.  Train and compare LASSO, GAM, and a boosting model to separate Memantine from Saline treated samples. Interpret the results, discussing correlations, over-learning, non-linear effects, and interactions.
3.  Evaluate if Memantine treatment induces any specific changes in protein levels of trisomic mice who were stimulated to learn.
""")

# --- Setup and Data Loading ---
add_markdown_cell("## 0. Setup and Data Loading")
add_code_cell("""\
# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.linear_model import LogisticRegressionCV # For LASSO
from pygam import LogisticGAM, s # For GAM
from sklearn.ensemble import GradientBoostingClassifier # For Boosting
from sklearn.metrics import roc_auc_score, confusion_matrix, roc_curve
from sklearn.inspection import PartialDependenceDisplay
from scipy.stats import ttest_ind # For t-tests
import statsmodels.stats.multitest as smm # For FDR correction

# Set plotting style
sns.set(style='whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)

# Load the dataset
# This assumes 'Cortex3.csv' is in the same directory as the notebook.
# If not, please adjust the path.
try:
    df = pd.read_csv("Cortex3.csv")
except FileNotFoundError:
    print("Error: Cortex3.csv not found. Please ensure the file is in the correct directory.")
    raise

print("Data Shape (samples, columns):", df.shape)
print("\\nFirst 5 rows of the data:")
print(df.head())

# Identify protein columns (ending with '_N') and categorical/metadata columns
protein_cols = [col for col in df.columns if col.endswith('_N')]
metadata_cols = ['Genotype', 'Treatment', 'Behavior']
print(f"\\nFound {len(protein_cols)} protein expression columns.")
print(f"Found {len(metadata_cols)} metadata columns: {metadata_cols}")
""")

# --- Question 1: Data Study and Description ---
add_markdown_cell("""\
## 1. Study and Describe the Data

This section explores the dataset to understand its structure, identify potential issues, and prepare for statistical modeling.
""")

add_markdown_cell("### 1.1. Initial Data Inspection")
add_code_cell("""\
# Purpose: Check for missing data and understand the data types of columns.
# This helps ensure data quality and informs necessary preprocessing steps.

# Check for missing values
print("\\nMissing values per column (showing columns with any missing values):")
missing_values = df.isnull().sum()
print(missing_values[missing_values > 0])
if missing_values.sum() == 0:
    print("No missing values found in the dataset.")

# Check data types
print("\\nData types summary:")
print(df.dtypes.value_counts())
print("\\nData types of metadata columns:")
print(df[metadata_cols].dtypes)
print("\\nData types of first 3 protein columns:")
print(df[protein_cols[:3]].dtypes)

# Convert categorical columns from 'object' to 'category' type for efficiency and semantic correctness.
for col in metadata_cols:
    if df[col].dtype == 'object':
        df[col] = df[col].astype('category')
print("\\nData types of metadata columns after potential conversion:")
print(df[metadata_cols].dtypes)
""")

add_markdown_cell("### 1.2. Descriptive Statistics")
add_code_cell("""\
# Purpose: Get a quantitative overview of the data.
# Descriptive statistics for numerical columns help understand their scale, central tendency, and spread.
# Value counts for categorical columns show the distribution of samples across different classes.

# Summary statistics for numerical (protein) columns
print("\\nSummary statistics for protein expression levels (first 5 proteins):")
print(df[protein_cols].describe().transpose().head())

# Value counts for categorical columns
print("\\nDistribution of samples across Genotype:")
print(df['Genotype'].value_counts())
print("\\nDistribution of samples across Treatment:")
print(df['Treatment'].value_counts())
print("\\nDistribution of samples across Behavior:")
print(df['Behavior'].value_counts())
""")

add_markdown_cell("### 1.3. Visual Data Exploration")
add_markdown_cell("#### 1.3.1. Distribution of Protein Expression Levels")
add_code_cell("""\
# Purpose: Visualize the distribution of individual protein expression levels.
# Histograms help identify skewness, multimodality, or other patterns.
# Boxplots help compare distributions across different groups (e.g., Treatment).

# Plot histograms for a few example protein expression levels
print("\\nHistograms for the first 3 protein expression levels:")
plt.figure(figsize=(18, 5))
for i, col in enumerate(protein_cols[:3]): # Plot first 3 proteins
    plt.subplot(1, 3, i+1)
    sns.histplot(df[col], kde=True, bins=15)
    plt.title(f'Distribution of {col}')
plt.tight_layout()
plt.show()

# Boxplots for a few example proteins grouped by 'Treatment'
print("\\nBoxplots for the first 3 protein expression levels by Treatment:")
plt.figure(figsize=(18, 5))
for i, col in enumerate(protein_cols[:3]):
    plt.subplot(1, 3, i+1)
    sns.boxplot(x='Treatment', y=col, data=df)
    plt.title(f'{col} by Treatment')
plt.tight_layout()
plt.show()
""")

add_markdown_cell("#### 1.3.2. Protein Expression Correlation Heatmap")
add_code_cell("""\
# Purpose: Visualize the correlation structure between protein expression levels.
# A heatmap helps identify groups of correlated proteins (multicollinearity),
# which can affect model coefficient stability and interpretation.

# Calculate and plot the correlation matrix for protein expression levels
corr_matrix = df[protein_cols].corr()
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix, cmap="coolwarm", center=0, annot=False, xticklabels=False, yticklabels=False)
plt.title("Protein Expression Correlation Heatmap")
plt.show()
""")

add_markdown_cell("#### 1.3.3. Principal Component Analysis (PCA)")
add_code_cell("""\
# Purpose: Perform dimensionality reduction to visualize sample-to-sample variation.
# PCA helps check for outliers, batch effects (if applicable), or gross separation by
# metadata variables like Treatment, Genotype, or Behavior in a lower-dimensional space.
# It also indicates how much variance is captured by the first few principal components.

# Scale protein expression data before PCA
X_proteins_scaled = StandardScaler().fit_transform(df[protein_cols])

# Perform PCA
from sklearn.decomposition import PCA
# Using enough components to explain a good portion of variance, e.g., 5-10, then visualize first 2.
pca_full = PCA(n_components=min(10, X_proteins_scaled.shape[1])) 
pca_full.fit(X_proteins_scaled)

pca_2comp = PCA(n_components=2)
pcs_2comp = pca_2comp.fit_transform(X_proteins_scaled)

pc_df = pd.DataFrame(data=pcs_2comp, columns=['PC1', 'PC2'])
pc_df_combined = pd.concat([pc_df, df[metadata_cols].reset_index(drop=True)], axis=1) # reset_index if df index is not standard

print(f"\\nExplained variance by PC1: {pca_2comp.explained_variance_ratio_[0]:.3f}")
print(f"Explained variance by PC2: {pca_2comp.explained_variance_ratio_[1]:.3f}")
print(f"Total explained variance by first 2 PCs: {np.sum(pca_2comp.explained_variance_ratio_):.3f}")
print(f"Cumulative explained variance by first 10 PCs (or fewer if <10 features): {np.cumsum(pca_full.explained_variance_ratio_)}")


plt.figure(figsize=(10, 7))
sns.scatterplot(x='PC1', y='PC2', hue='Treatment', style='Genotype', size='Behavior', data=pc_df_combined, s=100, alpha=0.8)
plt.title('PCA of Cortex Proteome (First Two Principal Components)')
plt.xlabel(f'PC1 ({pca_2comp.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca_2comp.explained_variance_ratio_[1]*100:.1f}%)')
plt.legend(title='Groups', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()
""")

add_markdown_cell("""\
### 1.4. Indications of Potential Issues for Statistical Modeling

Based on the data exploration:

1.  **High Dimensionality (Curse of Dimensionality):**
    * The dataset has 70 samples and 70 protein expression features (predictors). This $p \\approx n$ scenario is a classic high-dimensionality problem. [cite: 4]
    * **Issue:** This increases the risk of overfitting models. Models might learn noise specific to the training data and perform poorly on unseen data. It also makes feature selection crucial. [cite: 8] Scaling variables, as done for PCA and modeling, is also important due to varying ranges of protein expression levels. [cite: 63]

2.  **Multicollinearity:**
    * The correlation heatmap shows blocks of correlated proteins.
    * **Issue:** High multicollinearity can make model coefficients unstable and difficult to interpret for linear models like LASSO. LASSO performs feature selection but might arbitrarily pick one from a group of highly correlated features. Tree-based models like Gradient Boosting are less affected in terms of predictive accuracy, but feature importance might be distributed.

3.  **Small Sample Size per Class/Subgroup:**
    * The 'Treatment' variable has 38 Memantine and 32 Saline samples, which is reasonably balanced for overall classification.
    * However, for Question 3, the specific subgroup (Ts65Dn mice, C/S behavior) has only 9 Memantine and 7 Saline samples.
    * **Issue:** Small sample sizes reduce statistical power, making it harder to detect true effects. Results from analyses on small subgroups should be interpreted with caution. This is particularly relevant for the t-tests in Q3.

4.  **Potential Outliers and Distributions:**
    * Histograms for some proteins might show skewed distributions.
    * **Issue:** Outliers can unduly influence some models. Non-normal distributions might affect assumptions of certain statistical tests (like t-tests in Q3). For models like GAMs, the flexibility of splines can handle some non-normality.

5.  **Model Selection Complexity:**
    * With $p \\approx n$, selecting the "best" model requires careful validation to avoid overfitting and ensure generalizability. Techniques like cross-validation are essential. [cite: 11, 35]

The PCA plot showed some, but not complete, separation by Treatment and Genotype, suggesting that protein expression levels contain relevant information, but predictive modeling will be necessary to achieve good classification.
""")

# --- Question 2: Model Comparison (LASSO, GAM, Boosting) ---
add_markdown_cell("""\
## 2. Train and Compare LASSO, GAM, and Boosting Models

This section focuses on training and comparing three different models (LASSO, GAM, Gradient Boosting) to classify samples based on their treatment (Memantine vs. Saline) using protein expression data. We will use 5-fold cross-validation to estimate the test error (AUC) for model comparison, as discussed in lectures. [cite: 11, 35]
""")

add_markdown_cell("### 2.1. Data Preparation for Classification")
add_code_cell("""\
# Purpose: Prepare the data for classification tasks.
# This involves defining the feature matrix (X) and target variable (y),
# encoding the categorical target variable into numerical format,
# and scaling the features, which is crucial for LASSO. [cite: 63]

# Define feature matrix X (protein expressions) and target y (Treatment)
X = df[protein_cols].values
# Encode 'Treatment' labels: Memantine=1, Saline=0
# This ensures consistency with how many classifiers expect binary outcomes.
label_encoder = LabelEncoder()
y = label_encoder.fit_transform(df['Treatment']) 
# Confirm encoding: Memantine should be 1 as per original notebook comment
print(f"Original Treatment labels: {df['Treatment'].unique()}")
print(f"Encoded Treatment labels: {np.unique(y)} maps to classes: {label_encoder.classes_}")
if label_encoder.classes_[1] != 'Memantine': # Assuming 'Memantine' is the positive class
    print("Re-encoding 'Memantine' to 1 and 'Saline' to 0.")
    y = (df['Treatment'] == 'Memantine').astype(int)
    # Create a consistent mapping for later reference if needed
    # label_encoder.classes_ = np.array(['Saline', 'Memantine']) 

# Scale the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split data into training and testing sets 
# This split is for an initial evaluation of models after parameter tuning (for LASSO)
# and for generating PDPs. The primary model comparison will use 5-fold CV on the full (scaled) dataset.
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.3, stratify=y, random_state=42)

print(f"\\nX_train shape: {X_train.shape}, y_train shape: {y_train.shape}")
print(f"X_test shape: {X_test.shape}, y_test shape: {y_test.shape}")
print(f"Class distribution in y_train: Saline (0)={np.sum(y_train==0)}, Memantine (1)={np.sum(y_train==1)}")
print(f"Class distribution in y_test: Saline (0)={np.sum(y_test==0)}, Memantine (1)={np.sum(y_test==1)}")
""")

add_markdown_cell("### 2.2. Model Training and Evaluation on Test Split")
add_markdown_cell("#### 2.2.1. LASSO (Logistic Regression with L1 penalty)")
add_code_cell("""\
# Purpose: Train a LASSO model. LogisticRegressionCV is used, which internally
# performs cross-validation to find the best regularization parameter C (inverse of lambda).
# The L1 penalty promotes sparsity, performing feature selection. [cite: 48, 66]

lasso_model_cv = LogisticRegressionCV(
    Cs=np.logspace(-4, 2, 20),  # Explore a range of C values
    cv=5,                      # 5-fold internal cross-validation for C
    penalty='l1',
    solver='saga',             # 'saga' or 'liblinear' are suitable for L1
    scoring='roc_auc',
    max_iter=10000,
    random_state=42,
    n_jobs=-1
)
lasso_model_cv.fit(X_train, y_train) # Fit on the training data to find best C

print(f"Best C for LASSO (from training data CV): {lasso_model_cv.C_[0]:.4f}")

# Evaluate on the held-out test set
y_prob_lasso_test = lasso_model_cv.predict_proba(X_test)[:, 1]
auc_lasso_test = roc_auc_score(y_test, y_prob_lasso_test)
print(f"LASSO Test AUC: {auc_lasso_test:.3f}")

# Inspect coefficients from the model trained with the best C
lasso_coeffs = pd.Series(lasso_model_cv.coef_[0], index=protein_cols)
selected_features_lasso = lasso_coeffs[lasso_coeffs != 0]
print(f"\\nNumber of features selected by LASSO: {len(selected_features_lasso)}")
if len(selected_features_lasso) > 0:
    print("Top 5 selected features by LASSO (absolute coefficient value):")
    print(selected_features_lasso.abs().sort_values(ascending=False).head())
else:
    print("No features selected by LASSO with the chosen C values.")

# Plot coefficient paths for LogisticRegressionCV
# coefs_paths_ stores coefficients for each class, for each C, for each fold.
# We'll average over folds for each C.
# Shape of coefs_paths_: (n_classes, n_folds, n_features, n_cs) if refit=False, or (n_classes, n_features, n_cs) if refit=True (default for L1)
# For binary classification, coefs for class 1 are usually of interest.
coefs_path = lasso_model_cv.coefs_paths_[1] # Assuming class 1 (Memantine) coefs
if coefs_path.ndim == 3: # if per-fold paths are stored
    coefs_path_mean = np.mean(coefs_path, axis=0) # Average over folds
else: # if paths are for the model refit on full train data with best C for each C
    coefs_path_mean = coefs_path 

plt.figure(figsize=(12, 7))
# Each column in coefs_path_mean.T corresponds to a feature's path across Cs
# LogisticRegressionCV uses Cs (inverse of lambda). Plot against 1/Cs for lambda-like behavior.
lambdas_lasso = 1. / lasso_model_cv.Cs_
sorted_lambda_indices = np.argsort(lambdas_lasso) # Sort for plotting

for i in range(coefs_path_mean.shape[0]):
    plt.plot(lambdas_lasso[sorted_lambda_indices], coefs_path_mean[i, sorted_lambda_indices], label=f'{protein_cols[i]}' if i < 5 else None) # Label first 5

plt.xscale('log')
plt.xlabel('Lambda (Inverse of C, Log Scale)')
plt.ylabel('Coefficient Value')
plt.title('LASSO Coefficient Paths')
plt.axvline(1./lasso_model_cv.C_[0], color='red', linestyle='--', label=f'Optimal Lambda (1/C_best = {1./lasso_model_cv.C_[0]:.4f})')
# plt.legend(loc='upper right', bbox_to_anchor=(1.25, 1)) # Adjust legend if too many features
plt.legend()
plt.show()

""")

add_markdown_cell("#### 2.2.2. Generalized Additive Model (GAM)")
add_code_cell("""\
# Purpose: Train a GAM to capture potential non-linear relationships between
# protein expression and treatment outcome.
# Due to high dimensionality, features selected by LASSO are used for GAM.

if len(selected_features_lasso) == 0:
    print("LASSO selected 0 features. GAM training will be skipped or use a fallback.")
    # As a fallback, let's select top 5 features based on Gradient Boosting importance if available later,
    # or simply skip GAM if no robust way to select a small feature set.
    # For now, we'll print a message and handle it in the CV section.
    gam_feature_indices_for_fit = []
    gam_feature_names_for_fit = []
    gam_model_fit = None
    auc_gam_test = np.nan
else:
    gam_feature_names_for_fit = selected_features_lasso.index.tolist()
    gam_feature_indices_for_fit = [protein_cols.index(name) for name in gam_feature_names_for_fit]

    print(f"\\nGAM will be trained on {len(gam_feature_names_for_fit)} features selected by LASSO: {gam_feature_names_for_fit[:5]}...")

    # Construct terms for GAM: one spline term 's(i)' for each selected feature.
    gam_terms_fit = s(0) 
    for i in range(1, len(gam_feature_names_for_fit)):
        gam_terms_fit += s(i)

    gam_model_fit = LogisticGAM(gam_terms_fit)
    gam_model_fit.fit(X_train[:, gam_feature_indices_for_fit], y_train)

    # Evaluate GAM on the test set
    y_prob_gam_test = gam_model_fit.predict_proba(X_test[:, gam_feature_indices_for_fit])
    auc_gam_test = roc_auc_score(y_test, y_prob_gam_test)

print(f"GAM Test AUC (on LASSO-selected features): {auc_gam_test:.3f}")

# It's also good practice to check GAM summary for effective degrees of freedom (EDoF) of splines.
# gam_model_fit.summary() # This would print a text summary.
""")

add_markdown_cell("#### 2.2.3. Gradient Boosting Machine (GBM)")
add_code_cell("""\
# Purpose: Train a Gradient Boosting model, which can capture complex non-linear
# relationships and feature interactions.

gb_model_fit = GradientBoostingClassifier(random_state=42, n_estimators=100, max_depth=3, learning_rate=0.1)
# Note: n_estimators and other hyperparameters ideally should be tuned via CV (e.g., GridSearchCV).
# For this assignment, fixed reasonable values are used for initial evaluation.
gb_model_fit.fit(X_train, y_train)

# Evaluate on the test set
y_prob_gb_test = gb_model_fit.predict_proba(X_test)[:, 1]
auc_gb_test = roc_auc_score(y_test, y_prob_gb_test)
print(f"Gradient Boosting Test AUC: {auc_gb_test:.3f}")

# Feature importances
gb_importances = pd.Series(gb_model_fit.feature_importances_, index=protein_cols)
top10_gb_features = gb_importances.nlargest(10)

plt.figure(figsize=(8, 6))
sns.barplot(x=top10_gb_features.values, y=top10_gb_features.index)
plt.title('GBM Feature Importance (Top 10 from training split model)')
plt.xlabel('Importance')
plt.show()
""")

add_markdown_cell("### 2.3. Model Comparison using 5-Fold Cross-Validation AUC")
add_markdown_cell("""\
# Purpose: More robustly compare the generalization performance of LASSO, GAM, and GBM.
# We use 5-fold stratified cross-validation on the *entire scaled dataset* (X_scaled, y)
# and report the mean and standard deviation of the AUC scores.
# Stratification helps ensure class proportions are similar across folds.
""")
add_code_cell("""\
cv_stratified = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# --- LASSO CV ---
# For LASSO, we create a new LogisticRegressionCV instance to be used in cross_val_score.
# This ensures that for each fold of the outer CV, the hyperparameter C is tuned
# using an inner CV on that fold's training data. This is a more rigorous approach.
# Alternatively, one could fix C based on lasso_model_cv.C_[0] and CV that fixed model.
# The lecture notes [cite: 59] emphasize CV for lambda (or C) selection.
# The LogisticRegressionCV itself handles this internal CV for lambda.
# So, cross_val_score on a LogisticRegressionCV instance effectively does nested CV.
lasso_for_cv = LogisticRegressionCV(
    Cs=np.logspace(-4, 2, 20), cv=5, penalty='l1', solver='saga',
    scoring='roc_auc', max_iter=10000, random_state=123, n_jobs=-1
)
lasso_cv_scores = cross_val_score(lasso_for_cv, X_scaled, y, cv=cv_stratified, scoring='roc_auc', n_jobs=-1)
print(f"LASSO 5-fold CV AUC: {lasso_cv_scores.mean():.3f} ± {lasso_cv_scores.std():.3f}")


# --- GAM CV ---
# GAM CV is tricky if feature selection is part of the pipeline to be CV'd.
# For this assignment's scope, we will use the features selected *once* by LASSO on the initial training split.
# This means the feature set for GAM is fixed across folds of *this* CV.
# This is a simplification; a fully rigorous approach would re-run LASSO feature selection in each fold.
if len(gam_feature_indices_for_fit) > 0:
    # We need to re-initialize the GAM model for each fold because it's stateful after fitting.
    # cross_val_score clones the estimator.
    gam_for_cv = LogisticGAM(gam_terms_fit) # Use terms based on initially selected features
    gam_cv_scores = cross_val_score(gam_for_cv, X_scaled[:, gam_feature_indices_for_fit], y, cv=cv_stratified, scoring='roc_auc')
    print(f"GAM ({len(gam_feature_names_for_fit)} LASSO-selected features) 5-fold CV AUC: {gam_cv_scores.mean():.3f} ± {gam_cv_scores.std():.3f}")
else:
    print("GAM not cross-validated as no features were selected by initial LASSO run.")
    gam_cv_scores = np.array([np.nan]) # Placeholder


# --- Gradient Boosting CV ---
# We use the same fixed hyperparameters for GBM as in the initial fit for consistency in this comparison.
# Ideally, GBM hyperparameters (like n_estimators, max_depth, learning_rate) would also be tuned
# within each fold of a nested CV, or a simpler CV to find one best set of HPs.
gb_for_cv = GradientBoostingClassifier(random_state=42, n_estimators=100, max_depth=3, learning_rate=0.1)
gb_cv_scores = cross_val_score(gb_for_cv, X_scaled, y, cv=cv_stratified, scoring='roc_auc', n_jobs=-1)
print(f"Gradient Boosting 5-fold CV AUC: {gb_cv_scores.mean():.3f} ± {gb_cv_scores.std():.3f}")
print("Note for GBM: n_estimators can be tuned, e.g., using gbm.perf in R or validation curves in Python.")


# Plot CV scores for comparison
plt.figure(figsize=(10, 6))
cv_results_df = pd.DataFrame({
    'Model': ['LASSO', 'GAM', 'Gradient Boosting'],
    'Mean AUC': [lasso_cv_scores.mean(), gam_cv_scores.mean(), gb_cv_scores.mean()],
    'Std AUC': [lasso_cv_scores.std(), gam_cv_scores.std(), gb_cv_scores.std()]
}).dropna().sort_values(by='Mean AUC', ascending=False)

sns.barplot(x='Mean AUC', y='Model', data=cv_results_df, xerr=cv_results_df['Std AUC'], palette='viridis')
plt.title('Model Comparison (5-Fold Stratified CV AUC)')
plt.xlabel('Mean Area Under ROC Curve (AUC)')
plt.ylabel('Model')
plt.xlim(0.5, 1.0) # Assuming AUC is between 0.5 (random) and 1.0 (perfect)
for index, row in cv_results_df.iterrows():
    plt.text(row['Mean AUC'] + row['Std AUC'] + 0.01, index, f"{row['Mean AUC']:.3f}", color='black', ha="left", va='center')
plt.show()

# Store the best model from CV for later interpretation if needed
best_model_name = cv_results_df.iloc[0]['Model']
print(f"\\nBest performing model based on mean CV AUC: {best_model_name}")
""")

add_markdown_cell("### 2.4. Interpretation of Optimization Results")
add_markdown_cell("""\
#### 1. Do correlations between variables influence the results? How?

* **Yes, correlations significantly influence results, particularly for LASSO.**
    * **LASSO:** As discussed in lectures (e.g., `Lecture2.pdf`), LASSO's L1 penalty makes it sensitive to multicollinearity. [cite: 22] When predictors are highly correlated, LASSO tends to select one feature from the group and zero out the coefficients of others. This selection can be somewhat arbitrary, meaning a different fold in CV or a slight change in data might lead to a different feature being selected from the correlated cluster. The protein correlation heatmap indicated such correlated clusters. The LASSO coefficient paths plot visually demonstrates how coefficients change (many going to zero) as lambda (inverse of C) increases.
    * **GAM:** If highly correlated features were used as individual spline terms, interpreting their individual contributions would be difficult. By using features pre-selected by LASSO, this issue is somewhat reduced for the GAM step, but the initial LASSO selection itself is affected by correlations.
    * **Gradient Boosting:** Tree-based models like GBM are generally more robust to multicollinearity regarding predictive performance. However, feature importances can be spread across correlated features, potentially diluting the perceived importance of any single feature within such a group.

#### 2. Is over-learning a problem for finding the optimal model?

* **Yes, over-learning (overfitting) is a major concern, especially with $p \\approx n$ (70 features, 70 samples).** [cite: 4, 8]
    * **LASSO:** `LogisticRegressionCV` explicitly combats overfitting by tuning the regularization strength ($C$) through internal cross-validation. The L1 penalty encourages sparsity, leading to simpler models that are less likely to overfit. [cite: 46]
    * **GAM:** GAMs can overfit if splines are overly flexible. Using a reduced feature set (from LASSO) helps. `pygam` uses generalized cross-validation (GCV) by default to select smoothness for splines, which helps prevent overfitting. The 5-fold CV AUC provides an external check on generalization.
    * **Gradient Boosting:** GBMs can overfit with too many trees (`n_estimators`) or too complex trees (`max_depth`). The chosen parameters (`n_estimators=100`, `max_depth=3`) are common defaults but could be tuned. The 5-fold CV AUC for GBM is the primary measure here of its generalization ability.
    * **Overall:** Comparing the mean CV AUC scores and their standard deviations across the models (as plotted) is key. A model with high mean CV AUC and low standard deviation is preferred, indicating good and stable generalization. This aligns with the lecture's emphasis on using CV to estimate test error. [cite: 11, 35]

#### 3. Is there evidence for non-linear effects?

* **GAM:** GAMs are explicitly designed to capture non-linear relationships via splines. If the GAM (even on a reduced feature set) performs comparably to or better than LASSO in the 5-fold CV, it suggests that non-linearities captured by GAM are beneficial.
* **Gradient Boosting:** GBMs inherently model non-linearities. Partial Dependence Plots (PDPs) for top GBM features can visualize these.
* **LASSO:** As a linear model (on the logit scale), LASSO does not capture non-linear relationships in the original feature space. If significant non-linearities exist and are important for prediction, LASSO might underperform.

Let's examine PDPs for the top features from the GBM model fitted on the training split (X_train, y_train) to visualize potential non-linearities.
""")
add_code_cell("""\
# Plotting partial dependence for top features from GAM fitted on X_train
if gam_model_fit is not None and len(gam_feature_names_for_fit) > 0:
    print("\\nPartial Dependence Plots for GAM (top features from training split model):")
    # Determine number of plots (max 3 or number of features)
    n_pdp_gam = min(3, len(gam_feature_names_for_fit))
    fig_gam, axes_gam = plt.subplots(1, n_pdp_gam, figsize=(min(18, 6*n_pdp_gam), 5), squeeze=False) # Ensure axes_gam is 2D
    
    pdp_feature_indices_gam_plot = list(range(n_pdp_gam)) # 0-based for X_test[:, gam_feature_indices_for_fit]
    pdp_feature_names_gam_plot = gam_feature_names_for_fit[:n_pdp_gam]

    PartialDependenceDisplay.from_estimator(
        gam_model_fit,
        X_test[:, gam_feature_indices_for_fit], 
        features=pdp_feature_indices_gam_plot,
        feature_names=pdp_feature_names_gam_plot, 
        ax=axes_gam[0], # Pass the 1D array of axes
        line_kw={"color": "green"}
    )
    fig_gam.suptitle('Partial Dependence Plots for top GAM features (on test data)', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
else:
    print("GAM PDPs not plotted as GAM was not trained or had no features.")

# Plotting partial dependence for top 2 features from GBM fitted on X_train
if len(top10_gb_features) > 0:
    print("\\nPartial Dependence Plots for Gradient Boosting (top 2 features from training split model):")
    gb_pdp_feature_indices_plot = [protein_cols.index(name) for name in top10_gb_features.index[:min(2, len(top10_gb_features))]]
    n_pdp_gb = len(gb_pdp_feature_indices_plot)

    if n_pdp_gb > 0:
        fig_gb, axes_gb = plt.subplots(1, n_pdp_gb, figsize=(6*n_pdp_gb, 5), squeeze=False) # Ensure axes_gb is 2D
        PartialDependenceDisplay.from_estimator(
            gb_model_fit,
            X_train, 
            features=gb_pdp_feature_indices_plot,
            feature_names=protein_cols, 
            ax=axes_gb[0], # Pass the 1D array of axes
            line_kw={"color": "blue"}
        )
        fig_gb.suptitle('Partial Dependence Plots for top GBM features (on training data)', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()
""")
add_markdown_cell("""\
* **Interpretation of PDPs:** If the PDPs show non-flat, non-linear curves, it provides evidence that the model (GBM or GAM) has learned and is leveraging non-linear relationships for those features. The GAM summary (if `gam_model_fit.summary()` was called) would also indicate the effective degrees of freedom for each spline, where values > 1 suggest non-linearity.

#### 4. Is there evidence for important interactions between variables?

* **Gradient Boosting:** GBMs are adept at capturing interactions between features since decision trees naturally model them. If GBM significantly outperforms LASSO and the main-effects GAM (which doesn't explicitly model interactions unless terms like `te(X1,X2)` are added), this could suggest that feature interactions are important for predicting Treatment.
* **LASSO & Main-Effects GAM:** These models, as implemented here, primarily capture main effects. They do not inherently model interactions unless interaction terms are manually engineered and included.

The relative CV AUC scores will be the primary guide. If GBM is the top performer by a notable margin, it's plausible that its ability to model both non-linearities and interactions contributes to this success.
""")

# --- Question 3: Memantine effects in trisomic, stimulated mice ---
add_markdown_cell("""\
## 3. Evaluate Memantine's Effect on Protein Levels in Trisomic, Stimulated Mice

This section investigates whether Memantine treatment induces specific changes in protein expression levels in the subgroup of trisomic mice (Genotype='Ts65Dn') that were stimulated to learn (Behavior='C/S'). This involves statistical testing for differences in mean protein expression between Memantine and Saline groups within this specific cohort, with appropriate correction for multiple comparisons.
""")

add_code_cell("""\
# Purpose: Identify proteins differentially expressed due to Memantine treatment
# in a specific subgroup of mice. This involves:
# 1. Subsetting the data.
# 2. Performing independent t-tests for each protein.
# 3. Correcting for multiple hypothesis testing using FDR (Benjamini-Hochberg).

# Filter the DataFrame for the relevant subgroup
subset_df = df[(df['Genotype'] == 'Ts65Dn') & (df['Behavior'] == 'C/S')].copy()

print(f"Number of samples in the subset (Ts65Dn, C/S): {subset_df.shape[0]}")
print("\\nTreatment distribution in this subset:")
print(subset_df['Treatment'].value_counts())

# Perform t-tests for each protein
p_values_q3 = []
protein_names_q3 = []
mean_differences_q3 = []

mem_group_q3 = subset_df[subset_df['Treatment'] == 'Memantine'][protein_cols]
sal_group_q3 = subset_df[subset_df['Treatment'] == 'Saline'][protein_cols]

if mem_group_q3.shape[0] < 2 or sal_group_q3.shape[0] < 2: # Need at least 2 samples per group for t-test
    print("\\nNot enough samples in one or both treatment groups in the subset for robust t-tests.")
else:
    for protein in protein_cols:
        mem_data = mem_group_q3[protein].dropna()
        sal_data = sal_group_q3[protein].dropna()
        
        # Ensure sufficient data points in each group for the specific protein
        if len(mem_data) >= 2 and len(sal_data) >= 2:
            t_stat, p_val = ttest_ind(mem_data, sal_data, nan_policy='omit', equal_var=False) # Welch's t-test (unequal variance)
            p_values_q3.append(p_val)
            protein_names_q3.append(protein)
            mean_differences_q3.append(mem_data.mean() - sal_data.mean())
        else:
            # print(f"Skipping t-test for {protein} in Q3 due to insufficient data.")
            pass

    if not protein_names_q3:
        print("No proteins were eligible for t-tests in Q3 due to insufficient data in subgroups after dropping NaNs.")
    else:
        # Apply False Discovery Rate (FDR) correction
        reject_q3, pvals_corrected_q3, _, _ = smm.multipletests(p_values_q3, alpha=0.05, method='fdr_bh')
        
        results_q3_df = pd.DataFrame({
            'Protein': protein_names_q3,
            'Mean_Difference (Mem-Sal)': mean_differences_q3,
            'P_value': p_values_q3,
            'FDR_Corrected_P_value': pvals_corrected_q3,
            'Significant_FDR_0.05': reject_q3
        })
        
        significant_proteins_q3_df = results_q3_df[results_q3_df['Significant_FDR_0.05']]
        print(f"\\nNumber of proteins tested in Q3: {len(protein_names_q3)}")
        print(f"Number of significant proteins in Q3 (FDR < 0.05): {significant_proteins_q3_df.shape[0]}")
        
        if not significant_proteins_q3_df.empty:
            print("\\nSignificant proteins in Q3 (Ts65Dn, C/S; Memantine vs Saline):")
            print(significant_proteins_q3_df.sort_values(by='FDR_Corrected_P_value'))
            
            # Plot mean differences for significant proteins
            plt.figure(figsize=(10, max(6, 0.5 * len(significant_proteins_q3_df))))
            sns.barplot(x='Mean_Difference (Mem-Sal)', y='Protein', 
                        data=significant_proteins_q3_df.sort_values('Mean_Difference (Mem-Sal)'), 
                        palette='coolwarm_r')
            plt.title('Proteins Significantly Affected by Memantine in Ts65Dn/C/S Mice (FDR < 0.05)')
            plt.xlabel('Mean Expression Difference (Memantine - Saline)')
            plt.ylabel('Protein')
            plt.axvline(0, color='grey', linestyle='--')
            plt.show()
        else:
            print("No proteins showed statistically significant changes after FDR correction in this subgroup.")

        # Volcano plot for Q3 results
        # Ensure pvals_corrected_q3 is not empty before log10
        if pvals_corrected_q3.size > 0:
            log10_fdr_pvals_q3 = -np.log10(pvals_corrected_q3)
            # Handle potential -inf if p_corrected is 1, or inf if p_corrected is 0
            log10_fdr_pvals_q3[np.isinf(log10_fdr_pvals_q3) & (pvals_corrected_q3 == 0)] = np.nanmax(log10_fdr_pvals_q3[np.isfinite(log10_fdr_pvals_q3)]) + 1 # Replace inf with a large value
            log10_fdr_pvals_q3[np.isinf(log10_fdr_pvals_q3) & (pvals_corrected_q3 == 1)] = 0 # Replace -inf with 0


            volcano_q3_df = pd.DataFrame({
                'Protein': protein_names_q3,
                'Mean_Difference': mean_differences_q3,
                'neg_log10_FDR_Pval': log10_fdr_pvals_q3,
                'Significant': reject_q3
            })

            plt.figure(figsize=(10, 7))
            sns.scatterplot(
                x='Mean_Difference', 
                y='neg_log10_FDR_Pval', 
                data=volcano_q3_df, 
                hue='Significant', 
                palette={True: 'red', False: 'grey'},
                alpha=0.7
            )
            plt.axhline(-np.log10(0.05), color='blue', linestyle='--', label='FDR = 0.05 Threshold')
            # Optional: add lines for effect size thresholds if meaningful
            # plt.axvline(1, color='orange', linestyle='--', label='Effect Size = +1')
            # plt.axvline(-1, color='orange', linestyle='--')
            plt.title('Volcano Plot: Memantine vs. Saline in Ts65Dn, C/S Mice')
            plt.xlabel('Mean Difference (Memantine - Saline)')
            plt.ylabel('-log10(FDR Corrected P-value)')
            if volcano_q3_df['Significant'].any():
                 plt.legend(title='Significant (FDR < 0.05)')
            else:
                 plt.legend().remove()
            plt.show()
        else:
            print("Volcano plot not generated for Q3 as no proteins were tested or p-values were problematic.")
""")

add_markdown_cell("""\
### 3.1. Conclusion for Memantine's Effect

The analysis of the specific subgroup (trisomic 'Ts65Dn' mice stimulated to learn 'C/S') involved comparing protein expression between Memantine-treated (n=9) and Saline-treated (n=7) mice. Independent two-sample t-tests were performed for each protein, followed by Benjamini-Hochberg FDR correction to account for multiple comparisons.

**Results:**
* **(The Python code output above will state whether any proteins were found to be significantly different after FDR correction at an alpha of 0.05.)**

**Interpretation:**
* If significant proteins are listed, these are the protein levels specifically altered by Memantine treatment in this context of Down syndrome model mice undergoing a learning task. The direction of change (increase or decrease with Memantine) would be indicated by the sign of the mean difference.
* If no proteins are found to be significantly different (as was the case in the original notebook output), it suggests that within the statistical power afforded by this small subgroup size, Memantine does not induce detectable changes in the measured protein levels that pass the FDR < 0.05 threshold. The volcano plot would visually confirm this.

**Limitations:**
The primary limitation of this subgroup analysis is the very small sample size (9 Memantine vs. 7 Saline). This significantly reduces the statistical power to detect subtle or even moderate changes in protein expression. Therefore, a null result (no significant proteins) should be interpreted with caution and does not necessarily mean Memantine has no effect, but rather that any effect is not strong enough to be detected with high confidence in this small sample.
""")

# --- Final Conclusion (Overall Summary) ---
add_markdown_cell("""\
## 4. Overall Summary and Conclusions

This notebook performed an analysis of protein expression data from mice to investigate data characteristics, compare predictive models for treatment classification, and assess specific drug effects in a targeted subgroup.

**1. Data Description:**
The dataset (70 samples, 70 protein features) presented characteristics of a high-dimensional biomedical dataset ($p \\approx n$). Key observations included:
* No missing values.
* Significant multicollinearity among proteins.
* Moderate sample sizes for overall classification, but very small sample sizes for specific subgroup analyses (Q3).
These factors highlight risks like overfitting and the need for careful model validation (e.g., using CV and regularization). [cite: 4, 8, 11, 35, 63]

**2. Model Comparison (Treatment Classification: Memantine vs. Saline):**
LASSO, GAM, and Gradient Boosting (GBM) models were compared using 5-fold stratified cross-validation AUC on the full scaled dataset.
* The CV AUC scores were approximately:
    * LASSO: *(Refer to actual output from cell 2.3)*
    * GAM (on LASSO-selected features): *(Refer to actual output from cell 2.3)*
    * Gradient Boosting: *(Refer to actual output from cell 2.3)*
* The best performing model based on mean CV AUC was *(State the best model from output)*.
* **Correlations:** LASSO's feature selection is sensitive to multicollinearity. [cite: 22] The coefficient paths plot showed how features are selected/shrunk.
* **Over-learning:** Addressed via LASSO's regularization[cite: 46], GAM's internal smoothness estimation, and robust CV evaluation for all models.
* **Non-linear effects:** PDPs for GAM and GBM provided visual evidence. The relative performance of these models against LASSO further suggested *(comment on importance based on results)*.
* **Interactions:** GBM's potential to capture interactions might explain its performance if it was superior.

**3. Memantine Effect in Ts65Dn, C/S Mice:**
T-tests with FDR correction were used to identify proteins differentially expressed between Memantine and Saline treatments in trisomic mice stimulated to learn.
* *(Report specific findings from Q3 output: e.g., "No proteins showed statistically significant changes..." or list significant ones)*.
* This finding is significantly limited by the small sample size (n=9 Memantine, n=7 Saline) of this subgroup.

**Recommendations for Future Work:**
* More extensive hyperparameter tuning for GAM and GBM using nested cross-validation.
* Exploration of advanced feature selection methods beyond LASSO, or ensemble feature selection.
* Techniques to explicitly model or analyze feature interactions in GBM (e.g., SHAP values).
* If pursuing subgroup analyses like Q3, larger sample sizes would be crucial for more robust conclusions.
""")

# Save the notebook to a file
notebook_filename = "revised_assignment_notebook_v2.ipynb"
with open(notebook_filename, 'w') as f:
    nbformat.write(nb, f)

print(f"\\nRevised notebook saved as {notebook_filename}")