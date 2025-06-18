import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.linear_model import LogisticRegressionCV
from pygam import LogisticGAM, s
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.inspection import PartialDependenceDisplay
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smm
from sklearn.decomposition import PCA

def analyze_protein_expression(file_path="Assignment/Cortex3.csv"):
    """
    Analyzes protein expression data in control and Down Syndrome mice.
    Performs data exploration, trains and compares classification models,
    and evaluates specific treatment effects.

    Args:
        file_path (str): The path to the CSV dataset.
    """

    # Set plotting style
    sns.set(style='whitegrid')
    plt.rcParams['figure.figsize'] = (10, 6)

    # 0. Setup and Data Loading
    print("--- 0. Setup and Data Loading ---")
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: {file_path} not found. Please ensure the file is in the correct directory.")
        raise

    print("Data Shape (samples, columns):", df.shape)
    print("\nFirst 5 rows of the data:")
    print(df.head())

    # Identify protein columns (ending with '_N') and categorical/metadata columns
    protein_cols = [col for col in df.columns if col.endswith('_N')]
    metadata_cols = ['Genotype', 'Treatment', 'Behavior']
    print(f"\nFound {len(protein_cols)} protein expression columns.")
    print(f"Found {len(metadata_cols)} metadata columns: {metadata_cols}")

    # 1. Study and Describe the Data
    print("\n--- 1. Study and Describe the Data ---")

    # 1.1. Initial Data Inspection
    print("\n--- 1.1. Initial Data Inspection ---")
    # Check for missing values
    print("\nMissing values per column (showing columns with any missing values):")
    missing_values = df.isnull().sum()
    print(missing_values[missing_values > 0])
    if missing_values.sum() == 0:
        print("No missing values found in the dataset.")

    # Check data types
    print("\nData types summary:")
    print(df.dtypes.value_counts())
    print("\nData types of metadata columns:")
    print(df[metadata_cols].dtypes)
    print("\nData types of first 3 protein columns:")
    print(df[protein_cols[:3]].dtypes)

    # Convert categorical columns from 'object' to 'category' type
    for col in metadata_cols:
        if df[col].dtype == 'object':
            df[col] = df[col].astype('category')
    print("\nData types of metadata columns after potential conversion:")
    print(df[metadata_cols].dtypes)

    # 1.2. Descriptive Statistics
    print("\n--- 1.2. Descriptive Statistics ---")
    # Summary statistics for numerical (protein) columns
    print("\nSummary statistics for protein expression levels (first 5 proteins):")
    print(df[protein_cols].describe().transpose().head())

    # Value counts for categorical columns
    print("\nDistribution of samples across Genotype:")
    print(df['Genotype'].value_counts())
    print("\nDistribution of samples across Treatment:")
    print(df['Treatment'].value_counts())
    print("\nDistribution of samples across Behavior:")
    print(df['Behavior'].value_counts())

    # 1.3. Visual Data Exploration
    print("\n--- 1.3. Visual Data Exploration ---")

    # 1.3.1. Distribution of Protein Expression Levels
    print("\nHistograms for the first 3 protein expression levels:")
    plt.figure(figsize=(18, 5))
    for i, col in enumerate(protein_cols[:3]):
        plt.subplot(1, 3, i + 1)
        sns.histplot(df[col], kde=True, bins=15)
        plt.title(f'Distribution of {col}')
    plt.tight_layout()
    plt.show()

    print("\nBoxplots for the first 3 protein expression levels by Treatment:")
    plt.figure(figsize=(18, 5))
    for i, col in enumerate(protein_cols[:3]):
        plt.subplot(1, 3, i + 1)
        sns.boxplot(x='Treatment', y=col, data=df)
        plt.title(f'{col} by Treatment')
    plt.tight_layout()
    plt.show()

    # 1.3.2. Protein Expression Correlation Heatmap
    corr_matrix = df[protein_cols].corr()
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, cmap="coolwarm", center=0, annot=False, xticklabels=False, yticklabels=False)
    plt.title("Protein Expression Correlation Heatmap")
    plt.show()

    # 1.3.3. Principal Component Analysis (PCA)
    # Scale protein expression data before PCA
    X_proteins_scaled = StandardScaler().fit_transform(df[protein_cols])

    # Perform PCA
    pca_full = PCA(n_components=min(10, X_proteins_scaled.shape[1]))
    pca_full.fit(X_proteins_scaled)

    pca_2comp = PCA(n_components=2)
    pcs_2comp = pca_2comp.fit_transform(X_proteins_scaled)
    pc_df = pd.DataFrame(data=pcs_2comp, columns=['PC1', 'PC2'])
    pc_df_combined = pd.concat([pc_df, df[metadata_cols].reset_index(drop=True)], axis=1)

    print(f"\nExplained variance by PC1: {pca_2comp.explained_variance_ratio_[0]:.3f}")
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

    # 1.4. Indications of Potential Issues for Statistical Modeling
    print("\n--- 1.4. Indications of Potential Issues for Statistical Modeling ---")
    print("1. High Dimensionality (p ~ n): 70 samples, 70 features increases overfitting risk.")
    print("2. Multicollinearity: Correlation heatmap shows correlated protein blocks, affecting LASSO coefficient stability.")
    print("3. Small Sample Size per Class/Subgroup: Small numbers in specific subgroups (e.g., Ts65Dn, C/S, Memantine) reduce statistical power.")
    print("4. Potential Outliers and Distributions: Some skewed protein distributions may affect model assumptions or influence.")
    print("5. Model Selection Complexity: Careful cross-validation is essential to avoid overfitting and ensure generalizability.")


    # 2. Train and Compare LASSO, GAM, and Boosting Models
    print("\n--- 2. Train and Compare LASSO, GAM, and Boosting Models ---")

    # 2.1. Data Preparation for Classification
    print("\n--- 2.1. Data Preparation for Classification ---")
    # Define feature matrix X (protein expressions) and target y (Treatment)
    X = df[protein_cols].values
    y = (df['Treatment'] == 'Memantine').astype(int)

    label_encoder = LabelEncoder()
    label_encoder.classes_ = np.array(['Saline', 'Memantine'])
    print(f"Original Treatment labels: {df['Treatment'].unique()}")
    print(f"Encoded Treatment labels: {np.unique(y)} maps to classes: {label_encoder.classes_}")

    # Scale the features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.3, stratify=y, random_state=42)
    print(f"\nX_train shape: {X_train.shape}, y_train shape: {y_train.shape}")
    print(f"X_test shape: {X_test.shape}, y_test shape: {y_test.shape}")
    print(f"Class distribution in y_train: Saline (0)={np.sum(y_train==0)}, Memantine (1)={np.sum(y_train==1)}")
    print(f"Class distribution in y_test: Saline (0)={np.sum(y_test==0)}, Memantine (1)={np.sum(y_test==1)}")

    # 2.2. Model Training and Evaluation on Test Split
    print("\n--- 2.2. Model Training and Evaluation on Test Split ---")

    # 2.2.1. LASSO (Logistic Regression with L1 penalty)
    print("\n--- 2.2.1. LASSO (Logistic Regression with L1 penalty) ---")
    lasso_model_cv = LogisticRegressionCV(
        Cs=np.logspace(-4, 2, 20),
        cv=5,
        penalty='l1',
        solver='saga',
        scoring='roc_auc',
        max_iter=10000,
        random_state=42,
        n_jobs=-1
    )
    lasso_model_cv.fit(X_train, y_train)

    print(f"Best C for LASSO (from training data CV): {lasso_model_cv.C_[0]:.4f}")

    # Evaluate on the held-out test set
    y_prob_lasso_test = lasso_model_cv.predict_proba(X_test)[:, 1]
    auc_lasso_test = roc_auc_score(y_test, y_prob_lasso_test)
    print(f"LASSO Test AUC: {auc_lasso_test:.3f}")

    # Inspect coefficients from the model trained with the best C
    lasso_coeffs = pd.Series(lasso_model_cv.coef_[0], index=protein_cols)
    selected_features_lasso = lasso_coeffs[lasso_coeffs != 0]
    print(f"\nNumber of features selected by LASSO: {len(selected_features_lasso)}")
    if len(selected_features_lasso) > 0:
        print("Top 5 selected features by LASSO (absolute coefficient value):")
        print(selected_features_lasso.abs().sort_values(ascending=False).head())
    else:
        print("No features selected by LASSO with the chosen C values.")

    # Plot coefficient paths for LogisticRegressionCV
    coefs_path = lasso_model_cv.coefs_paths_[1]
    if coefs_path.ndim == 3:
        coefs_path_mean = np.mean(coefs_path, axis=0)
    else:
        coefs_path_mean = coefs_path

    plt.figure(figsize=(12, 7))
    lambdas_lasso = 1. / lasso_model_cv.Cs_
    sorted_lambda_indices = np.argsort(lambdas_lasso)

    for i in range(coefs_path_mean.shape[0]):
        plt.plot(lambdas_lasso[sorted_lambda_indices], coefs_path_mean[i, sorted_lambda_indices], label=f'{protein_cols[i]}' if i < 5 else None)
    plt.xscale('log')
    plt.xlabel('Lambda (Inverse of C, Log Scale)')
    plt.ylabel('Coefficient Value')
    plt.title('LASSO Coefficient Paths')
    plt.axvline(1./lasso_model_cv.C_[0], color='red', linestyle='--', label=f'Optimal Lambda (1/C_best = {1./lasso_model_cv.C_[0]:.4f})')
    plt.legend()
    plt.show()

    # 2.2.2. Generalized Additive Model (GAM)
    print("\n--- 2.2.2. Generalized Additive Model (GAM) ---")
    if len(selected_features_lasso) == 0:
        print("LASSO selected 0 features. GAM training will be skipped.")
        gam_feature_indices_for_fit = []
        gam_feature_names_for_fit = []
        gam_model_fit = None
        auc_gam_test = np.nan
    else:
        gam_feature_names_for_fit = selected_features_lasso.index.tolist()
        gam_feature_indices_for_fit = [protein_cols.index(name) for name in gam_feature_names_for_fit]

        print(f"\nGAM will be trained on {len(gam_feature_names_for_fit)} features selected by LASSO: {gam_feature_names_for_fit[:5]}...")

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
        # gam_model_fit.summary() # Uncomment to print GAM summary

    # 2.2.3. Gradient Boosting Machine (GBM)
    print("\n--- 2.2.3. Gradient Boosting Machine (GBM) ---")
    gb_model_fit = GradientBoostingClassifier(random_state=42, n_estimators=100, max_depth=3, learning_rate=0.1)
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

    # 2.3. Model Comparison using 5-Fold Cross-Validation AUC
    print("\n--- 2.3. Model Comparison using 5-Fold Cross-Validation AUC ---")
    cv_stratified = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # LASSO CV
    lasso_for_cv = LogisticRegressionCV(
        Cs=np.logspace(-4, 2, 20), cv=5, penalty='l1', solver='saga',
        scoring='roc_auc', max_iter=10000, random_state=123, n_jobs=-1
    )
    lasso_cv_scores = cross_val_score(lasso_for_cv, X_scaled, y, cv=cv_stratified, scoring='roc_auc', n_jobs=-1)
    print(f"LASSO 5-fold CV AUC: {lasso_cv_scores.mean():.3f} \u00B1 {lasso_cv_scores.std():.3f}")

    # GAM CV
    if len(gam_feature_indices_for_fit) > 0:
        gam_for_cv = LogisticGAM(gam_terms_fit)
        gam_cv_scores = cross_val_score(gam_for_cv, X_scaled[:, gam_feature_indices_for_fit], y, cv=cv_stratified, scoring='roc_auc')
        print(f"GAM ({len(gam_feature_names_for_fit)} LASSO-selected features) 5-fold CV AUC: {gam_cv_scores.mean():.3f} \u00B1 {gam_cv_scores.std():.3f}")
    else:
        print("GAM not cross-validated as no features were selected by initial LASSO run.")
        gam_cv_scores = np.array([np.nan])

    # Gradient Boosting CV
    gb_for_cv = GradientBoostingClassifier(random_state=42, n_estimators=100, max_depth=3, learning_rate=0.1)
    gb_cv_scores = cross_val_score(gb_for_cv, X_scaled, y, cv=cv_stratified, scoring='roc_auc', n_jobs=-1)
    print(f"Gradient Boosting 5-fold CV AUC: {gb_cv_scores.mean():.3f} \u00B1 {gb_cv_scores.std():.3f}")
    print("Note for GBM: n_estimators can be tuned, e.g., using gbm.perf in R or validation curves in Python.")

    # Plot CV scores for comparison
    plt.figure(figsize=(10, 6))
    cv_results_df = pd.DataFrame({
        'Model': ['LASSO', 'GAM', 'Gradient Boosting'],
        'Mean AUC': [lasso_cv_scores.mean(), gam_cv_scores.mean(), gb_cv_scores.mean()],
        'Std AUC': [lasso_cv_scores.std(), gam_cv_scores.std(), gb_cv_scores.std()]
    }).dropna().sort_values(by='Mean AUC', ascending=False)

    xerr = cv_results_df['Std AUC'].values
    sns.barplot(x='Mean AUC', y='Model', data=cv_results_df, xerr=xerr, palette='viridis')
    plt.title('Model Comparison (5-Fold Stratified CV AUC)')
    plt.xlabel('Mean Area Under ROC Curve (AUC)')
    plt.ylabel('Model')
    plt.xlim(0.5, 1.0)

    for idx, (mean_auc, std_auc) in enumerate(zip(cv_results_df['Mean AUC'], cv_results_df['Std AUC'])):
        plt.text(mean_auc + std_auc + 0.01, idx, f"{mean_auc:.3f}", color='black', ha="left", va='center')
    plt.show()

    best_model_name = cv_results_df.iloc[0]['Model']
    print(f"\nBest performing model based on mean CV AUC: {best_model_name}")

    # 2.4. Interpretation of Optimization Results
    print("\n--- 2.4. Interpretation of Optimization Results ---")
    print("1. Influence of correlations between variables:")
    print("   - LASSO: Highly sensitive to multicollinearity; tends to select one feature from a correlated group and zero out others. This selection can be arbitrary.")
    print("   - GAM: If highly correlated features are used, interpreting individual spline contributions can be difficult, although pre-selection by LASSO helps.")
    print("   - Gradient Boosting: Generally robust to multicollinearity for predictive accuracy, but feature importance can be distributed across correlated features.")

    print("\n2. Over-learning as a problem for finding the optimal model:")
    print("   - Overfitting is a significant concern given p \u2248 n (70 features, 70 samples).")
    print("   - LASSO: Combats overfitting through L1 regularization and internal cross-validation for C selection, promoting sparser models.")
    print("   - GAM: Can overfit if splines are too flexible. Using a reduced feature set and GCV for smoothness helps. 5-fold CV provides external check.")
    print("   - Gradient Boosting: Prone to overfitting with too many or too complex trees. Cross-validation for hyperparameter tuning (e.g., n_estimators, max_depth) is crucial.")
    print("   - Overall: Comparing mean CV AUC and standard deviations is key to identifying models with good and stable generalization ability.")

    print("\n3. Evidence for non-linear effects:")
    print("   - GAM and Gradient Boosting are designed to capture non-linear relationships. If their performance (CV AUC) is comparable or superior to LASSO, it suggests non-linearities are beneficial.")

    # Plotting partial dependence (spline plots) for top features from GAM
    if gam_model_fit is not None and len(gam_feature_names_for_fit) > 0:
        print("\nPartial Dependence (Spline) Plots for GAM (selected features):")
        n_plot = min(3, len(gam_feature_names_for_fit))
        fig, axes = plt.subplots(1, n_plot, figsize=(min(18, 5 * n_plot), 5))
        if n_plot == 1:
            axes = [axes]

        for i in range(n_plot):
            XX = gam_model_fit.generate_X_grid(term=i)
            pdep = gam_model_fit.partial_dependence(term=i, X=XX)
            axes[i].plot(XX[:, 0], pdep, color='b')
            axes[i].set_title(f"Spline for {gam_feature_names_for_fit[i]}")
            axes[i].set_xlabel(gam_feature_names_for_fit[i])
            axes[i].set_ylabel("Partial Dependence")
        fig.suptitle('GAM Spline Plots for Top Features', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()
    else:
        print("GAM spline plots not shown as no features were selected or GAM not fitted.")

    # Plotting partial dependence for top 2 features from GBM
    if len(top10_gb_features) > 0:
        print("\nPartial Dependence Plots for Gradient Boosting (top 2 features from training split model):")
        gb_pdp_feature_indices_plot = [protein_cols.index(name) for name in top10_gb_features.index[:min(2, len(top10_gb_features))]]
        n_pdp_gb = len(gb_pdp_feature_indices_plot)
        if n_pdp_gb > 0:
            fig_gb, axes_gb = plt.subplots(1, n_pdp_gb, figsize=(6 * n_pdp_gb, 5), squeeze=True)
            if n_pdp_gb == 1:
                axes_gb = [axes_gb]
            PartialDependenceDisplay.from_estimator(
                gb_model_fit,
                X_train,
                features=gb_pdp_feature_indices_plot,
                feature_names=protein_cols,
                ax=axes_gb,
                line_kw={"color": "blue"}
            )
            fig_gb.suptitle('Partial Dependence Plots for top GBM features (on training data)', fontsize=16)
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            plt.show()
    else:
        print("GBM partial dependence plots not shown as no top features were identified.")

    print("\n4. Evidence for important interactions between variables:")
    print("   - Gradient Boosting: GBMs naturally capture feature interactions. If GBM significantly outperforms LASSO and main-effects GAM, it suggests interactions are important.")
    print("   - LASSO & Main-Effects GAM: These models primarily capture main effects unless interaction terms are explicitly engineered.")
    print("   - The relative CV AUC scores are the primary indicator; if GBM is notably superior, it points to the importance of non-linearities and/or interactions.")


    # 3. Evaluate Memantine's Effect on Protein Levels in Trisomic, Stimulated Mice
    print("\n--- 3. Evaluate Memantine's Effect on Protein Levels in Trisomic, Stimulated Mice ---")
    # Filter the DataFrame for the relevant subgroup
    subset_df = df[(df['Genotype'] == 'Ts65Dn') & (df['Behavior'] == 'C/S')].copy()
    print(f"Number of samples in the subset (Ts65Dn, C/S): {subset_df.shape[0]}")
    print("\nTreatment distribution in this subset:")
    print(subset_df['Treatment'].value_counts())

    # Perform t-tests for each protein
    p_values_q3 = []
    protein_names_q3 = []
    mean_differences_q3 = []

    mem_group_q3 = subset_df[subset_df['Treatment'] == 'Memantine'][protein_cols]
    sal_group_q3 = subset_df[subset_df['Treatment'] == 'Saline'][protein_cols]

    if mem_group_q3.shape[0] < 2 or sal_group_q3.shape[0] < 2:
        print("\nNot enough samples in one or both treatment groups in the subset for robust t-tests.")
    else:
        for protein in protein_cols:
            mem_data = mem_group_q3[protein].dropna()
            sal_data = sal_group_q3[protein].dropna()

            if len(mem_data) >= 2 and len(sal_data) >= 2:
                t_stat, p_val = ttest_ind(mem_data, sal_data, nan_policy='omit', equal_var=False) # Welch's t-test
                p_values_q3.append(p_val)
                protein_names_q3.append(protein)
                mean_differences_q3.append(mem_data.mean() - sal_data.mean())
            # else: print(f"Skipping t-test for {protein} in Q3 due to insufficient data after dropping NaNs.")


    if not protein_names_q3:
        print("No proteins were eligible for t-tests in Q3 due to insufficient data in subgroups after dropping NaNs.")
    else:
        # Apply False Discovery Rate (FDR) correction
        reject_q3, pvals_corrected_q3, _, _ = smm.multipletests(p_values_q3, alpha=0.05, method='fdr_bh')

        results_q3_df = pd.DataFrame({
            'Protein': protein_names_q3,
            'Mean_Difference (Memantine - Saline)': mean_differences_q3,
            'P_Value': p_values_q3,
            'FDR_Corrected_P_Value': pvals_corrected_q3,
            'Significant_FDR_0.05': reject_q3
        })
        results_q3_df = results_q3_df.sort_values(by='FDR_Corrected_P_Value')
        print("\nResults of t-tests for protein expression differences in Ts65Dn, C/S mice:")
        print(results_q3_df.head(10))

        significant_proteins_q3 = results_q3_df[results_q3_df['Significant_FDR_0.05']]
        if not significant_proteins_q3.empty:
            print(f"\nNumber of proteins significantly different (FDR < 0.05): {len(significant_proteins_q3)}")
            print("Significantly different proteins:")
            print(significant_proteins_q3)

            # Optional: Visualize significant proteins
            print("\nVisualizing top 3 significantly different proteins (if any):")
            top_significant_proteins = significant_proteins_q3.head(3)['Protein'].tolist()
            if top_significant_proteins:
                plt.figure(figsize=(15, 5))
                for i, protein in enumerate(top_significant_proteins):
                    plt.subplot(1, len(top_significant_proteins), i + 1)
                    sns.boxplot(x='Treatment', y=protein, data=subset_df)
                    plt.title(f'{protein} Expression by Treatment\n(Ts65Dn, C/S Subset)')
                    plt.ylabel('Expression Level')
                plt.tight_layout()
                plt.show()
            else:
                print("No significant proteins to visualize.")
        else:
            print("\nNo proteins showed a statistically significant difference in expression between Memantine and Saline treated mice (FDR < 0.05) in the Ts65Dn, C/S subset.")

if __name__ == "__main__":
    analyze_protein_expression(file_path="Cortex3.csv")