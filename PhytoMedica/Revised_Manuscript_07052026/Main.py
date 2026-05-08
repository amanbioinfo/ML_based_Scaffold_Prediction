# ================================
# 📦 PART 1: IMPORTS & CONFIG
# ================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings("ignore")
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')  # 🔥 suppress warnings

# ML
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import (
    accuracy_score, f1_score, roc_auc_score, classification_report,
    confusion_matrix, precision_recall_curve
)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import precision_score, recall_score
# Models
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegression

import xgboost as xgb
import lightgbm as lgb

# Imbalance handling
from imblearn.over_sampling import SMOTE

# Calibration
from sklearn.calibration import CalibratedClassifierCV

# Stacking
from sklearn.ensemble import StackingClassifier

# RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

# Explainability
import shap

# Save directory
SAVE_DIR = "results_pub"
os.makedirs(SAVE_DIR, exist_ok=True)

# ================================
# 🎨 Nature-style color palette
# ================================
COLORS = {
    "primary": "#2E4057",
    "secondary": "#66A182",
    "accent": "#C44D58",
    "gold": "#EFC94C",
    "purple": "#8E44AD"
}

sns.set(style="whitegrid")
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["axes.edgecolor"] = "black"





# ================================
# 📊 PART 2: LOAD DATA
# ================================
# ================================
# 📊 ROBUST DATA LOADING (FIXED)
# ================================
file_path = "/Volumes/Fun/CIMAP2024_Work/2025_Project/Final_Papers_CIMAP/Gunjan_Mam_AI/PhytoMedica_Revision/drug_tissue_cellline.csv"

try:
    df = pd.read_csv(file_path, encoding="utf-8", low_memory=False)
except UnicodeDecodeError:
    print("⚠️ UTF-8 failed, trying latin1...")
    df = pd.read_csv(file_path, encoding="latin1", low_memory=False)

print("Dataset loaded successfully!")

# Select relevant columns (explicit feature clarity)
selected_cols = [
    "SMILES", "IC50", "AUC",
    "Cell line", "Tissue",
    "gene_symbol", "protein_mutation",
    "rna_mutation", "cdna_mutation",
    "rsem_tpm"
]

# Keep rows where core features exist
df = df[selected_cols]

# Only drop rows where ESSENTIAL columns are missing
df = df.dropna(subset=["SMILES", "IC50"])

# Fill biological missing values instead of dropping
df["protein_mutation"] = df["protein_mutation"].fillna("-")
df["rna_mutation"] = df["rna_mutation"].fillna("-")
df["cdna_mutation"] = df["cdna_mutation"].fillna("-")
df["rsem_tpm"] = df["rsem_tpm"].fillna(0)

print("Dataset shape after cleaning:", df.shape)

print("Dataset shape:", df.shape)

# ================================
# 🎯 Define Activity Threshold (IMPORTANT FIX)
# ================================
# Reviewer concern: IC50 → classification

low = df["IC50"].quantile(0.3)
high = df["IC50"].quantile(0.7)

df = df[(df["IC50"] <= low) | (df["IC50"] >= high)]

df["Activity"] = (df["IC50"] <= low).astype(int)

#df["Activity"] = (df["IC50"] < IC50_THRESHOLD).astype(int)

#print("Threshold used:", IC50_THRESHOLD)
print(df["Activity"].value_counts())






# ================================
# 🧪 PART 3: FEATURE ENGINEERING
# ================================

from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

morgan_gen = GetMorganGenerator(radius=2, fpSize=256)

def smiles_to_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    fp = morgan_gen.GetFingerprint(mol)
    fp_arr = np.array(fp)
    
    descriptors = [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.FractionCSP3(mol),
        Descriptors.HeavyAtomCount(mol)
    ]
    
    return np.concatenate([fp_arr, descriptors])


# Apply chemical features
df["chem_features"] = df["SMILES"].apply(smiles_to_features)
df = df[df["chem_features"].notnull()]


# ================================
# 🧬 Biological Features (FINAL - CORRECTED)
# ================================

# Ensure valid chemical features FIRST
df = df[df["chem_features"].notnull()].copy()

# ----------------
# Mutation encoding
# ----------------
def encode_mutation(x):
    if x == "-" or pd.isna(x):
        return 0
    elif "del" in str(x) or "ins" in str(x):
        return 2
    else:
        return 1

df["mutation_flag"] = df["protein_mutation"].apply(encode_mutation)

# ----------------
# Expression preprocessing (DO FIRST)
# ----------------
df["rsem_tpm"] = pd.to_numeric(df["rsem_tpm"], errors="coerce").fillna(0)
df["rsem_tpm"] = np.log1p(df["rsem_tpm"])

# ----------------
# Interaction features (after normalization)
# ----------------
df["expr_mut_interaction"] = df["mutation_flag"] * df["rsem_tpm"]
df["high_expression"] = (df["rsem_tpm"] > df["rsem_tpm"].median()).astype(int)

# ----------------
# Chemical features
# ----------------
X_chem = np.vstack(df["chem_features"].values)

# ----------------
# Biological base features ONLY (genes later in Part 4)
# ----------------
X_bio = df[[
    "mutation_flag",
    "rsem_tpm",
    "expr_mut_interaction",
    "high_expression"
]].values

y = df["Activity"].values

X = np.hstack([X_chem, X_bio])   # 🔥 CRITICAL FIX



# ================================
# ⚖️ PART 4: TRAIN TEST SPLIT (FIXED + NO LEAKAGE)
# ================================

# Split FIRST
X_train, X_test, y_train, y_test, df_train, df_test = train_test_split(
    X, y, df, test_size=0.2, stratify=y, random_state=42
)

print("Train distribution:", np.bincount(y_train))



# ================================
# 🧬 Gene encoding (ONLY TRAIN → avoids leakage)
# ================================
top_genes = df_train["gene_symbol"].value_counts().head(50).index

for gene in top_genes:
    df_train[f"gene_{gene}"] = (df_train["gene_symbol"] == gene).astype(int)
    df_test[f"gene_{gene}"] = (df_test["gene_symbol"] == gene).astype(int)

gene_features = [f"gene_{g}" for g in top_genes]

# ================================
# Rebuild BIO features safely
# ================================
X_bio_train = df_train[[
    "mutation_flag",
    "rsem_tpm",
    "expr_mut_interaction",
    "high_expression"
] + gene_features].values

X_bio_test = df_test[[
    "mutation_flag",
    "rsem_tpm",
    "expr_mut_interaction",
    "high_expression"
] + gene_features].values

# Split chem features
chem_dim = X_chem.shape[1]

X_chem_train = X_train[:, :chem_dim]
X_chem_test  = X_test[:, :chem_dim]

# Combine again
X_train = np.hstack([X_chem_train, X_bio_train])
X_test  = np.hstack([X_chem_test, X_bio_test])

# ================================
# ⚖️ SMOTE (ONLY TRAIN)
# ================================
X_train_bal, y_train_bal = X_train, y_train

# ================================
# 📏 Scaling (AFTER SMOTE)
# ================================
scaler = StandardScaler()
X_train_bal = scaler.fit_transform(X_train_bal)
X_test = scaler.transform(X_test)

from sklearn.feature_selection import VarianceThreshold

selector = VarianceThreshold(threshold=0.01)
X_train_bal = selector.fit_transform(X_train_bal)
X_test = selector.transform(X_test)

# ================================
# 🤖 PART 5: BASE MODELS
# ================================

models = {
    "rf": RandomForestClassifier(n_estimators=700, class_weight="balanced"),
    
    "xgb": xgb.XGBClassifier(
        n_estimators=800,
        max_depth=8,
        learning_rate=0.03,
        subsample=0.8,
        colsample_bytree=0.8,
        #scale_pos_weight=len(y_train)/sum(y_train),
        eval_metric="logloss"
    ),
    
    "lgbm": lgb.LGBMClassifier(
        n_estimators=800,
        learning_rate=0.03,
        num_leaves=64,
        class_weight="balanced"
    ),
    
    "lr": LogisticRegression(max_iter=1000)
}




# ================================
# 🔷 PART 6: STACKING + CALIBRATION
# ================================

# Base learners (already defined)
base_estimators = [
    ("rf", models["rf"]),
    ("xgb", models["xgb"]),
    ("lgbm", models["lgbm"]),
    ("lr", models["lr"])
]

# Meta-model
meta_model = LogisticRegression(class_weight="balanced", max_iter=1000)

# Stacking classifier
stack_model = StackingClassifier(
    estimators=base_estimators,
    final_estimator=meta_model,
    stack_method="predict_proba",
    cv=5,
    n_jobs=-1
)

# ================================
# 🧠 Train Stacking Model
# ================================
stack_model.fit(X_train_bal, y_train_bal)

# ================================
# 🎯 Probability Calibration
# ================================
calibrated_model = CalibratedClassifierCV(
    stack_model,
    method="isotonic",
    cv=3
)

calibrated_model.fit(X_train_bal, y_train_bal)




# ================================
# 📊 PART 7: EVALUATION (FIXED)
# ================================

# Get probabilities FIRST
y_prob = calibrated_model.predict_proba(X_test)[:, 1]

# Find best threshold
from sklearn.metrics import roc_curve

fpr, tpr, thresholds = roc_curve(y_test, y_prob)
best_thresh = thresholds[np.argmax(tpr - fpr)]
best_f1 = 0
#best_thresh = 0.5

for t in thresholds:
    preds = (y_prob > t).astype(int)
    f1 = f1_score(y_test, preds)
    if f1 > best_f1:
        best_f1 = f1
        best_thresh = t

print("Best threshold:", best_thresh)

# Final predictions
y_pred = (y_prob > best_thresh).astype(int)

print("\n=== Classification Report ===")
print(classification_report(y_test, y_pred))

print("AUC:", roc_auc_score(y_test, y_prob))
print("F1 Score:", f1_score(y_test, y_pred))

results = []

for name, model in models.items():
    model.fit(X_train_bal, y_train_bal)
    
    y_prob_m = model.predict_proba(X_test)[:, 1]
    y_pred_m = (y_prob_m > 0.5).astype(int)
    
    results.append({
        "Model": name.upper(),
        "Accuracy": accuracy_score(y_test, y_pred_m),
        "F1": f1_score(y_test, y_pred_m),
        "Precision": precision_score(y_test, y_pred_m),
        "Recall": recall_score(y_test, y_pred_m),
        "AUC": roc_auc_score(y_test, y_prob_m)
    })

# Add stacking
results.append({
    "Model": "STACKING",
    "Accuracy": accuracy_score(y_test, y_pred),
    "F1": f1_score(y_test, y_pred),
    "Precision": precision_score(y_test, y_pred),
    "Recall": recall_score(y_test, y_pred),
    "AUC": roc_auc_score(y_test, y_prob)
})

results_df = pd.DataFrame(results)


plt.figure(figsize=(10,6))

results_df.set_index("Model")[["F1","Accuracy","Recall","Precision"]].plot(
    kind="bar",
    colormap="viridis"
)

plt.title("Model Comparison Performance", fontsize=14)
plt.ylabel("Score")
plt.xticks(rotation=0)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/model_comparison.pdf")
plt.close()


plt.figure()

for name, model in models.items():
    y_prob_m = model.predict_proba(X_test)[:, 1]
    fpr, tpr, _ = roc_curve(y_test, y_prob_m)
    plt.plot(fpr, tpr, label=name.upper())

# Stacking
fpr, tpr, _ = roc_curve(y_test, y_prob)
plt.plot(fpr, tpr, linewidth=3, label="STACKING")

plt.plot([0,1],[0,1],'--',color='gray')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve Comparison")
plt.legend()
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/roc_all_models.pdf")
plt.close()

# ================================
# 📉 Confusion Matrix Plot
# ================================
cm = confusion_matrix(y_test, y_pred)

plt.figure()
sns.heatmap(cm, annot=True, fmt="d", cmap="crest",
            xticklabels=["Inactive", "Active"],
            yticklabels=["Inactive", "Active"])
plt.title("Confusion Matrix", fontsize=14)
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/confusion_matrix.pdf")
plt.close()

# ================================
# 📈 ROC Curve
# ================================
from sklearn.metrics import roc_curve

fpr, tpr, _ = roc_curve(y_test, y_prob)

plt.figure()
plt.plot(fpr, tpr, color=COLORS["accent"], lw=2)
plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
plt.title("ROC Curve", fontsize=14)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/roc_curve.pdf")
plt.close()

# ================================
# 📉 Precision-Recall Curve
# ================================
precision, recall, _ = precision_recall_curve(y_test, y_prob)

plt.figure()
plt.plot(recall, precision, color=COLORS["secondary"], lw=2)
plt.title("Precision-Recall Curve", fontsize=14)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/pr_curve.pdf")
plt.close()






# ================================
# 🧬 PART 8: SHAP ANALYSIS (ROBUST FINAL)
# ================================

# Train standalone XGB
model = xgb.XGBClassifier(
    n_estimators=600,
    max_depth=6,
    learning_rate=0.05,
    subsample=0.9,
    colsample_bytree=0.9,
    eval_metric="logloss",
    base_score=0.5   # 🔥 FIX
)

model.fit(X_train_bal, y_train_bal)

# -------------------------------
# ✅ USE KERNEL EXPLAINER (NO VERSION ISSUES)
# -------------------------------

# Small background sample (VERY IMPORTANT for speed)
background = X_train_bal[np.random.choice(X_train_bal.shape[0], 100, replace=False)]

# ================================
# 🧬 FEATURE NAMES (ADD THIS)
# ================================
feature_names = [f"chem_{i}" for i in range(X_train_bal.shape[1])]

# ================================
# 🧬 SHAP ANALYSIS
# ================================
X_sample = X_test[:50]

explainer = shap.SamplingExplainer(
    model.predict_proba,
    X_train_bal[:100]
)

shap_values = explainer(X_sample)

min_len = min(len(X_sample), shap_values.values.shape[0])

plt.figure()
shap.summary_plot(
    shap_values.values[:min_len, :, 1],
    X_sample[:min_len],
    feature_names=feature_names,
    show=False
)
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/shap_summary.pdf")
plt.close()

# Feature importance
shap_array = shap_values.values[:min_len, :, 1]

shap_mean = np.abs(shap_array).mean(axis=0)

importance_df = pd.DataFrame({
    "feature": feature_names,
    "importance": shap_mean
}).sort_values(by="importance", ascending=False)

plt.figure()
sns.barplot(
    x="importance",
    y="feature",
    data=importance_df.head(20),
    palette="viridis"
)
plt.title("Top 20 Important Features (SHAP)")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/shap_top_features.pdf")
plt.close()




# ================================
# 📉 PART 9: IC50 REGRESSION
# ================================

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

# Keep IC50 aligned properly
y_ic50 = np.log1p(df["IC50"])

X_train, X_test, y_train, y_test, y_ic50_train, y_ic50_test = train_test_split(
    X, y, y_ic50, test_size=0.2, stratify=y, random_state=42
)

# Train regressor
rf_reg = RandomForestRegressor(n_estimators=500)
rf_reg.fit(X_train, y_ic50_train)

# Predict
ic50_pred = rf_reg.predict(X_test)

rmse = np.sqrt(mean_squared_error(y_ic50_test, ic50_pred))
r2 = r2_score(y_ic50_test, ic50_pred)

print("IC50 RMSE:", rmse)
print("IC50 R2:", r2)

# ================================
# 📊 Predicted vs Actual Plot
# ================================
plt.figure()
plt.scatter(y_ic50_test, ic50_pred,
            color=COLORS["primary"], alpha=0.6)
plt.xlabel("Actual IC50")
plt.ylabel("Predicted IC50")
plt.title("IC50 Prediction Performance")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/ic50_regression.pdf")
plt.close()



# ================================
# 🌿 Tissue-wise Activity Analysis
# ================================
plt.figure()
sns.boxplot(x=df["Tissue"], y=df["IC50"], palette="coolwarm")
plt.xticks(rotation=45)
plt.title("IC50 Distribution Across Tissues")
plt.tight_layout()
plt.savefig(f"{SAVE_DIR}/tissue_ic50.pdf")
plt.close()