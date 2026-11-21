"""
Step 4: ChemBERTa (Transformer) 기반 분자 활성 분류
- 사전학습된 ChemBERTa 모델로 분자 임베딩 추출
- PCA 차원 축소 후 RandomForest 분류/회귀
- t-SNE 시각화
"""
import os
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_score, cross_val_predict, StratifiedKFold
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "02_filtering")
INPUT_FILE = os.path.join(INPUT_DIR, "egfr_filtered.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

MODEL_NAME = "seyonec/ChemBERTa-zinc-base-v1"
EMB_CACHE = os.path.join(OUTPUT_DIR, "chemberta_embeddings.npy")


def extract_embeddings(smiles_list, batch_size=64):
    """ChemBERTa로 분자 임베딩 추출 (CLS token)"""
    if os.path.exists(EMB_CACHE):
        print(f"Loading cached embeddings from {EMB_CACHE}")
        return np.load(EMB_CACHE)

    print(f"Loading ChemBERTa model: {MODEL_NAME}")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME)
    model.eval()

    all_embeddings = []
    total = len(smiles_list)

    for i in range(0, total, batch_size):
        batch = smiles_list[i:i+batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=512)
        with torch.no_grad():
            outputs = model(**inputs)
            cls_embeddings = outputs.last_hidden_state[:, 0, :].numpy()
        all_embeddings.append(cls_embeddings)
        if (i // batch_size) % 20 == 0:
            print(f"  Embedding: {min(i+batch_size, total)}/{total}")

    result = np.vstack(all_embeddings)
    np.save(EMB_CACHE, result)
    print(f"Embeddings saved: {EMB_CACHE}")
    return result


def main():
    print("Loading data...")
    df = pd.read_csv(INPUT_FILE)
    print(f"Total compounds: {len(df)}")

    smiles_list = df["canonical_smiles"].tolist()
    embeddings = extract_embeddings(smiles_list)
    print(f"Embedding shape: {embeddings.shape}")

    # PCA 차원 축소 (768 → 50)
    print("PCA dimensionality reduction (768 → 50)...")
    pca = PCA(n_components=50, random_state=42)
    emb_pca = pca.fit_transform(embeddings)
    print(f"  Explained variance: {pca.explained_variance_ratio_.sum():.3f}")

    # 활성 라벨: pChEMBL >= 7 (IC50 <= 100 nM) = active
    threshold = 7.0
    df["activity_label"] = (df["pchembl_value"] >= threshold).astype(int)
    labels = df["activity_label"].values
    n_active = labels.sum()
    print(f"Active (pChEMBL >= {threshold}): {n_active} / {len(labels)} ({n_active/len(labels)*100:.1f}%)")

    # 분류 — RandomForest (빠름)
    print("Training RandomForest classifier (5-fold CV)...")
    clf = RandomForestClassifier(n_estimators=300, max_depth=20, n_jobs=-1, random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(clf, emb_pca, labels, cv=cv, scoring="roc_auc", n_jobs=-1)
    print(f"  ROC-AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

    # Out-of-fold 예측으로 데이터 누수 방지
    probs = cross_val_predict(clf, emb_pca, labels, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
    df["chemberta_active_prob"] = probs
    # 최종 모델 학습 (시각화 및 추후 사용)
    clf.fit(emb_pca, labels)

    # 회귀 — pChEMBL 예측
    print("Training RandomForest regressor (5-fold CV)...")
    reg = RandomForestRegressor(n_estimators=300, max_depth=15, min_samples_leaf=5, n_jobs=-1, random_state=42)
    reg_scores = cross_val_score(reg, emb_pca, df["pchembl_value"].values, cv=5, scoring="r2", n_jobs=-1)
    print(f"  R²: {reg_scores.mean():.3f} ± {reg_scores.std():.3f}")

    # Out-of-fold 예측으로 데이터 누수 방지
    pchembl_pred = cross_val_predict(reg, emb_pca, df["pchembl_value"].values, cv=5, n_jobs=-1)
    df["chemberta_predicted_pchembl"] = pchembl_pred
    # 최종 모델 학습 (시각화 및 추후 사용)
    reg.fit(emb_pca, df["pchembl_value"].values)

    # === 시각화 ===
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("ChemBERTa Molecular Activity Analysis", fontsize=14, fontweight="bold")

    # 1) t-SNE (PCA 50차원에서)
    print("Computing t-SNE (n=1500)...")
    n_sample = min(1500, len(emb_pca))
    idx = np.random.RandomState(42).choice(len(emb_pca), n_sample, replace=False)
    tsne = TSNE(n_components=2, random_state=42, perplexity=30, max_iter=500, learning_rate="auto", init="pca")
    tsne_coords = tsne.fit_transform(emb_pca[idx])

    scatter = axes[0, 0].scatter(tsne_coords[:, 0], tsne_coords[:, 1],
                                  c=df["pchembl_value"].values[idx], cmap="RdYlGn",
                                  s=8, alpha=0.6, vmin=5, vmax=10)
    plt.colorbar(scatter, ax=axes[0, 0], label="pChEMBL Value")
    axes[0, 0].set_title("t-SNE of ChemBERTa Embeddings")
    axes[0, 0].set_xlabel("t-SNE 1")
    axes[0, 0].set_ylabel("t-SNE 2")

    # 2) CV ROC-AUC
    axes[0, 1].bar(range(1, 6), cv_scores, color="#2ca02c", edgecolor="white", width=0.6)
    axes[0, 1].axhline(cv_scores.mean(), color="red", linestyle="--", label=f"Mean: {cv_scores.mean():.3f}")
    axes[0, 1].set_xlabel("Fold")
    axes[0, 1].set_ylabel("ROC-AUC")
    axes[0, 1].set_title("5-Fold Cross-Validation (Classification)")
    axes[0, 1].legend()
    axes[0, 1].set_ylim(0.5, 1.0)

    # 3) 예측 vs 실험
    axes[1, 0].scatter(df["pchembl_value"], df["chemberta_predicted_pchembl"], alpha=0.2, s=5, c="#1f77b4")
    lims = [5, 11]
    axes[1, 0].plot(lims, lims, "r--", linewidth=1)
    axes[1, 0].set_xlabel("Experimental pChEMBL")
    axes[1, 0].set_ylabel("Predicted pChEMBL")
    axes[1, 0].set_title(f"Regression (CV R² = {reg_scores.mean():.3f})")
    axes[1, 0].set_xlim(lims)
    axes[1, 0].set_ylim(lims)

    # 4) 활성 확률 분포
    axes[1, 1].hist(probs[labels == 0], bins=50, alpha=0.6, color="#d62728", label=f"Inactive ({(labels==0).sum()})")
    axes[1, 1].hist(probs[labels == 1], bins=50, alpha=0.6, color="#2ca02c", label=f"Active ({(labels==1).sum()})")
    axes[1, 1].set_xlabel("Predicted Active Probability")
    axes[1, 1].set_ylabel("Count")
    axes[1, 1].set_title("Activity Classification")
    axes[1, 1].legend()

    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "fig04_chemberta_analysis.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {fig_path}")

    # 결과 저장
    output_file = os.path.join(OUTPUT_DIR, "egfr_chemberta_scores.csv")
    df.to_csv(output_file, index=False)

    print(f"\n=== ChemBERTa Analysis Summary ===")
    print(f"Classification ROC-AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    print(f"Regression R²: {reg_scores.mean():.3f} ± {reg_scores.std():.3f}")
    print(f"PCA explained variance: {pca.explained_variance_ratio_.sum():.3f}")
    print(f"\nTop 10 by active probability:")
    top10 = df.nlargest(10, "chemberta_active_prob")[
        ["molecule_chembl_id", "pchembl_value", "chemberta_active_prob", "chemberta_predicted_pchembl"]
    ]
    print(top10.to_string(index=False))
    print(f"\nOutput: {output_file}")


if __name__ == "__main__":
    main()
