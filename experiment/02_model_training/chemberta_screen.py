"""
Step 2: ChemBERTa 기반 EGFR 활성 예측 모델 학습 + Virtual Screening

1) EGFR 학습 데이터로 분류 모델 학습 (Active vs Inactive)
2) 학습된 모델로 스크리닝 라이브러리 47K개 활성 예측
3) 상위 ~1,000개 선별 → 다음 단계(ADMET)로 전달

모델: ChemBERTa (seyonec/ChemBERTa-zinc-base-v1)
  - ZINC 7700만 분자로 사전학습된 Transformer
  - CLS 토큰 768차원 임베딩 → PCA(50) → RandomForest 분류
"""
import os
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, cross_val_predict, StratifiedKFold
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import roc_auc_score, classification_report
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import joblib

PREPROCESS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "01_preprocessing")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

MODEL_NAME = "seyonec/ChemBERTa-zinc-base-v1"
TOP_N_SCREENING = 1000  # 스크리닝 상위 N개 선별


def extract_embeddings(smiles_list, batch_size=64, cache_path=None):
    """ChemBERTa CLS 토큰 임베딩 추출"""
    if cache_path and os.path.exists(cache_path):
        print(f"  Loading cached embeddings: {cache_path}")
        return np.load(cache_path)

    print(f"  Loading ChemBERTa: {MODEL_NAME}")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME)
    model.eval()

    all_emb = []
    total = len(smiles_list)
    for i in range(0, total, batch_size):
        batch = smiles_list[i:i+batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True,
                           truncation=True, max_length=512)
        with torch.no_grad():
            outputs = model(**inputs)
            cls = outputs.last_hidden_state[:, 0, :].numpy()
        all_emb.append(cls)
        if (i // batch_size) % 25 == 0:
            print(f"    {min(i+batch_size, total)}/{total}")

    result = np.vstack(all_emb)
    if cache_path:
        np.save(cache_path, result)
        print(f"  Cached: {cache_path}")
    return result


def train_and_evaluate(emb_pca, labels):
    """분류 모델 학습 + cross-validation 평가"""
    print("\n=== Model Training (5-fold CV) ===")
    clf = RandomForestClassifier(n_estimators=300, max_depth=20, n_jobs=-1, random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # Cross-validation 성능
    cv_scores = cross_val_score(clf, emb_pca, labels, cv=cv, scoring="roc_auc", n_jobs=-1)
    print(f"  ROC-AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

    # Out-of-fold 예측 (데이터 누수 방지)
    oof_probs = cross_val_predict(clf, emb_pca, labels, cv=cv,
                                   method="predict_proba", n_jobs=-1)[:, 1]

    # 최종 모델 학습
    clf.fit(emb_pca, labels)

    return clf, cv_scores, oof_probs


def plot_training_results(df_train, oof_probs, cv_scores, emb_pca, labels):
    """학습 결과 시각화"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("ChemBERTa EGFR Activity Model — Training Results",
                 fontsize=14, fontweight="bold")

    # 1) t-SNE
    print("  Computing t-SNE...")
    n_sample = min(2000, len(emb_pca))
    idx = np.random.RandomState(42).choice(len(emb_pca), n_sample, replace=False)
    tsne = TSNE(n_components=2, random_state=42, perplexity=30, learning_rate="auto", init="pca")
    coords = tsne.fit_transform(emb_pca[idx])

    colors = labels[idx]
    scatter = axes[0, 0].scatter(coords[:, 0], coords[:, 1], c=colors,
                                  cmap="RdYlGn", s=8, alpha=0.6)
    plt.colorbar(scatter, ax=axes[0, 0], label="Activity (0=Inactive, 1=Active)")
    axes[0, 0].set_title("t-SNE of ChemBERTa Embeddings")
    axes[0, 0].set_xlabel("t-SNE 1")
    axes[0, 0].set_ylabel("t-SNE 2")

    # 2) CV ROC-AUC
    axes[0, 1].bar(range(1, 6), cv_scores, color="#2ca02c", edgecolor="white", width=0.6)
    axes[0, 1].axhline(cv_scores.mean(), color="red", linestyle="--",
                        label=f"Mean: {cv_scores.mean():.3f}")
    axes[0, 1].set_xlabel("Fold")
    axes[0, 1].set_ylabel("ROC-AUC")
    axes[0, 1].set_title("5-Fold Cross-Validation")
    axes[0, 1].legend()
    axes[0, 1].set_ylim(0.5, 1.0)

    # 3) 활성 확률 분포
    axes[1, 0].hist(oof_probs[labels == 0], bins=50, alpha=0.6,
                     color="#d62728", label=f"Inactive ({(labels==0).sum()})")
    axes[1, 0].hist(oof_probs[labels == 1], bins=50, alpha=0.6,
                     color="#2ca02c", label=f"Active ({(labels==1).sum()})")
    axes[1, 0].set_xlabel("Predicted Active Probability (OOF)")
    axes[1, 0].set_ylabel("Count")
    axes[1, 0].set_title("Activity Classification (Out-of-Fold)")
    axes[1, 0].legend()

    # 4) pChEMBL vs 예측 확률
    axes[1, 1].scatter(df_train["pchembl_value"], oof_probs, alpha=0.2, s=5, c="#1f77b4")
    axes[1, 1].set_xlabel("Experimental pChEMBL Value")
    axes[1, 1].set_ylabel("Predicted Active Probability")
    axes[1, 1].set_title("Experimental Activity vs Prediction")
    axes[1, 1].axhline(0.5, color="red", linestyle="--", alpha=0.5)
    axes[1, 1].axvline(7.0, color="green", linestyle="--", alpha=0.5, label="Active threshold")
    axes[1, 1].legend()

    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "fig02_chemberta_training.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {out}")


def plot_screening_results(df_screen_top, df_screen_all):
    """스크리닝 결과 시각화"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Virtual Screening Results — ChemBERTa", fontsize=14, fontweight="bold")

    # 1) 전체 활성 확률 분포
    axes[0].hist(df_screen_all["predicted_active_prob"], bins=50,
                  color="#1f77b4", edgecolor="white", alpha=0.8)
    axes[0].axvline(df_screen_top["predicted_active_prob"].min(), color="red",
                     linestyle="--", label=f"Top {len(df_screen_top)} cutoff")
    axes[0].set_xlabel("Predicted Active Probability")
    axes[0].set_ylabel("Count")
    axes[0].set_title(f"Screening Library ({len(df_screen_all)} compounds)")
    axes[0].legend()

    # 2) 상위 화합물 MW vs logP
    if "MW" in df_screen_top.columns:
        axes[1].scatter(df_screen_top["MW"], df_screen_top["logP"],
                        c=df_screen_top["predicted_active_prob"],
                        cmap="RdYlGn", s=15, alpha=0.7)
        axes[1].set_xlabel("Molecular Weight")
        axes[1].set_ylabel("logP")
        axes[1].set_title(f"Top {len(df_screen_top)} Hits — Drug-likeness")
        plt.colorbar(axes[1].collections[0], ax=axes[1], label="Active Prob.")

    # 3) 상위 화합물 유사도 분포
    if "max_tanimoto_sim" in df_screen_top.columns:
        axes[2].hist(df_screen_top["max_tanimoto_sim"], bins=30,
                      color="#ff7f0e", edgecolor="white", alpha=0.8)
        axes[2].set_xlabel("Tanimoto Similarity to Known EGFR Drugs")
        axes[2].set_ylabel("Count")
        axes[2].set_title("Structural Novelty of Hits")

    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "fig03_screening_results.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {out}")


def main():
    # === 1. 학습 데이터 로드 ===
    print("Loading training data...")
    df_train = pd.read_csv(os.path.join(PREPROCESS_DIR, "training_filtered.csv"))

    # Active/Inactive만 사용 (ambiguous 제외)
    df_train_labeled = df_train[df_train["activity_label"].isin([0, 1])].copy()
    df_train_labeled = df_train_labeled.reset_index(drop=True)
    labels = df_train_labeled["activity_label"].values
    print(f"  Labeled compounds: {len(df_train_labeled)} (Active: {labels.sum()}, Inactive: {(labels==0).sum()})")

    # === 2. 학습 데이터 임베딩 ===
    print("\nExtracting training embeddings...")
    train_emb = extract_embeddings(
        df_train_labeled["canonical_smiles"].tolist(),
        cache_path=os.path.join(OUTPUT_DIR, "train_embeddings.npy")
    )

    # PCA
    print("  PCA (768 → 50)...")
    pca = PCA(n_components=50, random_state=42)
    train_pca = pca.fit_transform(train_emb)
    print(f"  Explained variance: {pca.explained_variance_ratio_.sum():.3f}")

    # === 3. 모델 학습 ===
    clf, cv_scores, oof_probs = train_and_evaluate(train_pca, labels)

    # 학습 결과 시각화
    plot_training_results(df_train_labeled, oof_probs, cv_scores, train_pca, labels)

    # 모델 + PCA 저장
    joblib.dump(clf, os.path.join(OUTPUT_DIR, "chemberta_clf.pkl"))
    joblib.dump(pca, os.path.join(OUTPUT_DIR, "chemberta_pca.pkl"))
    print("  Model saved.")

    # === 4. 스크리닝 라이브러리 임베딩 + 예측 ===
    print("\nLoading screening library...")
    df_screen = pd.read_csv(os.path.join(PREPROCESS_DIR, "screening_filtered.csv"))
    print(f"  Screening compounds: {len(df_screen)}")

    print("Extracting screening embeddings...")
    screen_emb = extract_embeddings(
        df_screen["canonical_smiles"].tolist(),
        cache_path=os.path.join(OUTPUT_DIR, "screen_embeddings.npy")
    )

    screen_pca = pca.transform(screen_emb)
    df_screen["predicted_active_prob"] = clf.predict_proba(screen_pca)[:, 1]

    # 상위 N개 선별
    df_screen = df_screen.sort_values("predicted_active_prob", ascending=False)
    df_screen_top = df_screen.head(TOP_N_SCREENING).copy().reset_index(drop=True)

    # 스크리닝 결과 시각화
    plot_screening_results(df_screen_top, df_screen)

    # 저장
    df_screen.to_csv(os.path.join(OUTPUT_DIR, "screening_all_scored.csv"), index=False)
    df_screen_top.to_csv(os.path.join(OUTPUT_DIR, "screening_top_hits.csv"), index=False)

    print(f"\n=== Virtual Screening Summary ===")
    print(f"Training ROC-AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    print(f"Screening library: {len(df_screen)} compounds")
    print(f"Top {TOP_N_SCREENING} hits selected (prob >= {df_screen_top['predicted_active_prob'].min():.3f})")
    print(f"Top 10 hits:")
    show = df_screen_top.head(10)[["molecule_chembl_id", "predicted_active_prob", "MW", "logP"]]
    print(show.to_string(index=False))
    print(f"\nOutput: {OUTPUT_DIR}/screening_top_hits.csv")


if __name__ == "__main__":
    main()
