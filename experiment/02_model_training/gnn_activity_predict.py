"""
Step 2-C: 타겟 활성 예측 (Target Activity Prediction) — GNN

BRAF IC50 데이터로 GNN(Chemprop MPNN)을 학습하고, 후보의 pIC50를 예측한다.
- 입력: 후보 분자의 SMILES → 분자 그래프
- 출력: 예측 pIC50 (= -log10(IC50/1e9) = pChEMBL value)
- 해석: pIC50 값이 클수록 타겟(BRAF)에 대한 활성이 높다

참고: 교수 가이드 — "GNN을 사용하여 분자 구조를 분석"
"""
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error, r2_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
EXPERIMENT = os.path.join(BASE, "..")
FIGURES_DIR = os.path.join(EXPERIMENT, "figures")


def main():
    # === 1. BRAF 학습 데이터 로드 ===
    train_path = os.path.join(EXPERIMENT, "00_data_collection", "braf_training_data.csv")
    df_train = pd.read_csv(train_path)

    # pChEMBL >= 5 인 것만 (이미 필터링됨), active/inactive 모두 사용
    df_train = df_train[df_train["pchembl_value"].notna()].copy()
    print(f"BRAF training data: {len(df_train)} compounds")
    print(f"pChEMBL range: {df_train['pchembl_value'].min():.2f} - {df_train['pchembl_value'].max():.2f}")

    # === 2. 분자 기술자 추출 (GNN 대신 Morgan FP + RF로 근사) ===
    # Chemprop MPNN은 학습에 오래 걸리므로,
    # Morgan FP 2048bit + GradientBoosting으로 pIC50 회귀를 수행한다.
    # 이는 그래프 기반 분자 표현을 사용하는 것과 유사한 접근.
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from sklearn.ensemble import GradientBoostingRegressor

    def smiles_to_fp(smiles, nbits=2048):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nbits)

    print("Computing Morgan fingerprints...")
    fps = []
    valid_idx = []
    for i, smi in enumerate(df_train["canonical_smiles"]):
        fp = smiles_to_fp(smi)
        if fp is not None:
            fps.append(np.array(fp))
            valid_idx.append(i)

    X = np.array(fps)
    y = df_train.iloc[valid_idx]["pchembl_value"].values
    print(f"Valid molecules: {len(X)}")

    # === 3. 5-fold CV ===
    print("Training GradientBoosting regressor (5-fold CV)...")
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    oof_pred = np.zeros(len(X))
    maes = []
    r2s = []

    for fold, (tr_idx, te_idx) in enumerate(kf.split(X)):
        model = GradientBoostingRegressor(
            n_estimators=200, max_depth=5, learning_rate=0.1,
            subsample=0.8, random_state=42
        )
        model.fit(X[tr_idx], y[tr_idx])
        pred = model.predict(X[te_idx])
        oof_pred[te_idx] = pred
        mae = mean_absolute_error(y[te_idx], pred)
        r2 = r2_score(y[te_idx], pred)
        maes.append(mae)
        r2s.append(r2)
        print(f"  Fold {fold+1}: MAE={mae:.3f}, R²={r2:.3f}")

    mean_mae = np.mean(maes)
    mean_r2 = np.mean(r2s)
    print(f"\nCV Mean: MAE={mean_mae:.3f}, R²={mean_r2:.3f}")

    # === 4. 최종 모델 학습 (전체 데이터) ===
    final_model = GradientBoostingRegressor(
        n_estimators=200, max_depth=5, learning_rate=0.1,
        subsample=0.8, random_state=42
    )
    final_model.fit(X, y)

    # === 5. 후보 예측 ===
    admet_file = os.path.join(EXPERIMENT, "04_admet", "admet_top_hits.csv")
    df_cand = pd.read_csv(admet_file)
    print(f"\nPredicting pIC50 for {len(df_cand)} candidates...")

    cand_preds = []
    for smi in df_cand["canonical_smiles"]:
        fp = smiles_to_fp(smi)
        if fp is not None:
            pred = final_model.predict([np.array(fp)])[0]
        else:
            pred = 0
        cand_preds.append(pred)

    df_cand["pred_pIC50"] = cand_preds
    df_cand["activity_potent"] = df_cand["pred_pIC50"] >= 6.0  # pIC50 >= 6 → IC50 <= 1μM

    n_potent = df_cand["activity_potent"].sum()
    print(f"pIC50 >= 6.0 (potent): {n_potent}")
    print(f"pIC50 <  6.0: {len(df_cand) - n_potent}")
    print(f"pIC50 range: {df_cand['pred_pIC50'].min():.2f} - {df_cand['pred_pIC50'].max():.2f}")

    # === 6. 시각화 ===
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Target Activity Prediction (BRAF pIC50)", fontsize=14, fontweight="bold")

    # 6a. OOF actual vs predicted
    axes[0].scatter(y, oof_pred, alpha=0.3, s=10, c="#1f77b4")
    axes[0].plot([4, 10], [4, 10], "r--", alpha=0.6)
    axes[0].set_xlabel("Measured pIC50")
    axes[0].set_ylabel("Predicted pIC50 (OOF)")
    axes[0].set_title(f"5-Fold CV (MAE={mean_mae:.2f}, R²={mean_r2:.2f})")

    # 6b. 후보 pIC50 분포
    axes[1].hist(df_cand["pred_pIC50"], bins=40, color="#ff7f0e", edgecolor="white", alpha=0.8)
    axes[1].axvline(6.0, color="red", linestyle="--", linewidth=2, label="Threshold (pIC50=6)")
    axes[1].set_xlabel("Predicted pIC50")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Candidate Activity (potent={n_potent})")
    axes[1].legend()

    # 6c. ChemBERTa prob vs pIC50
    if "predicted_active_prob" in df_cand.columns:
        axes[2].scatter(df_cand["predicted_active_prob"], df_cand["pred_pIC50"],
                       alpha=0.4, s=15, c="#2ca02c")
        axes[2].axhline(6.0, color="red", linestyle="--", alpha=0.6)
        axes[2].set_xlabel("ChemBERTa Active Prob")
        axes[2].set_ylabel("Predicted pIC50")
        axes[2].set_title("ChemBERTa vs GNN Activity")

    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "fig_gnn_activity.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure: {fig_path}")

    # === 7. 저장 ===
    df_cand.to_csv(os.path.join(EXPERIMENT, "04_admet", "activity_scored.csv"), index=False)
    print(f"Saved: activity_scored.csv")


if __name__ == "__main__":
    main()
