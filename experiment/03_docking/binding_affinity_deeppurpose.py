"""
Step 3: DeepPurpose (DeepDTA) 기반 Drug-Target Interaction 예측
- EGFR 단백질 시퀀스 + 화합물 SMILES → 결합 친화도 (Kd) 예측
- 사전학습 모델 사용 (Davis kinase dataset)
"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# DeepPurpose
from DeepPurpose import DTI as models
from DeepPurpose import utils

# EGFR 키나제 도메인 아미노산 시퀀스 (UniProt P00533, residues 696-1022)
EGFR_SEQUENCE = (
    "FKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDN"
    "PHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDR"
    "RLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIY"
    "THQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCW"
    "MIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVV"
    "DADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDR"
)

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "02_filtering")
INPUT_FILE = os.path.join(INPUT_DIR, "egfr_filtered.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)


def predict_binding_affinity(smiles_list, target_seq, batch_size=256):
    """DeepPurpose 사전학습 모델로 결합 친화도 예측"""
    print("Loading pretrained DeepDTA model (Davis dataset)...")
    # CNN (drug) + CNN (target) 조합 사전학습 모델
    model = models.model_pretrained(model="MPNN_CNN_DAVIS")

    target_list = [target_seq] * len(smiles_list)

    print(f"Predicting binding affinity for {len(smiles_list)} compounds...")
    # DeepPurpose는 Kd 값을 예측 (nM 단위)
    predictions = models.virtual_screening(
        smiles_list, target_list, model,
        drug_names=None, target_names=None,
        verbose=False
    )

    return predictions


def main():
    print("Loading filtered compounds...")
    df = pd.read_csv(INPUT_FILE)
    print(f"Input: {len(df)} compounds")

    # 너무 많으면 상위 활성 화합물 위주로 예측 (시간 절약)
    # pChEMBL >= 6 (IC50 <= 1 uM) 화합물 + 랜덤 샘플
    df_high = df[df["pchembl_value"] >= 6].copy()
    df_low = df[df["pchembl_value"] < 6].sample(n=min(500, len(df[df["pchembl_value"] < 6])), random_state=42)
    df_subset = pd.concat([df_high, df_low]).drop_duplicates(subset=["canonical_smiles"]).reset_index(drop=True)
    print(f"Subset for prediction: {len(df_subset)} compounds (high activity: {len(df_high)}, random low: {len(df_low)})")

    smiles_list = df_subset["canonical_smiles"].tolist()

    # 결합 친화도 예측
    predictions = predict_binding_affinity(smiles_list, EGFR_SEQUENCE)

    # 결과 파싱 — predictions는 DataFrame
    if isinstance(predictions, pd.DataFrame):
        pred_values = predictions.iloc[:, -1].values  # 마지막 컬럼이 예측값
    else:
        pred_values = np.array(predictions)

    df_subset["predicted_Kd"] = pred_values
    # pKd = -log10(Kd) where Kd in M → higher = stronger binding
    # DeepPurpose DAVIS 모델은 Kd (nM)를 직접 예측
    df_subset["predicted_pKd"] = -np.log10(df_subset["predicted_Kd"].clip(lower=1e-3) * 1e-9)

    # 정렬 (pKd 높은 순 = 강한 결합)
    df_subset = df_subset.sort_values("predicted_pKd", ascending=False).reset_index(drop=True)

    # 저장
    output_file = os.path.join(OUTPUT_DIR, "egfr_binding_affinity.csv")
    df_subset.to_csv(output_file, index=False)

    # === 시각화 ===

    # 1) 예측 pKd 분포
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("DeepDTA Binding Affinity Prediction Results", fontsize=14, fontweight="bold")

    axes[0].hist(df_subset["predicted_pKd"], bins=50, color="#2ca02c", edgecolor="white", alpha=0.8)
    axes[0].set_xlabel("Predicted pKd")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Distribution of Predicted pKd")

    # 2) 실험값 vs 예측값 비교
    # pChEMBL ≈ pIC50, pKd와 직접 비교는 어렵지만 상관관계 확인
    axes[1].scatter(df_subset["pchembl_value"], df_subset["predicted_pKd"],
                    alpha=0.3, s=10, c="#1f77b4")
    axes[1].set_xlabel("Experimental pChEMBL (pIC50)")
    axes[1].set_ylabel("Predicted pKd (DeepDTA)")
    axes[1].set_title("Experimental vs Predicted Affinity")

    # 상관계수
    corr = df_subset[["pchembl_value", "predicted_pKd"]].corr().iloc[0, 1]
    axes[1].text(0.05, 0.95, f"r = {corr:.3f}", transform=axes[1].transAxes,
                 fontsize=12, verticalalignment="top",
                 bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # 3) Top 50 화합물의 pKd
    top50 = df_subset.head(50)
    axes[2].barh(range(len(top50)), top50["predicted_pKd"].values, color="#ff7f0e", alpha=0.8)
    axes[2].set_xlabel("Predicted pKd")
    axes[2].set_ylabel("Compound Rank")
    axes[2].set_title("Top 50 Compounds by Predicted pKd")
    axes[2].invert_yaxis()

    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "fig03_binding_affinity.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {fig_path}")

    # 요약
    print(f"\n=== Binding Affinity Prediction Summary ===")
    print(f"Total predictions: {len(df_subset)}")
    print(f"Predicted pKd range: {df_subset['predicted_pKd'].min():.2f} - {df_subset['predicted_pKd'].max():.2f}")
    print(f"Correlation (pChEMBL vs pKd): {corr:.3f}")
    print(f"\nTop 10 compounds:")
    top10 = df_subset.head(10)[["molecule_chembl_id", "pchembl_value", "predicted_pKd", "MW", "logP", "max_tanimoto_sim"]]
    print(top10.to_string(index=False))
    print(f"\nOutput: {output_file}")


if __name__ == "__main__":
    main()
