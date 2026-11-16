"""
Step 5: ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) 예측
- admet_ai 패키지 사용 (Chemprop 기반 딥러닝 모델)
- 약물의 약동학적 성질 평가
"""
import os
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
from admet_ai import ADMETModel
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "03_docking")
BINDING_FILE = os.path.join(INPUT_DIR, "egfr_binding_affinity.csv")
CHEMBERTA_FILE = os.path.join(INPUT_DIR, "egfr_chemberta_scores.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# 주요 ADMET 속성 (보고서에 쓸 것)
KEY_PROPERTIES = {
    # Absorption
    "HIA_Hou": "Human Intestinal Absorption",
    "Caco2_Wang": "Caco-2 Permeability",
    "Bioavailability_Ma": "Oral Bioavailability",
    "Pgp_Broccatelli": "P-gp Inhibitor",
    # Distribution
    "BBB_Martins": "BBB Penetration",
    "VDss_Lombardo": "Volume of Distribution",
    # Metabolism
    "CYP2D6_Veith": "CYP2D6 Inhibitor",
    "CYP3A4_Veith": "CYP3A4 Inhibitor",
    "CYP2C9_Veith": "CYP2C9 Inhibitor",
    # Excretion
    "Half_Life_Obach": "Half-Life",
    "Clearance_Hepatocyte_AZ": "Hepatocyte Clearance",
    # Toxicity
    "hERG": "hERG Inhibition",
    "AMES": "AMES Mutagenicity",
    "DILI": "Drug-Induced Liver Injury",
    "LD50_Zhu": "Acute Toxicity (LD50)",
    "Skin_Reaction": "Skin Sensitization",
}


def predict_admet(smiles_list):
    """ADMET 예측"""
    print("Loading ADMET-AI model...")
    model = ADMETModel()

    print(f"Predicting ADMET properties for {len(smiles_list)} compounds...")
    results = model.predict(smiles=smiles_list)
    return results


def create_admet_radar(df_top, fig_path):
    """상위 후보 화합물의 ADMET radar chart"""
    # 약물 적합성 점수 (높을수록 좋음)
    radar_props = {
        "HIA_Hou": ("Absorption", True),  # True = 높을수록 좋음
        "Bioavailability_Ma": ("Bioavail.", True),
        "BBB_Martins": ("BBB", False),  # 항암제는 BBB 불필요
        "hERG": ("hERG Safe", False),  # 낮을수록 좋음 (non-inhibitor)
        "AMES": ("Non-Mutagen", False),
        "DILI": ("Liver Safe", False),
    }

    available = [k for k in radar_props if k in df_top.columns]
    if len(available) < 3:
        print("Not enough ADMET properties for radar chart, skipping.")
        return

    n_compounds = min(5, len(df_top))
    if n_compounds == 0:
        print("No compounds for radar chart.")
        return
    fig, axes = plt.subplots(1, n_compounds, figsize=(4*n_compounds, 4),
                              subplot_kw=dict(polar=True))
    if n_compounds == 1:
        axes = [axes]

    for idx, (_, row) in enumerate(df_top.head(n_compounds).iterrows()):
        ax = axes[idx]
        values = []
        labels = []
        for prop in available:
            label, is_good_high = radar_props[prop]
            val = row.get(prop, 0.5)
            if not is_good_high:
                val = 1.0 - val  # 반전 (낮을수록 좋은 것 → 높은 점수로)
            values.append(float(val) if pd.notna(val) else 0.5)
            labels.append(label)

        values.append(values[0])  # 닫기
        angles = np.linspace(0, 2*np.pi, len(labels), endpoint=False).tolist()
        angles.append(angles[0])

        ax.plot(angles, values, "o-", linewidth=2)
        ax.fill(angles, values, alpha=0.25)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylim(0, 1)
        ax.set_title(row["molecule_chembl_id"], fontsize=9, fontweight="bold", pad=15)

    plt.suptitle("ADMET Profile — Top Candidates", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Radar chart saved: {fig_path}")


def create_admet_heatmap(df_top, fig_path):
    """ADMET 히트맵"""
    display_props = [k for k in KEY_PROPERTIES if k in df_top.columns]
    if not display_props:
        print("No ADMET properties found for heatmap.")
        return

    n_show = min(20, len(df_top))
    data = df_top.head(n_show)[display_props].values.astype(float)
    labels_y = df_top.head(n_show)["molecule_chembl_id"].values
    labels_x = [KEY_PROPERTIES[p] for p in display_props]

    fig, ax = plt.subplots(figsize=(max(12, len(display_props)*0.8), max(6, n_show*0.4)))
    im = ax.imshow(data, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(labels_x)))
    ax.set_xticklabels(labels_x, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_show))
    ax.set_yticklabels(labels_y, fontsize=8)
    plt.colorbar(im, ax=ax, label="Score (0=bad, 1=good)", shrink=0.8)
    ax.set_title("ADMET Property Heatmap — Top Candidates", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved: {fig_path}")


def main():
    # 두 소스 병합 (binding affinity + chemberta)
    print("Loading and merging data...")
    df_bind = pd.read_csv(BINDING_FILE)
    df_chem = pd.read_csv(CHEMBERTA_FILE)

    # ChemBERTa 결과를 기준으로 binding affinity 병합
    df = df_chem.copy()
    bind_cols = ["molecule_chembl_id", "predicted_Kd", "predicted_pKd"]
    bind_cols = [c for c in bind_cols if c in df_bind.columns]
    if "predicted_pKd" in df_bind.columns:
        bind_data = df_bind[["molecule_chembl_id", "predicted_pKd"]].drop_duplicates(subset=["molecule_chembl_id"])
        df = df.merge(bind_data, on="molecule_chembl_id", how="left")

    # 상위 화합물 선택 (높은 활성 확률 + 높은 pChEMBL)
    df["combined_score"] = (
        0.4 * df["pchembl_value"] / df["pchembl_value"].max() +
        0.3 * df["chemberta_active_prob"] +
        0.3 * (df.get("predicted_pKd", df["pchembl_value"]) / df.get("predicted_pKd", df["pchembl_value"]).max())
    )
    df = df.sort_values("combined_score", ascending=False).reset_index(drop=True)

    # 상위 500개만 ADMET 예측 (시간 절약)
    n_predict = min(500, len(df))
    print(f"Selecting top {n_predict} compounds for ADMET prediction...")
    df_top = df.head(n_predict).copy()

    # ADMET 예측
    admet_results = predict_admet(df_top["canonical_smiles"].tolist())

    # 결과 병합
    if isinstance(admet_results, pd.DataFrame):
        admet_results = admet_results.reset_index(drop=True)
        df_top = df_top.reset_index(drop=True)
        df_top = pd.concat([df_top, admet_results], axis=1)

    # ADMET 필터링
    # 좋은 약물: HIA 높음, hERG 낮음, AMES 낮음, DILI 낮음
    # ADMET 종합 점수 계산 (가중 평균)
    admet_score = pd.Series(0.0, index=df_top.index)
    n_props = 0
    good_high = ["HIA_Hou", "Caco2_Wang", "Bioavailability_Ma"]  # 높을수록 좋음
    good_low = ["hERG", "AMES", "DILI", "CYP2D6_Veith"]  # 낮을수록 좋음
    for prop in good_high:
        if prop in df_top.columns:
            admet_score += df_top[prop].fillna(0.5)
            n_props += 1
    for prop in good_low:
        if prop in df_top.columns:
            admet_score += (1 - df_top[prop].fillna(0.5))
            n_props += 1
    if n_props > 0:
        admet_score /= n_props
    df_top["admet_score"] = admet_score

    # 상위 50%만 통과 (적응적 필터)
    admet_threshold = df_top["admet_score"].median()
    admet_pass = df_top["admet_score"] >= admet_threshold

    df_passed = df_top[admet_pass].copy()
    print(f"ADMET filter: {len(df_top)} → {len(df_passed)} compounds passed")

    # 시각화
    create_admet_heatmap(df_passed, os.path.join(FIGURES_DIR, "fig05_admet_heatmap.png"))
    create_admet_radar(df_passed, os.path.join(FIGURES_DIR, "fig06_admet_radar.png"))

    # ADMET 분포 시각화
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("ADMET Property Distributions", fontsize=14, fontweight="bold")
    plot_props = ["HIA_Hou", "Bioavailability_Ma", "BBB_Martins", "hERG", "AMES", "DILI"]
    colors = ["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728", "#9467bd", "#e377c2"]
    for ax, prop, color in zip(axes.flat, plot_props, colors):
        if prop in df_top.columns:
            ax.hist(df_top[prop].dropna(), bins=30, color=color, edgecolor="white", alpha=0.8)
            ax.set_title(KEY_PROPERTIES.get(prop, prop))
            ax.set_xlabel("Score")
            ax.set_ylabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig07_admet_distributions.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # 저장
    output_all = os.path.join(OUTPUT_DIR, "egfr_admet_all.csv")
    output_passed = os.path.join(OUTPUT_DIR, "egfr_admet_passed.csv")
    df_top.to_csv(output_all, index=False)
    df_passed.to_csv(output_passed, index=False)

    print(f"\n=== ADMET Prediction Summary ===")
    print(f"Input compounds: {n_predict}")
    print(f"ADMET passed: {len(df_passed)}")
    for prop in plot_props:
        if prop in df_top.columns:
            print(f"  {KEY_PROPERTIES.get(prop, prop)}: mean={df_top[prop].mean():.3f}")
    print(f"\nTop 10 ADMET-passed compounds:")
    show_cols = ["molecule_chembl_id", "pchembl_value", "combined_score"]
    show_cols += [p for p in ["HIA_Hou", "hERG", "AMES", "DILI"] if p in df_passed.columns]
    print(df_passed.head(10)[show_cols].to_string(index=False))
    print(f"\nOutput: {output_all}")
    print(f"Passed: {output_passed}")


if __name__ == "__main__":
    main()
