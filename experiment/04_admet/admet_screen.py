"""
Step 3: ADMET-AI 약동학 평가 + 물리화학적 성질 평가

ChemBERTa 상위 1,000개 → ADMET 16속성 예측 → 상위 300개 선별
  - 물리화학적 성질: HIA(흡수), Caco-2(투과), Bioavailability(생체이용률)
  - 독성: hERG(심장), AMES(변이원성), DILI(간독성)
  - 대사: CYP 억제
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

INPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "02_model_training", "screening_top_hits.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

TOP_N = 300  # ADMET 통과 상위 N개 → GNINA 도킹 대상

# ADMET 속성 분류
ADMET_PROPERTIES = {
    # 물리화학적 성질 (Physicochemical)
    "HIA_Hou": ("Human Intestinal Absorption", "physico", True),
    "Caco2_Wang": ("Caco-2 Permeability", "physico", True),
    "Bioavailability_Ma": ("Oral Bioavailability", "physico", True),
    "Pgp_Broccatelli": ("P-gp Inhibitor", "physico", False),
    # 분포 (Distribution)
    "BBB_Martins": ("BBB Penetration", "dist", False),  # 항암제는 불필요
    "VDss_Lombardo": ("Volume of Distribution", "dist", None),
    # 대사 (Metabolism)
    "CYP2D6_Veith": ("CYP2D6 Inhibitor", "metab", False),
    "CYP3A4_Veith": ("CYP3A4 Inhibitor", "metab", False),
    "CYP2C9_Veith": ("CYP2C9 Inhibitor", "metab", False),
    # 배설 (Excretion)
    "Half_Life_Obach": ("Half-Life", "excr", None),
    "Clearance_Hepatocyte_AZ": ("Hepatocyte Clearance", "excr", None),
    # 독성 (Toxicity)
    "hERG": ("hERG Inhibition", "tox", False),
    "AMES": ("AMES Mutagenicity", "tox", False),
    "DILI": ("Drug-Induced Liver Injury", "tox", False),
    "LD50_Zhu": ("Acute Toxicity (LD50)", "tox", None),
    "Skin_Reaction": ("Skin Sensitization", "tox", False),
}


def compute_admet_score(df):
    """ADMET 종합 점수 계산"""
    score = pd.Series(0.0, index=df.index)
    n = 0
    for prop, (_, _, good_high) in ADMET_PROPERTIES.items():
        if prop not in df.columns or good_high is None:
            continue
        vals = df[prop].fillna(0.5)
        if good_high:
            score += vals
        else:
            score += (1 - vals)
        n += 1
    if n > 0:
        score /= n
    return score


def plot_admet_heatmap(df_top, fig_path):
    """ADMET 히트맵"""
    display_props = [k for k in ADMET_PROPERTIES if k in df_top.columns]
    if not display_props:
        return

    n_show = min(20, len(df_top))
    data = df_top.head(n_show)[display_props].values.astype(float)
    labels_y = df_top.head(n_show)["molecule_chembl_id"].values
    labels_x = [ADMET_PROPERTIES[p][0] for p in display_props]

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
    print(f"  Heatmap saved: {fig_path}")


def plot_admet_radar(df_top, fig_path):
    """상위 5개 후보 ADMET radar chart"""
    radar_props = {
        "HIA_Hou": ("Absorption", True),
        "Bioavailability_Ma": ("Bioavail.", True),
        "hERG": ("hERG Safe", False),
        "AMES": ("Non-Mutagen", False),
        "DILI": ("Liver Safe", False),
        "Caco2_Wang": ("Permeability", True),
    }
    available = [k for k in radar_props if k in df_top.columns]
    if len(available) < 3:
        return

    n_compounds = min(5, len(df_top))
    fig, axes = plt.subplots(1, n_compounds, figsize=(4*n_compounds, 4),
                              subplot_kw=dict(polar=True))
    if n_compounds == 1:
        axes = [axes]

    for idx, (_, row) in enumerate(df_top.head(n_compounds).iterrows()):
        ax = axes[idx]
        values, labels = [], []
        for prop in available:
            label, is_good_high = radar_props[prop]
            val = float(row.get(prop, 0.5))
            if not is_good_high:
                val = 1.0 - val
            values.append(val if pd.notna(val) else 0.5)
            labels.append(label)

        values.append(values[0])
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
    print(f"  Radar chart saved: {fig_path}")


def plot_admet_distributions(df, fig_path):
    """주요 ADMET 속성 분포"""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("ADMET Property Distributions — Screening Hits",
                 fontsize=14, fontweight="bold")
    plot_props = ["HIA_Hou", "Bioavailability_Ma", "Caco2_Wang", "hERG", "AMES", "DILI"]
    colors = ["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728", "#9467bd", "#e377c2"]
    for ax, prop, color in zip(axes.flat, plot_props, colors):
        if prop in df.columns:
            ax.hist(df[prop].dropna(), bins=30, color=color, edgecolor="white", alpha=0.8)
            ax.set_title(ADMET_PROPERTIES.get(prop, (prop,))[0])
            ax.set_xlabel("Score")
            ax.set_ylabel("Count")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Distributions saved: {fig_path}")


def main():
    print("Loading screening hits...")
    df = pd.read_csv(INPUT_FILE)
    print(f"  Input: {len(df)} compounds")

    # ADMET 예측
    print("Predicting ADMET properties...")
    model = ADMETModel()
    admet_results = model.predict(smiles=df["canonical_smiles"].tolist())

    if isinstance(admet_results, pd.DataFrame):
        admet_results = admet_results.reset_index(drop=True)
        df = df.reset_index(drop=True)
        df = pd.concat([df, admet_results], axis=1)

    # ADMET 종합 점수
    df["admet_score"] = compute_admet_score(df)

    # 종합 점수 (ChemBERTa 활성 + ADMET)
    df["combined_score"] = (
        0.5 * df["predicted_active_prob"] +
        0.5 * df["admet_score"]
    )
    df = df.sort_values("combined_score", ascending=False).reset_index(drop=True)

    # 상위 N개 선별
    df_top = df.head(TOP_N).copy()

    # 시각화
    print("Generating figures...")
    plot_admet_heatmap(df_top, os.path.join(FIGURES_DIR, "fig04_admet_heatmap.png"))
    plot_admet_radar(df_top, os.path.join(FIGURES_DIR, "fig05_admet_radar.png"))
    plot_admet_distributions(df, os.path.join(FIGURES_DIR, "fig06_admet_distributions.png"))

    # 저장
    df.to_csv(os.path.join(OUTPUT_DIR, "admet_all_scored.csv"), index=False)
    df_top.to_csv(os.path.join(OUTPUT_DIR, "admet_top_hits.csv"), index=False)

    print(f"\n=== ADMET Screening Summary ===")
    print(f"Input: {len(df)} compounds")
    print(f"Top {TOP_N} selected (combined score >= {df_top['combined_score'].min():.3f})")
    for prop in ["HIA_Hou", "hERG", "AMES", "DILI"]:
        if prop in df.columns:
            print(f"  {ADMET_PROPERTIES[prop][0]}: mean={df[prop].mean():.3f}")
    print(f"\nTop 10:")
    show_cols = ["molecule_chembl_id", "predicted_active_prob", "admet_score", "combined_score"]
    print(df_top.head(10)[show_cols].to_string(index=False))
    print(f"\nOutput: {OUTPUT_DIR}/admet_top_hits.csv")


if __name__ == "__main__":
    main()
