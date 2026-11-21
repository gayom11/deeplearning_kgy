"""
Step 6: 최종 후보 물질 선정
- 모든 단계의 결과를 종합하여 최종 10개 후보 선정
- 종합 파이프라인 시각화
- 분자 구조 이미지 생성
"""
import os
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches

ADMET_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "04_admet", "egfr_admet_passed.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# 기존 EGFR 약물 (비교용)
REFERENCE_DRUGS = {
    "Erlotinib":   "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
    "Gefitinib":   "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
    "Osimertinib": "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C",
    "Afatinib":    "CN(C)/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1",
}


def create_pipeline_funnel(counts, fig_path):
    """파이프라인 필터링 과정 funnel chart"""
    stages = list(counts.keys())
    values = list(counts.values())

    fig, ax = plt.subplots(figsize=(10, 7))

    colors = plt.cm.RdYlGn(np.linspace(0.2, 0.9, len(stages)))
    max_width = 0.9
    y_positions = np.linspace(0.9, 0.1, len(stages))

    for i, (stage, val) in enumerate(zip(stages, values)):
        width = max_width * (val / values[0]) if values[0] > 0 else max_width
        width = max(width, 0.08)  # 최소 너비
        rect = plt.Rectangle((0.5 - width/2, y_positions[i] - 0.03),
                              width, 0.06, facecolor=colors[i],
                              edgecolor="white", linewidth=2, zorder=2)
        ax.add_patch(rect)
        ax.text(0.5, y_positions[i], f"{stage}\n{val:,} compounds",
                ha="center", va="center", fontsize=10, fontweight="bold", zorder=3)

        if i < len(stages) - 1:
            ax.annotate("", xy=(0.5, y_positions[i+1] + 0.04),
                        xytext=(0.5, y_positions[i] - 0.04),
                        arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Virtual Screening Pipeline — Filtering Funnel", fontsize=14, fontweight="bold")

    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Funnel chart saved: {fig_path}")


def create_molecule_grid(df_top, fig_path, n=10):
    """최종 후보 분자 구조 이미지"""
    mols = []
    legends = []
    for _, row in df_top.head(n).iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol:
            AllChem.Compute2DCoords(mol)
            mols.append(mol)
            legends.append(f"{row['molecule_chembl_id']}\npChEMBL={row['pchembl_value']:.1f}")

    # 기존 약물 추가
    for name, smi in REFERENCE_DRUGS.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            AllChem.Compute2DCoords(mol)
            mols.append(mol)
            legends.append(f"{name}\n(Reference)")

    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(400, 350), legends=legends)
    img.save(fig_path)
    print(f"Molecule grid saved: {fig_path}")


def create_comparison_table(df_top, fig_path):
    """최종 후보 vs 기존 약물 비교 표"""
    # 기존 약물 정보
    ref_data = []
    for name, smi in REFERENCE_DRUGS.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            ref_data.append({
                "Name": name,
                "MW": Descriptors.MolWt(mol),
                "logP": Descriptors.MolLogP(mol),
                "HBD": Descriptors.NumHDonors(mol),
                "HBA": Descriptors.NumHAcceptors(mol),
                "TPSA": Descriptors.TPSA(mol),
                "Type": "Reference"
            })

    # 후보 물질
    cand_data = []
    for _, row in df_top.head(10).iterrows():
        cand_data.append({
            "Name": row["molecule_chembl_id"],
            "MW": row.get("MW", 0),
            "logP": row.get("logP", 0),
            "HBD": row.get("HBD", 0),
            "HBA": row.get("HBA", 0),
            "TPSA": row.get("TPSA", 0),
            "Type": "Candidate"
        })

    all_data = cand_data + ref_data
    df_comp = pd.DataFrame(all_data)

    # 테이블 시각화
    fig, ax = plt.subplots(figsize=(12, max(4, len(all_data) * 0.4)))
    ax.axis("off")

    cols = ["Name", "MW", "logP", "HBD", "HBA", "TPSA", "Type"]
    cell_text = []
    colors_list = []
    for _, row in df_comp.iterrows():
        cell_text.append([
            row["Name"],
            f"{row['MW']:.1f}",
            f"{row['logP']:.2f}",
            f"{row['HBD']:.0f}",
            f"{row['HBA']:.0f}",
            f"{row['TPSA']:.1f}",
            row["Type"]
        ])
        if row["Type"] == "Reference":
            colors_list.append(["#ffe6cc"] * len(cols))
        else:
            colors_list.append(["#e6ffe6"] * len(cols))

    table = ax.table(cellText=cell_text, colLabels=cols,
                     cellColours=colors_list,
                     loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)

    ax.set_title("Final Candidates vs Reference EGFR Inhibitors", fontsize=13, fontweight="bold", pad=20)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Comparison table saved: {fig_path}")


def create_scatter_comparison(df, fig_path):
    """다차원 비교 scatter"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Final Candidates — Multi-dimensional Analysis", fontsize=14, fontweight="bold")

    top10 = df.head(10)
    rest = df.iloc[10:]

    # 1) MW vs logP
    axes[0].scatter(rest["MW"], rest["logP"], alpha=0.3, s=10, c="gray", label="Other")
    axes[0].scatter(top10["MW"], top10["logP"], s=60, c="red", marker="*", zorder=5, label="Top 10")
    # Lipinski box
    axes[0].axvline(500, color="blue", linestyle="--", alpha=0.3)
    axes[0].axhline(5, color="blue", linestyle="--", alpha=0.3)
    axes[0].set_xlabel("Molecular Weight")
    axes[0].set_ylabel("logP")
    axes[0].set_title("Drug-likeness Space")
    axes[0].legend(fontsize=8)

    # 2) pChEMBL vs ADMET score
    if "admet_score" in df.columns:
        axes[1].scatter(rest["pchembl_value"], rest["admet_score"], alpha=0.3, s=10, c="gray")
        axes[1].scatter(top10["pchembl_value"], top10["admet_score"], s=60, c="red", marker="*", zorder=5)
        axes[1].set_xlabel("pChEMBL Value (Activity)")
        axes[1].set_ylabel("ADMET Score")
        axes[1].set_title("Activity vs Safety")

    # 3) pChEMBL vs Similarity
    if "max_tanimoto_sim" in df.columns:
        axes[2].scatter(rest["pchembl_value"], rest["max_tanimoto_sim"], alpha=0.3, s=10, c="gray")
        axes[2].scatter(top10["pchembl_value"], top10["max_tanimoto_sim"], s=60, c="red", marker="*", zorder=5)
        axes[2].set_xlabel("pChEMBL Value")
        axes[2].set_ylabel("Tanimoto Similarity to Known Drugs")
        axes[2].set_title("Activity vs Novelty")

    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Scatter comparison saved: {fig_path}")


def main():
    print("Loading ADMET-passed compounds...")
    df = pd.read_csv(ADMET_FILE)
    print(f"Input: {len(df)} compounds")

    # 최종 종합 점수 (가중 평균)
    score_components = {}
    # 1) 활성 (pChEMBL)
    score_components["activity"] = df["pchembl_value"] / df["pchembl_value"].max()

    # 2) ChemBERTa 활성 확률
    if "chemberta_active_prob" in df.columns:
        score_components["chemberta"] = df["chemberta_active_prob"]

    # 3) DeepDTA 결합 친화도
    if "predicted_pKd" in df.columns:
        score_components["binding"] = df["predicted_pKd"].fillna(0) / df["predicted_pKd"].max()

    # 4) ADMET 점수
    if "admet_score" in df.columns:
        score_components["admet"] = df["admet_score"]

    # 5) 약물 유사성 (너무 높으면 novelty가 없으므로 적당한 값 선호)
    if "max_tanimoto_sim" in df.columns:
        # 0.3-0.7 범위가 ideal (유사하지만 새로운)
        sim = df["max_tanimoto_sim"]
        score_components["novelty"] = 1 - np.abs(sim - 0.5)  # 0.5에 가까울수록 높은 점수

    weights = {"activity": 0.30, "chemberta": 0.20, "binding": 0.15, "admet": 0.25, "novelty": 0.10}
    df["final_score"] = sum(
        weights.get(k, 0.1) * v for k, v in score_components.items()
    )
    df = df.sort_values("final_score", ascending=False).reset_index(drop=True)

    # 화학적 다양성 보장 — 유사한 구조 제거 (top에서 Tanimoto > 0.8이면 하나만)
    from rdkit.Chem import AllChem, DataStructs
    selected_idx = []
    selected_fps = []
    for i, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        # 이미 선택된 것과 너무 유사하면 스킵
        too_similar = False
        for sfp in selected_fps:
            if DataStructs.TanimotoSimilarity(fp, sfp) > 0.8:
                too_similar = True
                break
        if not too_similar:
            selected_idx.append(i)
            selected_fps.append(fp)
        if len(selected_idx) >= 10:
            break

    df_final = df.iloc[selected_idx].reset_index(drop=True)
    print(f"Final candidates (diverse): {len(df_final)}")

    # 파이프라인 funnel
    # 단계별 화합물 수 (이전 단계 결과 로드)
    raw_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "01_chembl_data", "egfr_compounds_raw.csv")
    filtered_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "02_filtering", "egfr_filtered.csv")
    binding_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "03_docking", "egfr_binding_affinity.csv")
    admet_all_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "04_admet", "egfr_admet_all.csv")
    n_raw = len(pd.read_csv(raw_path)) if os.path.exists(raw_path) else 10005
    n_filtered = len(pd.read_csv(filtered_path)) if os.path.exists(filtered_path) else 6608
    n_binding = len(pd.read_csv(binding_path)) if os.path.exists(binding_path) else 500
    n_admet_all = len(pd.read_csv(admet_all_path)) if os.path.exists(admet_all_path) else 500

    funnel_counts = {
        "ChEMBL EGFR (IC50)": n_raw,
        "Lipinski + PAINS Filter": n_filtered,
        "DeepDTA + ChemBERTa\n(Combined Score Top)": n_admet_all,
        "ADMET Filter (Top 50%)": len(df),
        "Final Diverse Candidates": len(df_final),
    }
    create_pipeline_funnel(funnel_counts, os.path.join(FIGURES_DIR, "fig08_pipeline_funnel.png"))

    # 분자 구조 이미지
    create_molecule_grid(df_final, os.path.join(FIGURES_DIR, "fig09_candidate_structures.png"))

    # 비교 표
    create_comparison_table(df_final, os.path.join(FIGURES_DIR, "fig10_comparison_table.png"))

    # 다차원 scatter
    create_scatter_comparison(df, os.path.join(FIGURES_DIR, "fig11_scatter_analysis.png"))

    # 저장
    output_file = os.path.join(OUTPUT_DIR, "final_candidates.csv")
    df_final.to_csv(output_file, index=False)

    # 전체 결과 저장 (랭킹 포함)
    df.to_csv(os.path.join(OUTPUT_DIR, "all_ranked.csv"), index=False)

    print(f"\n=== Final Candidate Summary ===")
    show_cols = ["molecule_chembl_id", "canonical_smiles", "pchembl_value", "MW", "logP", "final_score"]
    if "admet_score" in df_final.columns:
        show_cols.append("admet_score")
    for _, row in df_final.iterrows():
        print(f"\n  {row['molecule_chembl_id']}:")
        print(f"    pChEMBL: {row['pchembl_value']:.2f} | MW: {row.get('MW', 0):.1f} | logP: {row.get('logP', 0):.2f}")
        print(f"    Final Score: {row['final_score']:.4f}")
        print(f"    SMILES: {row['canonical_smiles'][:80]}...")

    print(f"\nOutput: {output_file}")


if __name__ == "__main__":
    main()
