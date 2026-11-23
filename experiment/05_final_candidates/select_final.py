"""
Step 5: 최종 후보 선정 — 개별 기준 통과 방식

각 모델의 결과를 하나의 점수로 임의 통합하지 않는다.
대신 모든 기준을 개별적으로 통과한 화합물만 최종 후보로 선정한다.

기준:
1. ChemBERTa 활성 확률 >= 0.8 (BRAF 활성 예측)
2. ADMET 주요 항목 통과 (HIA >= 0.5, AMES <= 0.5)
3. 용해도 logS >= -4
4. GNN pIC50 >= 6.0 (IC50 <= 1μM)
5. GNINA CNN affinity 상위 (도킹 성공)
6. Tanimoto < 0.7 (기존 BRAF 약물과 구조 차별성)
"""
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.DataStructs import TanimotoSimilarity
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
EXPERIMENT = os.path.join(BASE, "..")
FIGURES_DIR = os.path.join(EXPERIMENT, "figures")
DOCKING_DIR = os.path.join(EXPERIMENT, "03_docking")

REFERENCE_DRUGS = {
    "Vemurafenib":  "CCS(=O)(=O)Nc1ccc(F)c(C(=O)c2cc(-c3ccc(Cl)cc3)c[nH]2)c1F",
    "Dabrafenib":   "CC(C)(C)c1nc(-c2cccc(NS(=O)(=O)c3c(F)cccc3F)c2F)c(-c2ccnc(N)n2)s1",
    "Encorafenib":  "COC(=O)Nc1ccc(Oc2cc(NC(=O)c3ccc(F)cn3)ncn2)c(Cl)c1",
}


def parse_docked_sdf(sdf_path):
    """GNINA docked.sdf 파싱"""
    from rdkit.Chem import SDMolSupplier
    suppl = SDMolSupplier(sdf_path)
    records = []
    for mol in suppl:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        records.append({
            "molecule_chembl_id": name,
            "minimizedAffinity": float(props.get("minimizedAffinity", 0)),
            "CNNaffinity": float(props.get("CNNaffinity", 0)),
            "CNNscore": float(props.get("CNNscore", 0)),
        })
    df = pd.DataFrame(records)
    # 각 화합물 최고 CNN affinity만
    df = df.sort_values("CNNaffinity", ascending=False).drop_duplicates("molecule_chembl_id", keep="first")
    return df


def compute_tanimoto_to_refs(smiles):
    """기존 BRAF 약물 대비 최대 Tanimoto"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 1.0
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    max_sim = 0
    for name, smi in REFERENCE_DRUGS.items():
        ref_mol = Chem.MolFromSmiles(smi)
        if ref_mol:
            ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
            sim = TanimotoSimilarity(fp, ref_fp)
            max_sim = max(max_sim, sim)
    return max_sim


def main():
    # === 1. 데이터 로드 ===
    # 용해도
    sol = pd.read_csv(os.path.join(EXPERIMENT, "04_admet", "solubility_scored.csv"))
    # GNN pIC50
    act = pd.read_csv(os.path.join(EXPERIMENT, "04_admet", "activity_scored.csv"))
    # GNINA 도킹
    docked_path = os.path.join(DOCKING_DIR, "docked.sdf")
    dock = parse_docked_sdf(docked_path)

    print(f"=== Data Loaded ===")
    print(f"Solubility scored: {len(sol)}")
    print(f"Activity scored:   {len(act)}")
    print(f"Docked:            {len(dock)}")

    # === 2. 병합 ===
    # sol과 act는 같은 300개 기반이므로 merge
    df = sol.merge(act[["molecule_chembl_id", "pred_pIC50", "activity_potent"]],
                   on="molecule_chembl_id", how="left")
    df = df.merge(dock[["molecule_chembl_id", "CNNaffinity", "minimizedAffinity", "CNNscore"]],
                  on="molecule_chembl_id", how="left")

    print(f"Merged: {len(df)} compounds")
    print(f"  with docking: {df['CNNaffinity'].notna().sum()}")

    # === 3. 개별 기준 판정 ===
    df["pass_chemberta"] = df["predicted_active_prob"] >= 0.8
    df["pass_solubility"] = df["pred_logS"] >= -4
    df["pass_pIC50"] = df["pred_pIC50"] >= 6.0
    df["pass_docking"] = df["CNNaffinity"].notna() & (df["CNNaffinity"] > 5.0)

    # ADMET 개별 항목
    admet_pass = True
    if "HIA_Hou" in df.columns:
        df["pass_HIA"] = df["HIA_Hou"] >= 0.5
        admet_pass = admet_pass & df["pass_HIA"]
    if "AMES" in df.columns:
        df["pass_AMES"] = df["AMES"] <= 0.5  # AMES 양성 = 유해
        admet_pass = admet_pass & df["pass_AMES"]
    df["pass_admet"] = admet_pass

    # 전체 통과
    df["pass_all"] = (df["pass_chemberta"] & df["pass_solubility"] &
                      df["pass_pIC50"] & df["pass_docking"] & df["pass_admet"])

    print(f"\n=== Individual Criteria ===")
    for col in [c for c in df.columns if c.startswith("pass_")]:
        print(f"  {col}: {df[col].sum()}")

    # === 4. Tanimoto 다양성 필터 ===
    candidates = df[df["pass_all"]].copy()
    candidates["max_sim_to_drugs"] = candidates["canonical_smiles"].apply(compute_tanimoto_to_refs)
    candidates = candidates[candidates["max_sim_to_drugs"] < 0.7]

    # CNN affinity 순으로 정렬
    candidates = candidates.sort_values("CNNaffinity", ascending=False).head(10)
    print(f"\nFinal candidates (all pass + Tanimoto < 0.7): {len(candidates)}")

    # === 5. 저장 ===
    save_cols = ["molecule_chembl_id", "canonical_smiles", "predicted_active_prob",
                 "pred_pIC50", "pred_logS", "CNNaffinity", "minimizedAffinity",
                 "admet_score", "max_sim_to_drugs", "MW", "logP"]
    save_cols = [c for c in save_cols if c in candidates.columns]
    candidates[save_cols].to_csv(os.path.join(BASE, "final_candidates.csv"), index=False)
    df.to_csv(os.path.join(BASE, "all_ranked.csv"), index=False)

    print(f"\nTop 10:")
    for _, r in candidates.iterrows():
        print(f"  {r['molecule_chembl_id']} | ChemBERTa={r['predicted_active_prob']:.3f} "
              f"| pIC50={r['pred_pIC50']:.2f} | logS={r['pred_logS']:.2f} "
              f"| CNN={r['CNNaffinity']:.2f} | Vina={r['minimizedAffinity']:.2f}")

    # === 6. 시각화 ===
    # 6a. 후보 구조 그림
    mols = []
    legends = []
    for _, r in candidates.iterrows():
        mol = Chem.MolFromSmiles(r["canonical_smiles"])
        if mol:
            AllChem.Compute2DCoords(mol)
            mols.append(mol)
            legends.append(r["molecule_chembl_id"])

    if mols:
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 250), legends=legends)
        img.save(os.path.join(FIGURES_DIR, "fig08_candidate_structures.png"))

    # 6b. 개별 모델 결과 scatter
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Final Candidates: Individual Model Results", fontsize=14, fontweight="bold")

    # ChemBERTa vs pIC50
    ax = axes[0, 0]
    ax.scatter(df["predicted_active_prob"], df["pred_pIC50"], alpha=0.3, s=10, c="gray")
    if len(candidates) > 0:
        ax.scatter(candidates["predicted_active_prob"], candidates["pred_pIC50"],
                   s=40, c="red", zorder=5, label=f"Final ({len(candidates)})")
    ax.axvline(0.8, color="blue", linestyle="--", alpha=0.5)
    ax.axhline(6.0, color="blue", linestyle="--", alpha=0.5)
    ax.set_xlabel("ChemBERTa Active Prob")
    ax.set_ylabel("Predicted pIC50")
    ax.set_title("Activity: ChemBERTa vs GNN")
    ax.legend()

    # logS vs CNN
    ax = axes[0, 1]
    docked_df = df[df["CNNaffinity"].notna()]
    ax.scatter(docked_df["pred_logS"], docked_df["CNNaffinity"], alpha=0.3, s=10, c="gray")
    if len(candidates) > 0:
        ax.scatter(candidates["pred_logS"], candidates["CNNaffinity"],
                   s=40, c="red", zorder=5, label=f"Final ({len(candidates)})")
    ax.axvline(-4, color="blue", linestyle="--", alpha=0.5)
    ax.set_xlabel("Predicted logS")
    ax.set_ylabel("CNN Affinity")
    ax.set_title("Solubility vs Docking")
    ax.legend()

    # MW vs logP
    ax = axes[1, 0]
    ax.scatter(df["MW"], df["logP"], alpha=0.3, s=10, c="gray")
    if len(candidates) > 0:
        ax.scatter(candidates["MW"], candidates["logP"],
                   s=40, c="red", zorder=5, label=f"Final ({len(candidates)})")
    ax.set_xlabel("Molecular Weight")
    ax.set_ylabel("logP")
    ax.set_title("Drug-likeness Space")
    ax.legend()

    # Pass counts bar
    ax = axes[1, 1]
    criteria = ["pass_chemberta", "pass_pIC50", "pass_solubility", "pass_docking", "pass_admet", "pass_all"]
    labels = ["ChemBERTa\n≥0.8", "pIC50\n≥6.0", "logS\n≥-4", "Docking\n>5.0", "ADMET", "ALL\nPass"]
    counts = [df[c].sum() for c in criteria]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#e377c2"]
    ax.bar(labels, counts, color=colors, edgecolor="white")
    for i, v in enumerate(counts):
        ax.text(i, v + 2, str(v), ha="center", fontsize=10, fontweight="bold")
    ax.set_ylabel("Compounds Passing")
    ax.set_title("Individual Criteria Pass Counts")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig09_scatter_analysis.png"), dpi=150, bbox_inches="tight")
    plt.close()

    print(f"\nFigures saved to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
