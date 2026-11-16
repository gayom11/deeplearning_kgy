"""
Step 2: Lipinski Rule of Five + 추가 Drug-likeness 필터링
- MW <= 500, logP <= 5, HBD <= 5, HBA <= 10
- PAINS 필터 (pan-assay interference compounds 제거)
- 분자 유사도 계산 (Morgan Fingerprint)
"""
import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Draw, AllChem, DataStructs
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "01_chembl_data")
INPUT_FILE = os.path.join(INPUT_DIR, "egfr_compounds_raw.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# 기존 EGFR 약물 (비교 기준)
REFERENCE_DRUGS = {
    "Erlotinib":   "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
    "Gefitinib":   "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
    "Osimertinib": "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C",
    "Afatinib":    "CN(C)/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1",
}


def compute_descriptors(smiles):
    """SMILES로부터 분자 기술자 계산"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "MW": Descriptors.MolWt(mol),
        "logP": Descriptors.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotBonds": Lipinski.NumRotatableBonds(mol),
        "RingCount": Lipinski.RingCount(mol),
        "AromaticRings": Descriptors.NumAromaticRings(mol),
        "HeavyAtoms": Descriptors.HeavyAtomCount(mol),
    }


def lipinski_filter(row):
    """Lipinski Rule of Five 통과 여부"""
    violations = 0
    if row["MW"] > 500: violations += 1
    if row["logP"] > 5: violations += 1
    if row["HBD"] > 5: violations += 1
    if row["HBA"] > 10: violations += 1
    return violations <= 1  # 1개까지 위반 허용


def pains_filter(smiles_list):
    """PAINS 필터 — pan-assay interference compounds 제거"""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    is_clean = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            is_clean.append(False)
        else:
            is_clean.append(not catalog.HasMatch(mol))
    return is_clean


def compute_similarity_to_refs(smiles_list):
    """기존 EGFR 약물 대비 Tanimoto 유사도 (Morgan FP, radius=2)"""
    ref_fps = {}
    for name, smi in REFERENCE_DRUGS.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            ref_fps[name] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

    max_sims = []
    best_refs = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            max_sims.append(0.0)
            best_refs.append("")
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        sims = {name: DataStructs.TanimotoSimilarity(fp, rfp) for name, rfp in ref_fps.items()}
        best = max(sims, key=sims.get)
        max_sims.append(sims[best])
        best_refs.append(best)
    return max_sims, best_refs


def plot_filtering_summary(df_before, df_after):
    """필터링 전후 비교 시각화"""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("Lipinski & PAINS Filtering Results", fontsize=14, fontweight="bold")

    props = ["MW", "logP", "HBD", "HBA", "TPSA", "pchembl_value"]
    labels = ["Molecular Weight", "logP", "H-Bond Donors", "H-Bond Acceptors", "TPSA (Å²)", "pChEMBL Value"]
    colors_before = "#d62728"
    colors_after = "#2ca02c"

    for ax, prop, label in zip(axes.flat, props, labels):
        if prop in df_before.columns:
            ax.hist(df_before[prop].dropna(), bins=40, alpha=0.5, color=colors_before, label=f"Before ({len(df_before)})")
        if prop in df_after.columns:
            ax.hist(df_after[prop].dropna(), bins=40, alpha=0.7, color=colors_after, label=f"After ({len(df_after)})")
        ax.set_xlabel(label)
        ax.set_ylabel("Count")
        ax.legend(fontsize=8)

    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "fig01_lipinski_filtering.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {out}")


def plot_similarity_distribution(df):
    """유사도 분포 시각화"""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["max_tanimoto_sim"], bins=50, color="#1f77b4", edgecolor="white", alpha=0.8)
    ax.axvline(x=0.4, color="red", linestyle="--", linewidth=1.5, label="Similarity threshold (0.4)")
    ax.set_xlabel("Max Tanimoto Similarity to Known EGFR Drugs")
    ax.set_ylabel("Count")
    ax.set_title("Molecular Similarity to Reference EGFR Inhibitors")
    ax.legend()
    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "fig02_similarity_distribution.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {out}")


def main():
    print("Loading ChEMBL data...")
    df = pd.read_csv(INPUT_FILE)
    print(f"Input compounds: {len(df)}")

    # 분자 기술자 계산
    print("Computing molecular descriptors...")
    desc_list = []
    valid_idx = []
    for i, smi in enumerate(df["canonical_smiles"]):
        desc = compute_descriptors(smi)
        if desc is not None:
            desc_list.append(desc)
            valid_idx.append(i)

    df = df.iloc[valid_idx].reset_index(drop=True)
    desc_df = pd.DataFrame(desc_list)
    df = pd.concat([df, desc_df], axis=1)
    print(f"Valid molecules: {len(df)}")

    df_before = df.copy()

    # Lipinski 필터
    print("Applying Lipinski Rule of Five...")
    mask_lipinski = df.apply(lipinski_filter, axis=1)
    df = df[mask_lipinski].reset_index(drop=True)
    print(f"After Lipinski filter: {len(df)}")

    # PAINS 필터
    print("Applying PAINS filter...")
    mask_pains = pains_filter(df["canonical_smiles"].tolist())
    df = df[mask_pains].reset_index(drop=True)
    print(f"After PAINS filter: {len(df)}")

    # 유사도 계산
    print("Computing similarity to known EGFR drugs...")
    max_sims, best_refs = compute_similarity_to_refs(df["canonical_smiles"].tolist())
    df["max_tanimoto_sim"] = max_sims
    df["most_similar_drug"] = best_refs

    # 시각화
    plot_filtering_summary(df_before, df)
    plot_similarity_distribution(df)

    # 저장
    output_file = os.path.join(OUTPUT_DIR, "egfr_filtered.csv")
    df.to_csv(output_file, index=False)

    print(f"\n=== Filtering Summary ===")
    print(f"Before filtering: {len(df_before)}")
    print(f"After Lipinski + PAINS: {len(df)}")
    print(f"Removal rate: {(1 - len(df)/len(df_before))*100:.1f}%")
    print(f"Top similar compounds:")
    top = df.nlargest(10, "max_tanimoto_sim")[["molecule_chembl_id", "pchembl_value", "MW", "logP", "max_tanimoto_sim", "most_similar_drug"]]
    print(top.to_string(index=False))
    print(f"\nOutput: {output_file}")


if __name__ == "__main__":
    main()
