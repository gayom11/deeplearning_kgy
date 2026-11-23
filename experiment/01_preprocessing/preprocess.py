"""
Step 1: 전처리 — Lipinski + PAINS 필터링 + 분자 기술자 계산

학습 데이터와 스크리닝 라이브러리 모두에 동일한 전처리 적용:
1. RDKit으로 SMILES 유효성 검증
2. 분자 기술자 계산 (MW, logP, HBD, HBA, TPSA 등)
3. Lipinski Rule of Five 필터 (1회 위반 허용)
4. PAINS 필터 (pan-assay interference compounds 제거)
5. 기존 BRAF 약물 대비 Tanimoto 유사도 계산
"""
import os
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem, DataStructs
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "00_data_collection")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# 기존 BRAF 약물 (참조용)
REFERENCE_DRUGS = {
    "Vemurafenib":  "CCS(=O)(=O)Nc1ccc(F)c(C(=O)c2cc(-c3ccc(Cl)cc3)c[nH]2)c1F",
    "Dabrafenib":   "CC(C)(C)c1nc(-c2cccc(NS(=O)(=O)c3c(F)cccc3F)c2F)c(-c2ccnc(N)n2)s1",
    "Encorafenib":  "COC(=O)Nc1ccc(Oc2cc(NC(=O)c3ccc(F)cn3)ncn2)c(Cl)c1",
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
    """Lipinski Rule of Five (1회 위반 허용)"""
    violations = 0
    if row["MW"] > 500: violations += 1
    if row["logP"] > 5: violations += 1
    if row["HBD"] > 5: violations += 1
    if row["HBA"] > 10: violations += 1
    return violations <= 1


def pains_filter(smiles_list):
    """PAINS 필터"""
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
    """기존 BRAF 약물 대비 Tanimoto 유사도 (Morgan FP)"""
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


def process_dataset(df, name):
    """하나의 데이터셋에 전처리 적용"""
    print(f"\n--- Processing {name} ({len(df)} compounds) ---")

    # 분자 기술자 계산
    print("  Computing molecular descriptors...")
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
    print(f"  Valid molecules: {len(df)}")

    n_before = len(df)

    # Lipinski 필터
    mask_lip = df.apply(lipinski_filter, axis=1)
    df = df[mask_lip].reset_index(drop=True)
    print(f"  After Lipinski: {len(df)}")

    # PAINS 필터
    mask_pains = pains_filter(df["canonical_smiles"].tolist())
    df = df[mask_pains].reset_index(drop=True)
    print(f"  After PAINS: {len(df)}")

    # 유사도 계산
    print("  Computing similarity to known BRAF drugs...")
    max_sims, best_refs = compute_similarity_to_refs(df["canonical_smiles"].tolist())
    df["max_tanimoto_sim"] = max_sims
    df["most_similar_drug"] = best_refs

    print(f"  {name}: {n_before} → {len(df)} ({(1-len(df)/n_before)*100:.1f}% removed)")
    return df


def plot_preprocessing_summary(df_train, df_screen):
    """전처리 결과 요약 시각화"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Preprocessing Results: Training Data vs Screening Library",
                 fontsize=14, fontweight="bold")

    props = ["MW", "logP", "HBD", "HBA", "TPSA", "max_tanimoto_sim"]
    labels = ["Molecular Weight", "logP", "H-Bond Donors",
              "H-Bond Acceptors", "TPSA (Å²)", "Tanimoto Sim. to BRAF Drugs"]

    for ax, prop, label in zip(axes.flat, props, labels):
        if prop in df_train.columns:
            ax.hist(df_train[prop].dropna(), bins=40, alpha=0.5,
                    color="#1f77b4", label=f"Training ({len(df_train)})")
        if prop in df_screen.columns:
            ax.hist(df_screen[prop].dropna(), bins=40, alpha=0.5,
                    color="#ff7f0e", label=f"Screening ({len(df_screen)})")
        ax.set_xlabel(label)
        ax.set_ylabel("Count")
        ax.legend(fontsize=8)

    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "fig01_preprocessing.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {out}")


def main():
    # 데이터 로드
    df_train = pd.read_csv(os.path.join(DATA_DIR, "braf_training_data.csv"))
    df_screen = pd.read_csv(os.path.join(DATA_DIR, "screening_library.csv"))

    # 전처리 적용
    df_train = process_dataset(df_train, "Training Data")
    df_screen = process_dataset(df_screen, "Screening Library")

    # 시각화
    plot_preprocessing_summary(df_train, df_screen)

    # 저장
    df_train.to_csv(os.path.join(OUTPUT_DIR, "training_filtered.csv"), index=False)
    df_screen.to_csv(os.path.join(OUTPUT_DIR, "screening_filtered.csv"), index=False)

    print(f"\n=== Preprocessing Summary ===")
    print(f"Training:  {len(df_train)} compounds")
    if "activity_label" in df_train.columns:
        n_a = (df_train["activity_label"] == 1).sum()
        n_i = (df_train["activity_label"] == 0).sum()
        print(f"  Active: {n_a}, Inactive: {n_i}, Ambiguous: {len(df_train)-n_a-n_i}")
    print(f"Screening: {len(df_screen)} compounds")
    print(f"Output: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
