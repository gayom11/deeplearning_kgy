"""
Step 7: 분자 구조 생성 (REINVENT-inspired)
- Murcko Scaffold 분해 → 스캐폴드 기반 아날로그 생성
- SELFIES 기반 랜덤 돌연변이 → 새로운 화학 공간 탐색
- RDKit 기반 bioisosteric replacement
- 생성된 분자를 기존 파이프라인으로 재평가
"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import random
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, Draw, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
import selfies as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CANDIDATES_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "final_candidates.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")

random.seed(42)
np.random.seed(42)

# SELFIES 토큰 (일반적인 약물 화학 그룹)
COMMON_TOKENS = [
    "[C]", "[=C]", "[N]", "[=N]", "[O]", "[=O]", "[S]", "[F]", "[Cl]",
    "[Br]", "[Ring1]", "[Branch1]", "[Branch2]", "[#C]", "[NH]", "[OH]",
    "[C@@H1]", "[C@H1]", "[Ring2]", "[=N]", "[P]",
]

# Bioisosteric 치환 그룹
BIOISOSTERES = {
    "F": ["Cl", "Br", "CF3"],
    "Cl": ["F", "Br", "CF3"],
    "OCH3": ["OH", "OCHF2", "SCH3"],
    "NH2": ["NHCH3", "OH", "NHCOCH3"],
    "COOH": ["CONH2", "SO2NH2", "CONHOH"],
    "OH": ["NH2", "SH", "F"],
}


def get_murcko_scaffold(smiles):
    """Murcko 스캐폴드 추출"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return None


def selfies_mutate(smiles, n_mutations=3, n_attempts=50):
    """SELFIES 기반 분자 돌연변이"""
    generated = set()
    try:
        selfies_str = sf.encoder(smiles)
        if selfies_str is None:
            return generated
    except Exception:
        return generated

    tokens = list(sf.split_selfies(selfies_str))

    for _ in range(n_attempts):
        mutated = tokens.copy()

        for _ in range(random.randint(1, n_mutations)):
            mutation_type = random.choice(["replace", "insert", "delete"])

            if mutation_type == "replace" and len(mutated) > 0:
                pos = random.randint(0, len(mutated) - 1)
                mutated[pos] = random.choice(COMMON_TOKENS)
            elif mutation_type == "insert":
                pos = random.randint(0, len(mutated))
                mutated.insert(pos, random.choice(COMMON_TOKENS))
            elif mutation_type == "delete" and len(mutated) > 2:
                pos = random.randint(0, len(mutated) - 1)
                mutated.pop(pos)

        new_selfies = "".join(mutated)
        try:
            new_smiles = sf.decoder(new_selfies)
            if new_smiles is None:
                continue
            mol = Chem.MolFromSmiles(new_smiles)
            if mol is None:
                continue
            # Validity checks
            canon = Chem.MolToSmiles(mol)
            mw = Descriptors.MolWt(mol)
            if 200 <= mw <= 600 and canon != smiles:  # Drug-like MW range
                generated.add(canon)
        except Exception:
            continue

    return generated


def scaffold_hop(smiles, n_variants=20):
    """스캐폴드 기반 변이 — 사이드체인 변경"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return set()

    generated = set()

    # R-group 위치에서 작은 변경
    # Halogen 교체
    for old, replacements in [("F", ["Cl", "Br"]), ("Cl", ["F", "Br"]), ("Br", ["F", "Cl"])]:
        for new in replacements:
            new_smi = smiles.replace(old, new, 1)
            m = Chem.MolFromSmiles(new_smi)
            if m:
                generated.add(Chem.MolToSmiles(m))

    # 메틸/에틸 추가/제거 시도
    patterns = [
        ("(C)", ""),        # 메틸 제거
        ("(CC)", "(C)"),    # 에틸→메틸
        ("(OC)", "(O)"),    # 메톡시→하이드록시
        ("(F)", "(Cl)"),    # F→Cl
    ]
    for old_pat, new_pat in patterns:
        if old_pat in smiles:
            new_smi = smiles.replace(old_pat, new_pat, 1)
            m = Chem.MolFromSmiles(new_smi)
            if m:
                canon = Chem.MolToSmiles(m)
                mw = Descriptors.MolWt(m)
                if 200 <= mw <= 600:
                    generated.add(canon)

    return generated


def evaluate_generated(generated_smiles, ref_smiles):
    """생성된 분자의 Drug-likeness 평가"""
    results = []
    ref_fps = []
    for smi in ref_smiles:
        m = Chem.MolFromSmiles(smi)
        if m:
            ref_fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))

    for smi in generated_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)

        # Lipinski check
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        if violations > 1:
            continue

        # Novelty (vs original candidates)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        max_sim = max(DataStructs.TanimotoSimilarity(fp, rfp) for rfp in ref_fps) if ref_fps else 0

        scaffold = get_murcko_scaffold(smi)

        results.append({
            "smiles": smi,
            "MW": mw,
            "logP": logp,
            "HBD": hbd,
            "HBA": hba,
            "TPSA": tpsa,
            "max_sim_to_candidates": max_sim,
            "scaffold": scaffold or "",
            "source": "generated",
        })

    return pd.DataFrame(results)


def main():
    print("Loading final candidates...")
    df = pd.read_csv(CANDIDATES_FILE)
    print(f"Seed molecules: {len(df)}")

    seed_smiles = df["canonical_smiles"].tolist()

    # 스캐폴드 분석
    scaffolds = {}
    for smi in seed_smiles:
        scaffold = get_murcko_scaffold(smi)
        if scaffold:
            scaffolds[scaffold] = scaffolds.get(scaffold, 0) + 1

    print(f"Unique scaffolds: {len(scaffolds)}")
    for s, c in sorted(scaffolds.items(), key=lambda x: -x[1])[:5]:
        print(f"  {s}: {c}")

    # === 분자 생성 ===
    all_generated = set()

    # 1) SELFIES 돌연변이
    print("\n1) SELFIES-based mutation...")
    for smi in seed_smiles:
        mutants = selfies_mutate(smi, n_mutations=3, n_attempts=100)
        all_generated.update(mutants)
    print(f"   Generated: {len(all_generated)} unique molecules")

    # 2) 스캐폴드 호핑
    print("2) Scaffold hopping...")
    n_before = len(all_generated)
    for smi in seed_smiles:
        variants = scaffold_hop(smi, n_variants=20)
        all_generated.update(variants)
    print(f"   Added: {len(all_generated) - n_before} variants")

    # 중복 및 원본 제거
    all_generated -= set(seed_smiles)
    print(f"\nTotal unique generated molecules: {len(all_generated)}")

    # 평가
    print("Evaluating generated molecules...")
    df_gen = evaluate_generated(list(all_generated), seed_smiles)
    print(f"Drug-like generated molecules: {len(df_gen)}")

    # 유사도 기반 분류
    # Novel: sim < 0.5, Analog: 0.5-0.8, Very similar: > 0.8
    if len(df_gen) > 0:
        df_gen["category"] = pd.cut(
            df_gen["max_sim_to_candidates"],
            bins=[0, 0.4, 0.6, 0.8, 1.0],
            labels=["Novel", "Moderate", "Similar", "Very Similar"]
        )

    # === 시각화 ===

    # 1) 생성 분자 vs 원본 화학 공간 비교
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Generated Molecules — Chemical Space Exploration", fontsize=14, fontweight="bold")

    if len(df_gen) > 0:
        # MW vs logP
        axes[0].scatter(df_gen["MW"], df_gen["logP"], alpha=0.3, s=15, c="#2ca02c", label=f"Generated ({len(df_gen)})")
    for _, row in df.iterrows():
        axes[0].scatter(row.get("MW", 400), row.get("logP", 3), s=100, c="red", marker="*", zorder=5)
    axes[0].scatter([], [], s=100, c="red", marker="*", label="Seed candidates")
    axes[0].set_xlabel("Molecular Weight")
    axes[0].set_ylabel("logP")
    axes[0].set_title("Chemical Space")
    axes[0].legend(fontsize=8)

    # Similarity 분포
    if len(df_gen) > 0:
        axes[1].hist(df_gen["max_sim_to_candidates"], bins=30, color="#ff7f0e", edgecolor="white", alpha=0.8)
    axes[1].set_xlabel("Tanimoto Similarity to Seed Candidates")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Novelty Distribution")

    # Category 파이차트
    if len(df_gen) > 0 and "category" in df_gen.columns:
        counts = df_gen["category"].value_counts()
        if len(counts) > 0:
            axes[2].pie(counts.values, labels=counts.index, autopct="%1.0f%%",
                       colors=["#2ca02c", "#ff7f0e", "#1f77b4", "#d62728"])
            axes[2].set_title("Novelty Categories")

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig12_generated_molecules.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # 2) 상위 생성 분자 구조 이미지
    if len(df_gen) > 0:
        # Novel이면서 drug-like인 것 우선
        df_gen_sorted = df_gen.sort_values("max_sim_to_candidates", ascending=True)
        top_novel = df_gen_sorted.head(min(10, len(df_gen_sorted)))

        mols = []
        legends = []
        for _, row in top_novel.iterrows():
            m = Chem.MolFromSmiles(row["smiles"])
            if m:
                AllChem.Compute2DCoords(m)
                mols.append(m)
                legends.append(f"Sim={row['max_sim_to_candidates']:.2f}\nMW={row['MW']:.0f}")

        if mols:
            img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(350, 300), legends=legends)
            img.save(os.path.join(FIGURES_DIR, "fig13_novel_generated.png"))
            print(f"Novel molecule grid saved")

    # 저장
    output_file = os.path.join(OUTPUT_DIR, "generated_molecules.csv")
    if len(df_gen) > 0:
        df_gen.to_csv(output_file, index=False)

    print(f"\n=== Molecular Generation Summary ===")
    print(f"Seed molecules: {len(seed_smiles)}")
    print(f"Unique scaffolds: {len(scaffolds)}")
    print(f"Total generated: {len(all_generated)}")
    print(f"Drug-like generated: {len(df_gen)}")
    if len(df_gen) > 0 and "category" in df_gen.columns:
        print(f"Novelty breakdown:")
        print(df_gen["category"].value_counts().to_string())
    print(f"\nOutput: {output_file}")


if __name__ == "__main__":
    main()
