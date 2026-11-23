"""
Step 1-B: Scaffold 기반 대표 화합물 선별

Lipinski+PAINS 통과 화합물을 Bemis-Murcko scaffold로 묶고,
그룹당 최대 3개만 남긴다. → 한 계열에 후보가 몰리지 않도록.

참고: 교수 가이드 Page 8 — "공통 골격(scaffold)이 같은 후보끼리
묶고 그룹당 2~3개만 남긴다"
"""
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "screening_filtered.csv")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")

MAX_PER_SCAFFOLD = 3


def get_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaf)
    except Exception:
        return None


def main():
    df = pd.read_csv(INPUT_FILE)
    print(f"Input: {len(df)} compounds")

    # Scaffold 추출
    print("Extracting Bemis-Murcko scaffolds...")
    df["scaffold"] = df["canonical_smiles"].apply(get_scaffold)
    df_valid = df[df["scaffold"].notna()].copy()
    print(f"Valid scaffolds: {len(df_valid)}")

    n_scaffolds = df_valid["scaffold"].nunique()
    print(f"Unique scaffolds: {n_scaffolds}")

    # 그룹당 최대 3개 선별
    # 선택 규칙: Ro5 위반 수 적은 순 → MW 작은 순 → ChEMBL ID 순
    # (약물 가능성 판단은 예측·도킹 단계에서 수행)
    SELECTION_RULE = "Ro5 violations ASC → MW ASC → ID ASC"
    print(f"Selecting max {MAX_PER_SCAFFOLD} per scaffold...")
    print(f"  Selection rule: {SELECTION_RULE}")
    selected = []
    excluded_scaffold = []

    for scaf, group in df_valid.groupby("scaffold"):
        if len(group) <= MAX_PER_SCAFFOLD:
            for _, row in group.iterrows():
                selected.append({**row, "selection_rule": "group_size_leq_3"})
        else:
            # Ro5 위반 수 적은 순 → MW 작은 순 → ID 순
            group_sorted = group.sort_values(
                ["HBD", "MW", "molecule_chembl_id"],  # HBD as proxy for violations
                ascending=[True, True, True]
            )
            # 실제 위반 수 계산
            def count_violations(r):
                v = 0
                if r["MW"] > 500: v += 1
                if r["logP"] > 5: v += 1
                if r["HBD"] > 5: v += 1
                if r["HBA"] > 10: v += 1
                return v
            group_sorted["ro5_violations"] = group_sorted.apply(count_violations, axis=1)
            group_sorted = group_sorted.sort_values(
                ["ro5_violations", "MW", "molecule_chembl_id"],
                ascending=[True, True, True]
            )
            chosen = group_sorted.head(MAX_PER_SCAFFOLD)
            for _, row in chosen.iterrows():
                selected.append({**row, "selection_rule": SELECTION_RULE})
            for _, row in group_sorted.iloc[MAX_PER_SCAFFOLD:].iterrows():
                excluded_scaffold.append({
                    "molecule_chembl_id": row["molecule_chembl_id"],
                    "scaffold": scaf,
                    "group_size": len(group),
                    "ro5_violations": row.get("ro5_violations", ""),
                    "reason": f"scaffold_group_limit: kept {MAX_PER_SCAFFOLD}/{len(group)} by {SELECTION_RULE}"
                })

    df_selected = pd.DataFrame(selected).reset_index(drop=True)
    print(f"After scaffold selection: {len(df_selected)}")

    # 시각화
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Scaffold-based Selection", fontsize=14, fontweight="bold")

    # 1) scaffold 크기 분포
    scaf_sizes = df_valid.groupby("scaffold").size()
    axes[0].hist(scaf_sizes, bins=50, color="#1f77b4", edgecolor="white", alpha=0.8)
    axes[0].set_xlabel("Compounds per Scaffold")
    axes[0].set_ylabel("Number of Scaffolds")
    axes[0].set_title(f"Scaffold Group Sizes ({n_scaffolds} scaffolds)")
    axes[0].axvline(MAX_PER_SCAFFOLD, color="red", linestyle="--", label=f"Max {MAX_PER_SCAFFOLD}")
    axes[0].legend()

    # 2) 선별 전후 MW 분포
    axes[1].hist(df_valid["MW"], bins=40, alpha=0.5, color="#d62728", label=f"Before ({len(df_valid)})")
    axes[1].hist(df_selected["MW"], bins=40, alpha=0.5, color="#2ca02c", label=f"After ({len(df_selected)})")
    axes[1].set_xlabel("Molecular Weight")
    axes[1].set_ylabel("Count")
    axes[1].set_title("MW Distribution")
    axes[1].legend()

    # 3) 상위 10 scaffold
    top_scaf = scaf_sizes.nlargest(10)
    axes[2].barh(range(len(top_scaf)), top_scaf.values, color="#ff7f0e")
    axes[2].set_yticks(range(len(top_scaf)))
    axes[2].set_yticklabels([s[:30]+"..." if len(s)>30 else s for s in top_scaf.index], fontsize=7)
    axes[2].set_xlabel("Compound Count")
    axes[2].set_title("Top 10 Scaffolds by Size")
    axes[2].invert_yaxis()

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig01b_scaffold_selection.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # 저장
    df_selected.to_csv(os.path.join(OUTPUT_DIR, "screening_scaffold_selected.csv"), index=False)
    pd.DataFrame(excluded_scaffold).to_csv(
        os.path.join(OUTPUT_DIR, "scaffold_excluded.csv"), index=False)

    print(f"\n=== Scaffold Selection Summary ===")
    print(f"Before: {len(df_valid)}")
    print(f"Unique scaffolds: {n_scaffolds}")
    print(f"After (max {MAX_PER_SCAFFOLD}/group): {len(df_selected)}")
    print(f"Excluded: {len(excluded_scaffold)}")
    print(f"Reduction: {(1-len(df_selected)/len(df_valid))*100:.1f}%")


if __name__ == "__main__":
    main()
