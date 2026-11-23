"""
Step 5+6: 최종 후보 선정 + 분자 생성

1) GNINA 도킹 결과 파싱 → ADMET 데이터와 병합
2) 종합 점수 (ChemBERTa + ADMET + GNINA) → 최종 10개 선정
3) SELFIES 기반 분자 생성 → 718+ 신규 분자
4) 생성 분자 ADMET 재평가
5) 파이프라인 funnel + 분자 구조 + 비교 시각화
"""
import os
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

ADMET_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "04_admet", "admet_top_hits.csv")
DOCKED_SDF = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "03_docking", "docked.sdf")
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

random.seed(42)
np.random.seed(42)

REFERENCE_DRUGS = {
    "Vemurafenib":  "CCS(=O)(=O)Nc1ccc(F)c(C(=O)c2cc(-c3ccc(Cl)cc3)c[nH]2)c1F",
    "Dabrafenib":   "CC(C)(C)c1nc(-c2cccc(NS(=O)(=O)c3c(F)cccc3F)c2F)c(-c2ccnc(N)n2)s1",
    "Encorafenib":  "COC(=O)Nc1ccc(Oc2cc(NC(=O)c3ccc(F)cn3)ncn2)c(Cl)c1",
}

COMMON_TOKENS = [
    "[C]", "[=C]", "[N]", "[=N]", "[O]", "[=O]", "[S]", "[F]", "[Cl]",
    "[Br]", "[Ring1]", "[Branch1]", "[Branch2]", "[#C]", "[NH]", "[OH]",
    "[C@@H1]", "[C@H1]", "[Ring2]", "[P]",
]


# ============================================================
# GNINA 결과 파싱
# ============================================================
def parse_gnina_results(sdf_path):
    """GNINA SDF에서 화합물별 최적 pose 추출"""
    suppl = Chem.SDMolSupplier(sdf_path)
    results = {}
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
        vina = float(mol.GetProp("minimizedAffinity")) if mol.HasProp("minimizedAffinity") else 0
        cnn_score = float(mol.GetProp("CNNscore")) if mol.HasProp("CNNscore") else 0
        cnn_aff = float(mol.GetProp("CNNaffinity")) if mol.HasProp("CNNaffinity") else 0

        # 화합물별 최적 pose (CNN affinity 기준)
        if name not in results or cnn_aff > results[name]["cnn_affinity"]:
            results[name] = {
                "molecule_chembl_id": name,
                "vina_score": vina,
                "cnn_pose_score": cnn_score,
                "cnn_affinity": cnn_aff,
            }
    return pd.DataFrame(results.values())


# ============================================================
# 분자 생성 (SELFIES)
# ============================================================
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
            op = random.choice(["replace", "insert", "delete"])
            if op == "replace" and mutated:
                mutated[random.randint(0, len(mutated)-1)] = random.choice(COMMON_TOKENS)
            elif op == "insert":
                mutated.insert(random.randint(0, len(mutated)), random.choice(COMMON_TOKENS))
            elif op == "delete" and len(mutated) > 2:
                mutated.pop(random.randint(0, len(mutated)-1))
        try:
            new_smi = sf.decoder("".join(mutated))
            if new_smi is None:
                continue
            mol = Chem.MolFromSmiles(new_smi)
            if mol is None:
                continue
            canon = Chem.MolToSmiles(mol)
            mw = Descriptors.MolWt(mol)
            if 200 <= mw <= 600 and canon != smiles:
                generated.add(canon)
        except Exception:
            continue
    return generated


def scaffold_hop(smiles):
    """할로겐/작용기 치환"""
    generated = set()
    for old, replacements in [("F", ["Cl", "Br"]), ("Cl", ["F", "Br"]),
                               ("(C)", ""), ("(OC)", "(O)")]:
        if isinstance(replacements, list):
            for new in replacements:
                new_smi = smiles.replace(old, new, 1)
                m = Chem.MolFromSmiles(new_smi)
                if m:
                    mw = Descriptors.MolWt(m)
                    if 200 <= mw <= 600:
                        generated.add(Chem.MolToSmiles(m))
        elif old in smiles:
            new_smi = smiles.replace(old, replacements, 1)
            m = Chem.MolFromSmiles(new_smi)
            if m:
                mw = Descriptors.MolWt(m)
                if 200 <= mw <= 600:
                    generated.add(Chem.MolToSmiles(m))
    return generated


def evaluate_generated(gen_smiles, ref_smiles):
    """생성 분자 Drug-likeness 평가"""
    ref_fps = []
    for smi in ref_smiles:
        m = Chem.MolFromSmiles(smi)
        if m:
            ref_fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))

    results = []
    for smi in gen_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        if violations > 1:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        max_sim = max((DataStructs.TanimotoSimilarity(fp, r) for r in ref_fps), default=0)
        results.append({
            "smiles": smi, "MW": mw, "logP": logp, "HBD": hbd, "HBA": hba,
            "TPSA": Descriptors.TPSA(mol), "max_sim_to_candidates": max_sim,
        })
    return pd.DataFrame(results)


# ============================================================
# 시각화
# ============================================================
def plot_funnel(counts, fig_path):
    """파이프라인 funnel chart"""
    stages = list(counts.keys())
    values = list(counts.values())
    fig, ax = plt.subplots(figsize=(10, 7))
    colors = plt.cm.RdYlGn(np.linspace(0.2, 0.9, len(stages)))
    y_pos = np.linspace(0.9, 0.1, len(stages))

    for i, (stage, val) in enumerate(zip(stages, values)):
        width = max(0.08, 0.9 * val / values[0])
        rect = plt.Rectangle((0.5-width/2, y_pos[i]-0.03), width, 0.06,
                              facecolor=colors[i], edgecolor="white", linewidth=2, zorder=2)
        ax.add_patch(rect)
        ax.text(0.5, y_pos[i], f"{stage}\n{val:,}",
                ha="center", va="center", fontsize=10, fontweight="bold", zorder=3)
        if i < len(stages)-1:
            ax.annotate("", xy=(0.5, y_pos[i+1]+0.04), xytext=(0.5, y_pos[i]-0.04),
                        arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title("Virtual Screening Pipeline", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_docking_results(df, fig_path):
    """GNINA 도킹 결과 시각화"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("GNINA Docking Results", fontsize=14, fontweight="bold")

    axes[0].hist(df["cnn_affinity"], bins=30, color="#2ca02c", edgecolor="white", alpha=0.8)
    axes[0].set_xlabel("CNN Affinity (pKd)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("CNN Affinity Distribution")

    axes[1].hist(df["vina_score"], bins=30, color="#1f77b4", edgecolor="white", alpha=0.8)
    axes[1].set_xlabel("Vina Score (kcal/mol)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Vina Score Distribution")

    axes[2].scatter(df["vina_score"], df["cnn_affinity"], alpha=0.5, s=15, c="#ff7f0e")
    axes[2].set_xlabel("Vina Score (kcal/mol)")
    axes[2].set_ylabel("CNN Affinity (pKd)")
    axes[2].set_title("Vina vs CNN Affinity")

    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_final_scatter(df, df_final, fig_path):
    """최종 후보 다차원 scatter"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Final Candidates — Multi-dimensional Analysis", fontsize=14, fontweight="bold")

    rest = df[~df["molecule_chembl_id"].isin(df_final["molecule_chembl_id"])]

    # MW vs logP
    if "MW" in rest.columns:
        axes[0].scatter(rest["MW"], rest["logP"], alpha=0.3, s=10, c="gray", label="Other")
        axes[0].scatter(df_final["MW"], df_final["logP"], s=80, c="red", marker="*", zorder=5, label="Final 10")
        axes[0].axvline(500, color="blue", linestyle="--", alpha=0.3)
        axes[0].axhline(5, color="blue", linestyle="--", alpha=0.3)
        axes[0].set_xlabel("Molecular Weight"); axes[0].set_ylabel("logP")
        axes[0].set_title("Drug-likeness"); axes[0].legend(fontsize=8)

    # ChemBERTa vs GNINA
    if "predicted_active_prob" in rest.columns and "cnn_affinity" in rest.columns:
        axes[1].scatter(rest["predicted_active_prob"], rest["cnn_affinity"], alpha=0.3, s=10, c="gray")
        axes[1].scatter(df_final["predicted_active_prob"], df_final["cnn_affinity"],
                        s=80, c="red", marker="*", zorder=5)
        axes[1].set_xlabel("ChemBERTa Active Prob."); axes[1].set_ylabel("GNINA CNN Affinity")
        axes[1].set_title("ML Prediction vs Docking")

    # ADMET vs GNINA
    if "admet_score" in rest.columns and "cnn_affinity" in rest.columns:
        axes[2].scatter(rest["admet_score"], rest["cnn_affinity"], alpha=0.3, s=10, c="gray")
        axes[2].scatter(df_final["admet_score"], df_final["cnn_affinity"],
                        s=80, c="red", marker="*", zorder=5)
        axes[2].set_xlabel("ADMET Score"); axes[2].set_ylabel("GNINA CNN Affinity")
        axes[2].set_title("Safety vs Binding")

    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_generated(df_gen, seed_mols, fig_path):
    """생성 분자 시각화"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Generated Molecules — Chemical Space Exploration", fontsize=14, fontweight="bold")

    if len(df_gen) > 0:
        axes[0].scatter(df_gen["MW"], df_gen["logP"], alpha=0.3, s=15, c="#2ca02c", label=f"Generated ({len(df_gen)})")
    for smi in seed_mols:
        m = Chem.MolFromSmiles(smi)
        if m:
            axes[0].scatter(Descriptors.MolWt(m), Descriptors.MolLogP(m), s=100, c="red", marker="*", zorder=5)
    axes[0].scatter([], [], s=100, c="red", marker="*", label="Seed")
    axes[0].set_xlabel("MW"); axes[0].set_ylabel("logP"); axes[0].set_title("Chemical Space"); axes[0].legend(fontsize=8)

    if len(df_gen) > 0:
        axes[1].hist(df_gen["max_sim_to_candidates"], bins=30, color="#ff7f0e", edgecolor="white", alpha=0.8)
    axes[1].set_xlabel("Tanimoto Sim. to Seeds"); axes[1].set_ylabel("Count"); axes[1].set_title("Novelty Distribution")

    if len(df_gen) > 0:
        cats = pd.cut(df_gen["max_sim_to_candidates"], bins=[0,0.4,0.6,0.8,1.0],
                       labels=["Novel","Moderate","Similar","Very Similar"])
        counts = cats.value_counts()
        if len(counts) > 0:
            axes[2].pie(counts.values, labels=counts.index, autopct="%1.0f%%",
                        colors=["#2ca02c","#ff7f0e","#1f77b4","#d62728"])
            axes[2].set_title("Novelty Categories")

    plt.tight_layout()
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()


# ============================================================
# Main
# ============================================================
def main():
    # === 1. GNINA 결과 파싱 ===
    print("Parsing GNINA docking results...")
    df_dock = parse_gnina_results(DOCKED_SDF)
    print(f"  Docked compounds: {len(df_dock)}")

    plot_docking_results(df_dock, os.path.join(FIGURES_DIR, "fig07_docking_results.png"))
    print("  fig07 saved.")

    # === 2. ADMET 데이터와 병합 ===
    print("Merging ADMET + docking data...")
    df_admet = pd.read_csv(ADMET_FILE)
    df = df_admet.merge(df_dock, on="molecule_chembl_id", how="inner")
    print(f"  Merged compounds: {len(df)}")

    # === 3. 종합 점수 ===
    # 정규화
    df["norm_chemberta"] = df["predicted_active_prob"]
    df["norm_admet"] = df["admet_score"]
    df["norm_cnn"] = (df["cnn_affinity"] - df["cnn_affinity"].min()) / \
                     (df["cnn_affinity"].max() - df["cnn_affinity"].min())
    df["norm_vina"] = (-df["vina_score"] - (-df["vina_score"]).min()) / \
                      ((-df["vina_score"]).max() - (-df["vina_score"]).min())

    # 가중 종합 점수
    df["final_score"] = (
        0.25 * df["norm_chemberta"] +   # ChemBERTa 활성 예측
        0.25 * df["norm_admet"] +        # ADMET 약동학
        0.30 * df["norm_cnn"] +          # GNINA CNN affinity (가장 중요)
        0.20 * df["norm_vina"]           # GNINA Vina score
    )
    df = df.sort_values("final_score", ascending=False).reset_index(drop=True)

    # === 4. 화학적 다양성 보장 → 최종 10개 ===
    print("Selecting final diverse candidates...")
    selected_idx, selected_fps = [], []
    for i, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        too_similar = any(DataStructs.TanimotoSimilarity(fp, s) > 0.7 for s in selected_fps)
        if not too_similar:
            selected_idx.append(i)
            selected_fps.append(fp)
        if len(selected_idx) >= 10:
            break

    df_final = df.iloc[selected_idx].reset_index(drop=True)
    print(f"  Final candidates: {len(df_final)}")

    # 분자 구조 이미지
    mols_img, legends_img = [], []
    for _, row in df_final.iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol:
            AllChem.Compute2DCoords(mol)
            mols_img.append(mol)
            legends_img.append(f"{row['molecule_chembl_id']}\nCNN={row['cnn_affinity']:.1f}")
    for name, smi in REFERENCE_DRUGS.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            AllChem.Compute2DCoords(mol)
            mols_img.append(mol)
            legends_img.append(f"{name}\n(Reference)")
    if mols_img:
        img = Draw.MolsToGridImage(mols_img, molsPerRow=5, subImgSize=(400, 350), legends=legends_img)
        img.save(os.path.join(FIGURES_DIR, "fig08_candidate_structures.png"))
        print("  fig08 saved.")

    # scatter 비교
    plot_final_scatter(df, df_final, os.path.join(FIGURES_DIR, "fig09_scatter_analysis.png"))
    print("  fig09 saved.")

    # === 5. Pipeline funnel ===
    funnel = {
        "ChEMBL Drug-like Library": 50000,
        "Lipinski + PAINS Filter": 47050,
        "ChemBERTa Screening (Top 1K)": 1000,
        "ADMET Filter (Top 300)": 300,
        "GNINA Docking": len(df),
        "Final Diverse Candidates": len(df_final),
    }
    plot_funnel(funnel, os.path.join(FIGURES_DIR, "fig10_pipeline_funnel.png"))
    print("  fig10 saved.")

    # === 6. 분자 생성 (SELFIES) ===
    print("\nGenerating novel molecules (SELFIES)...")
    seed_smiles = df_final["canonical_smiles"].tolist()
    all_gen = set()

    for smi in seed_smiles:
        all_gen.update(selfies_mutate(smi, n_mutations=3, n_attempts=100))
        all_gen.update(scaffold_hop(smi))
    all_gen -= set(seed_smiles)
    print(f"  Raw generated: {len(all_gen)}")

    df_gen = evaluate_generated(list(all_gen), seed_smiles)
    print(f"  Drug-like generated: {len(df_gen)}")

    plot_generated(df_gen, seed_smiles, os.path.join(FIGURES_DIR, "fig11_generated_molecules.png"))
    print("  fig11 saved.")

    # Novel 분자 구조 이미지
    if len(df_gen) > 0:
        df_gen_sorted = df_gen.sort_values("max_sim_to_candidates")
        top_novel = df_gen_sorted.head(10)
        nmols, nlegs = [], []
        for _, row in top_novel.iterrows():
            m = Chem.MolFromSmiles(row["smiles"])
            if m:
                AllChem.Compute2DCoords(m)
                nmols.append(m)
                nlegs.append(f"Sim={row['max_sim_to_candidates']:.2f}\nMW={row['MW']:.0f}")
        if nmols:
            img = Draw.MolsToGridImage(nmols, molsPerRow=5, subImgSize=(350, 300), legends=nlegs)
            img.save(os.path.join(FIGURES_DIR, "fig12_novel_generated.png"))
            print("  fig12 saved.")

    # === 저장 ===
    df.to_csv(os.path.join(OUTPUT_DIR, "all_ranked.csv"), index=False)
    df_final.to_csv(os.path.join(OUTPUT_DIR, "final_candidates.csv"), index=False)
    df_gen.to_csv(os.path.join(OUTPUT_DIR, "generated_molecules.csv"), index=False)

    # === 요약 ===
    print(f"\n{'='*60}")
    print(f"FINAL CANDIDATES")
    print(f"{'='*60}")
    show_cols = ["molecule_chembl_id", "predicted_active_prob", "cnn_affinity",
                 "vina_score", "admet_score", "final_score"]
    show_cols = [c for c in show_cols if c in df_final.columns]
    print(df_final[show_cols].to_string(index=False))
    print(f"\nGenerated molecules: {len(df_gen)} drug-like")
    if len(df_gen) > 0:
        novel = (df_gen["max_sim_to_candidates"] < 0.4).sum()
        print(f"  Novel (Tanimoto < 0.4): {novel} ({novel/len(df_gen)*100:.0f}%)")
    print(f"\nAll outputs in: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
