"""
Step 2-B: 용해도 회귀 (Solubility Regression)

ESOL(Delaney) 데이터로 MLP를 학습하고, 스크리닝 후보의 logS를 예측한다.
- 입력: SMILES → 6개 분자 기술자 (MW, logP, HBD, HBA, TPSA, RotBonds)
- 출력: 예측 logS (log10(mol/L))
- 판정: logS >= -4 → 용해도 기준 통과

참고: 교수 가이드 — "logS 값이 클수록 물에 잘 녹는다"
"""
import os
import csv
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupShuffleSplit
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
EXPERIMENT = os.path.join(BASE, "..")
FIGURES_DIR = os.path.join(EXPERIMENT, "figures")
DELANEY_URL = "https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv"


def feat(mol):
    """6개 분자 기술자 추출"""
    return [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
        Lipinski.NumRotatableBonds(mol),
    ]


def main():
    # === 1. ESOL 학습 데이터 ===
    delaney_path = os.path.join(BASE, "delaney-processed.csv")
    if not os.path.exists(delaney_path):
        print("Downloading ESOL dataset...")
        import requests
        resp = requests.get(DELANEY_URL, timeout=30)
        with open(delaney_path, "w", encoding="utf-8") as f:
            f.write(resp.text)

    train = list(csv.DictReader(open(delaney_path, encoding="utf-8")))
    mols = [Chem.MolFromSmiles(r["smiles"]) for r in train]
    X = np.array([feat(m) for m in mols])
    y = np.array([float(r["measured log solubility in mols per litre"]) for r in train])
    groups = [MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m.GetNumAtoms() > 0 else "" for m in mols]

    # scaffold 기반 train/test 분리
    tr, te = next(GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=13).split(X, y, groups))
    model = make_pipeline(StandardScaler(), MLPRegressor(
        hidden_layer_sizes=(32, 16), early_stopping=True, max_iter=1000, random_state=13
    ))
    model.fit(X[tr], y[tr])
    mae = mean_absolute_error(y[te], model.predict(X[te]))
    print(f"Solubility MLP: train {len(tr)}, test {len(te)}, MAE={mae:.2f} log units")

    # === 2. ADMET 통과 후보에 대해 예측 ===
    admet_file = os.path.join(EXPERIMENT, "04_admet", "admet_top_hits.csv")
    if not os.path.exists(admet_file):
        print("ERROR: admet_top_hits.csv not found")
        return

    df = pd.read_csv(admet_file)
    print(f"\nPredicting solubility for {len(df)} ADMET-passed compounds...")

    preds = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol:
            pred = model.predict([feat(mol)])[0]
        else:
            pred = -999
        preds.append(pred)

    df["pred_logS"] = preds
    df["solubility_pass"] = df["pred_logS"] >= -4

    n_pass = df["solubility_pass"].sum()
    n_fail = len(df) - n_pass
    print(f"logS >= -4 (pass): {n_pass}")
    print(f"logS <  -4 (fail): {n_fail}")

    # === 3. 시각화 ===
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Solubility Prediction (MLP Regressor)", fontsize=14, fontweight="bold")

    # 3a. 학습 데이터 actual vs predicted
    y_pred_test = model.predict(X[te])
    axes[0].scatter(y[te], y_pred_test, alpha=0.5, s=20, c="#1f77b4")
    axes[0].plot([-8, 2], [-8, 2], "r--", alpha=0.6)
    axes[0].set_xlabel("Measured logS")
    axes[0].set_ylabel("Predicted logS")
    axes[0].set_title(f"Test Set (n={len(te)}, MAE={mae:.2f})")

    # 3b. 후보 logS 분포
    axes[1].hist(df["pred_logS"], bins=40, color="#2ca02c", edgecolor="white", alpha=0.8)
    axes[1].axvline(-4, color="red", linestyle="--", linewidth=2, label="Threshold (-4)")
    axes[1].set_xlabel("Predicted logS")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Candidate Solubility (pass={n_pass}, fail={n_fail})")
    axes[1].legend()

    # 3c. MW vs logS (pass/fail)
    passed = df[df["solubility_pass"]]
    failed = df[~df["solubility_pass"]]
    axes[2].scatter(failed["MW"], failed["pred_logS"], alpha=0.4, s=15, c="#d62728", label=f"Fail ({n_fail})")
    axes[2].scatter(passed["MW"], passed["pred_logS"], alpha=0.4, s=15, c="#2ca02c", label=f"Pass ({n_pass})")
    axes[2].axhline(-4, color="red", linestyle="--", alpha=0.6)
    axes[2].set_xlabel("Molecular Weight")
    axes[2].set_ylabel("Predicted logS")
    axes[2].set_title("MW vs Solubility")
    axes[2].legend()

    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "fig_solubility.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure: {fig_path}")

    # === 4. 저장 ===
    df.to_csv(os.path.join(EXPERIMENT, "04_admet", "solubility_scored.csv"), index=False)
    print(f"\nSaved: solubility_scored.csv")
    print(f"\n=== Solubility Summary ===")
    print(f"Model: MLP (32-16), ESOL {len(train)} compounds")
    print(f"MAE: {mae:.2f} log units")
    print(f"Candidates: {len(df)} → pass={n_pass}, fail={n_fail}")


if __name__ == "__main__":
    main()
