"""
Step 0-A: ChEMBL에서 BRAF (CHEMBL5145) 바이오활성 데이터 수집 (모델 학습용)

타겟: BRAF (B-Raf proto-oncogene, serine/threonine kinase)
질병: 흑색종 (melanoma)
"""
import requests
import pandas as pd
import time
import os

TARGET_CHEMBL_ID = "CHEMBL5145"  # BRAF
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "braf_training_data.csv")
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


def fetch_all_activities():
    all_activities = []
    offset = 0
    limit = 1000
    params = {
        "target_chembl_id": TARGET_CHEMBL_ID,
        "standard_type": "IC50",
        "standard_relation": "=",
        "limit": limit,
    }
    print(f"Fetching BRAF ({TARGET_CHEMBL_ID}) IC50 data from ChEMBL...")
    while True:
        params["offset"] = offset
        try:
            resp = requests.get(BASE_URL, params=params, timeout=60)
            resp.raise_for_status()
        except requests.RequestException as e:
            print(f"  Error at offset {offset}: {e}")
            time.sleep(5)
            continue
        data = resp.json()
        activities = data.get("activities", [])
        if not activities:
            break
        all_activities.extend(activities)
        total = data["page_meta"]["total_count"]
        print(f"  {len(all_activities)}/{total}")
        if data["page_meta"]["next"] is None:
            break
        offset += limit
        time.sleep(0.3)
    print(f"Total: {len(all_activities)}")
    return all_activities


def process_activities(activities):
    records = []
    for act in activities:
        smiles = act.get("canonical_smiles")
        pchembl = act.get("pchembl_value")
        if smiles and pchembl is not None:
            records.append({
                "molecule_chembl_id": act.get("molecule_chembl_id"),
                "canonical_smiles": smiles,
                "pchembl_value": float(pchembl),
                "standard_value_nM": float(act["standard_value"]) if act.get("standard_value") else None,
            })
    df = pd.DataFrame(records)
    df = df[df["pchembl_value"] >= 5.0].copy()
    df = df.sort_values("pchembl_value", ascending=False)
    df = df.drop_duplicates(subset=["canonical_smiles"], keep="first")
    df["activity_label"] = -1
    df.loc[df["pchembl_value"] >= 7.0, "activity_label"] = 1
    df.loc[df["pchembl_value"] < 6.0, "activity_label"] = 0
    return df.reset_index(drop=True)


def main():
    activities = fetch_all_activities()
    df = process_activities(activities)
    df.to_csv(OUTPUT_FILE, index=False)
    n_a = (df["activity_label"] == 1).sum()
    n_i = (df["activity_label"] == 0).sum()
    print(f"\n=== BRAF Training Data ===")
    print(f"Total: {len(df)}, Active: {n_a}, Inactive: {n_i}, Ambiguous: {len(df)-n_a-n_i}")
    print(f"Output: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
