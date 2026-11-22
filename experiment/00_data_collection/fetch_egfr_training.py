"""
Step 0-A: ChEMBL에서 EGFR (CHEMBL203) 바이오활성 데이터 수집 (모델 학습용)

이 데이터는 모델 학습에 사용된다:
- Active: pChEMBL >= 7.0 (IC50 <= 100 nM)
- Inactive: pChEMBL < 6.0 (IC50 > 1 uM)
- 중간값(6.0~7.0)은 학습에서 제외 (ambiguous zone)

Virtual screening 대상이 아님! (이미 활성이 알려진 화합물)
"""
import requests
import pandas as pd
import time
import os

TARGET_CHEMBL_ID = "CHEMBL203"  # EGFR
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "egfr_training_data.csv")

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


def fetch_all_activities():
    """ChEMBL API에서 EGFR IC50 활성 데이터 전체 다운로드"""
    all_activities = []
    offset = 0
    limit = 1000

    params = {
        "target_chembl_id": TARGET_CHEMBL_ID,
        "standard_type": "IC50",
        "standard_relation": "=",
        "limit": limit,
    }

    print("Fetching EGFR (CHEMBL203) IC50 data from ChEMBL...")

    while True:
        params["offset"] = offset
        try:
            resp = requests.get(BASE_URL, params=params, timeout=60)
            resp.raise_for_status()
        except requests.RequestException as e:
            print(f"  Request error at offset {offset}: {e}")
            time.sleep(5)
            continue

        data = resp.json()
        activities = data.get("activities", [])
        if not activities:
            break

        all_activities.extend(activities)
        total = data["page_meta"]["total_count"]
        print(f"  {len(all_activities)}/{total} records...")

        if data["page_meta"]["next"] is None:
            break
        offset += limit
        time.sleep(0.3)

    print(f"Total fetched: {len(all_activities)}")
    return all_activities


def process_activities(activities):
    """활성 데이터 정리"""
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

    # 동일 SMILES 중복 제거 (가장 높은 pChEMBL 유지)
    df = df.sort_values("pchembl_value", ascending=False)
    df = df.drop_duplicates(subset=["canonical_smiles"], keep="first")

    # 활성 라벨 부여
    df["activity_label"] = -1  # ambiguous
    df.loc[df["pchembl_value"] >= 7.0, "activity_label"] = 1   # active (IC50 <= 100 nM)
    df.loc[df["pchembl_value"] < 6.0, "activity_label"] = 0    # inactive (IC50 > 1 uM)

    df = df.reset_index(drop=True)
    return df


def main():
    activities = fetch_all_activities()
    df = process_activities(activities)
    df.to_csv(OUTPUT_FILE, index=False)

    n_active = (df["activity_label"] == 1).sum()
    n_inactive = (df["activity_label"] == 0).sum()
    n_ambiguous = (df["activity_label"] == -1).sum()

    print(f"\n=== EGFR Training Data Summary ===")
    print(f"Total unique compounds: {len(df)}")
    print(f"Active (pChEMBL >= 7.0):   {n_active}")
    print(f"Inactive (pChEMBL < 6.0):  {n_inactive}")
    print(f"Ambiguous (6.0~7.0):       {n_ambiguous} (excluded from training)")
    print(f"pChEMBL range: {df['pchembl_value'].min():.2f} - {df['pchembl_value'].max():.2f}")
    print(f"Output: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
