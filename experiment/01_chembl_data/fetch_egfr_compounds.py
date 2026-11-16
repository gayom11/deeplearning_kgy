"""
Step 1: ChEMBL에서 EGFR (CHEMBL203) 타겟의 활성 화합물 데이터 다운로드
- IC50 데이터만 필터링
- pChEMBL value >= 5 (IC50 <= 10 uM) 인 화합물만 수집
"""
import requests
import pandas as pd
import json
import time
import os

TARGET_CHEMBL_ID = "CHEMBL203"  # EGFR
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "egfr_compounds_raw.csv")

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"


def fetch_all_activities():
    """ChEMBL API에서 EGFR IC50 활성 데이터 전체 다운로드"""
    all_activities = []
    offset = 0
    limit = 1000  # API 최대 한도

    params = {
        "target_chembl_id": TARGET_CHEMBL_ID,
        "standard_type": "IC50",
        "standard_relation": "=",
        "limit": limit,
        "offset": offset,
    }

    print(f"Fetching EGFR (CHEMBL203) IC50 data from ChEMBL...")

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
        print(f"  Fetched {len(all_activities)}/{total} records...")

        if data["page_meta"]["next"] is None:
            break

        offset += limit
        time.sleep(0.2)  # API 부하 방지

    print(f"Total fetched: {len(all_activities)} records")
    return all_activities


def process_activities(activities):
    """활성 데이터 정리 및 필터링"""
    records = []
    for act in activities:
        smiles = act.get("canonical_smiles")
        pchembl = act.get("pchembl_value")
        ic50 = act.get("standard_value")
        mol_id = act.get("molecule_chembl_id")
        mol_name = act.get("molecule_pref_name")

        if smiles and pchembl is not None:
            records.append({
                "molecule_chembl_id": mol_id,
                "molecule_name": mol_name or "",
                "canonical_smiles": smiles,
                "standard_type": act.get("standard_type"),
                "standard_value_nM": float(ic50) if ic50 else None,
                "standard_units": act.get("standard_units"),
                "pchembl_value": float(pchembl),
                "assay_chembl_id": act.get("assay_chembl_id"),
                "assay_type": act.get("assay_type"),
                "document_year": act.get("document_year"),
            })

    df = pd.DataFrame(records)

    # pChEMBL >= 5 (IC50 <= 10 uM) 필터
    df = df[df["pchembl_value"] >= 5.0].copy()

    # 동일 분자 중복 제거 (가장 높은 pChEMBL 값 유지)
    df = df.sort_values("pchembl_value", ascending=False)
    df = df.drop_duplicates(subset=["canonical_smiles"], keep="first")

    df = df.reset_index(drop=True)
    return df


def main():
    activities = fetch_all_activities()

    # Raw JSON 백업 저장
    raw_path = os.path.join(OUTPUT_DIR, "egfr_activities_raw.json")
    with open(raw_path, "w", encoding="utf-8") as f:
        json.dump(activities, f, indent=2)
    print(f"Raw JSON saved: {raw_path}")

    df = process_activities(activities)
    df.to_csv(OUTPUT_FILE, index=False)

    print(f"\n=== Summary ===")
    print(f"Total unique compounds (pChEMBL >= 5): {len(df)}")
    print(f"pChEMBL range: {df['pchembl_value'].min():.2f} - {df['pchembl_value'].max():.2f}")
    print(f"IC50 range: {df['standard_value_nM'].min():.1f} - {df['standard_value_nM'].max():.1f} nM")
    print(f"Output: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
