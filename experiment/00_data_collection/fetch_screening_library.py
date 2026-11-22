"""
Step 0-B: ChEMBL에서 Virtual Screening 대상 화합물 라이브러리 수집

스크리닝 라이브러리 = EGFR에 대해 테스트된 적 없는 drug-like 화합물
- ChEMBL 전체 화합물에서 drug-like 필터 적용
- EGFR 학습 데이터에 포함된 화합물은 제외
- 최종 ~50,000개 수집 (CPU 환경에서 처리 가능한 규모)

이 화합물들에 대해 학습된 모델로 EGFR 활성을 "예측"하는 것이 virtual screening.
"""
import requests
import pandas as pd
import time
import os
import random

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "screening_library.csv")
TRAINING_FILE = os.path.join(OUTPUT_DIR, "egfr_training_data.csv")

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
TARGET_SIZE = 50000  # 목표 스크리닝 라이브러리 크기


def load_egfr_smiles():
    """EGFR 학습 데이터의 SMILES 로드 (제외용)"""
    if not os.path.exists(TRAINING_FILE):
        print("WARNING: EGFR training data not found. Run fetch_egfr_training.py first.")
        return set()
    df = pd.read_csv(TRAINING_FILE)
    return set(df["canonical_smiles"].tolist())


def fetch_druglike_molecules(egfr_smiles, target_count=TARGET_SIZE):
    """ChEMBL에서 drug-like 화합물 다운로드 (EGFR 화합물 제외)"""
    all_molecules = []
    offset = 0
    limit = 1000
    max_offset = 500000  # API 안전장치

    # Drug-like 필터 조건
    params = {
        "molecule_properties__mw_freebase__gte": 200,
        "molecule_properties__mw_freebase__lte": 500,
        "molecule_properties__alogp__gte": -1,
        "molecule_properties__alogp__lte": 5,
        "molecule_properties__hbd__lte": 5,
        "molecule_properties__hba__lte": 10,
        "molecule_structures__canonical_smiles__isnull": False,
        "molecule_type": "Small molecule",
        "limit": limit,
        "format": "json",
    }

    print(f"Fetching drug-like molecules from ChEMBL (target: {target_count})...")

    while len(all_molecules) < target_count and offset < max_offset:
        params["offset"] = offset
        try:
            resp = requests.get(BASE_URL, params=params, timeout=120)
            resp.raise_for_status()
        except requests.RequestException as e:
            print(f"  Request error at offset {offset}: {e}")
            time.sleep(10)
            continue

        data = resp.json()
        molecules = data.get("molecules", [])
        if not molecules:
            break

        for mol in molecules:
            structs = mol.get("molecule_structures")
            if not structs:
                continue
            smiles = structs.get("canonical_smiles")
            if not smiles or smiles in egfr_smiles:
                continue  # EGFR 학습 데이터 제외

            props = mol.get("molecule_properties") or {}
            chembl_id = mol.get("molecule_chembl_id", "")

            all_molecules.append({
                "molecule_chembl_id": chembl_id,
                "canonical_smiles": smiles,
                "mw_freebase": props.get("mw_freebase"),
                "alogp": props.get("alogp"),
                "hbd": props.get("hbd"),
                "hba": props.get("hba"),
                "psa": props.get("psa"),
                "num_ro5_violations": props.get("num_ro5_violations"),
            })

        total_available = data.get("page_meta", {}).get("total_count", "?")
        print(f"  Collected {len(all_molecules)} / target {target_count} (API total: {total_available})")

        if data.get("page_meta", {}).get("next") is None:
            break
        offset += limit
        time.sleep(0.3)

    return all_molecules


def main():
    # EGFR 학습 데이터 SMILES 로드 (제외용)
    egfr_smiles = load_egfr_smiles()
    print(f"EGFR training compounds to exclude: {len(egfr_smiles)}")

    # Drug-like 화합물 다운로드
    molecules = fetch_druglike_molecules(egfr_smiles, TARGET_SIZE)

    df = pd.DataFrame(molecules)

    # SMILES 중복 제거
    df = df.drop_duplicates(subset=["canonical_smiles"]).reset_index(drop=True)

    # 목표 크기로 샘플링 (초과 시)
    if len(df) > TARGET_SIZE:
        df = df.sample(n=TARGET_SIZE, random_state=42).reset_index(drop=True)

    df.to_csv(OUTPUT_FILE, index=False)

    print(f"\n=== Screening Library Summary ===")
    print(f"Total compounds: {len(df)}")
    print(f"MW range: {df['mw_freebase'].astype(float).min():.1f} - {df['mw_freebase'].astype(float).max():.1f}")
    print(f"ALogP range: {df['alogp'].astype(float).min():.2f} - {df['alogp'].astype(float).max():.2f}")
    print(f"EGFR compounds excluded: {len(egfr_smiles)}")
    print(f"Output: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
