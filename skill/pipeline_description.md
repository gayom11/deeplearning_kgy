# Virtual Screening 파이프라인 상세

## 전체 흐름
```
ChEMBL (CHEMBL203, EGFR)
│ 19,719 IC50 records → 10,005 unique compounds (pChEMBL ≥ 5)
▼
[Step 1] RDKit — Lipinski + PAINS Filter
│ MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 (1회 위반 허용)
│ → 6,608 compounds
▼
[Step 2] Morgan Fingerprint — Tanimoto Similarity
│ 기존 EGFR 약물 4종 대비 유사도 계산 (ECFP4, 2048-bit)
▼
[Step 3] DeepPurpose (MPNN_CNN_DAVIS) — 결합 친화도 예측
│ EGFR kinase domain 시퀀스 + SMILES → pKd 예측
▼
[Step 4] ChemBERTa — Transformer 임베딩 + 활성 분류
│ 768-dim → PCA(50-dim) → RandomForest → ROC-AUC 0.831
▼
[Step 5] ADMET-AI (Chemprop) — 약동학 평가
│ Combined Score 상위 500개 → HIA/hERG/AMES/DILI 등 16속성 → 상위 50% 선별 (250개)
▼
[Step 6] 종합 점수 & 화학적 다양성 필터
│ 가중 평균 + Tanimoto 중복 제거 (> 0.8) → 최종 10개
▼
[Step 7] SELFIES 돌연변이 + 스캐폴드 호핑
│ 10개 시드 → 718개 drug-like 신규 분자 (473개 novel)
```

## 코드 위치
| 단계 | 파일 |
|------|------|
| ChEMBL 다운로드 | experiment/01_chembl_data/fetch_egfr_compounds.py |
| Lipinski + PAINS | experiment/02_filtering/lipinski_filter.py |
| DeepDTA | experiment/03_docking/binding_affinity_deeppurpose.py |
| ChemBERTa | experiment/03_docking/chemberta_classification.py |
| ADMET | experiment/04_admet/admet_prediction.py |
| 최종 선정 | experiment/05_final_candidates/select_candidates.py |
| 분자 생성 | experiment/05_final_candidates/molecular_generation.py |

## 데이터 위치
| 파일 | 내용 | 행 수 |
|------|------|-------|
| experiment/01_chembl_data/egfr_compounds_raw.csv | ChEMBL 원본 | 10,005 |
| experiment/02_filtering/egfr_filtered.csv | Lipinski 통과 | 6,608 |
| experiment/03_docking/egfr_binding_affinity.csv | DeepDTA 결과 | 5,810 |
| experiment/03_docking/egfr_chemberta_scores.csv | ChemBERTa 결과 | 6,608 |
| experiment/04_admet/egfr_admet_passed.csv | ADMET 통과 | 250 |
| experiment/05_final_candidates/final_candidates.csv | 최종 10개 | 10 |
| experiment/05_final_candidates/generated_molecules.csv | 생성 분자 | 718 |
