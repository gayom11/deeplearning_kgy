# Virtual Screening 파이프라인

## 전체 흐름
```
ChEMBL 2,921,148개 화합물
│ Drug-like 필터 (MW 200-500, logP -1~5) → 50,000개 샘플링
│ EGFR 학습 데이터 10,005개 제외
▼
[Step 1] RDKit — Lipinski + PAINS Filter
│ MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 (1회 위반 허용)
│ 50,000 → 47,050개
▼
[Step 2] ChemBERTa — EGFR 활성 예측 (Virtual Screening 핵심)
│ 학습: EGFR 활성 화합물 5,002개 (Active 3,704 / Inactive 1,298)
│ 모델: ChemBERTa CLS 임베딩(768d) → PCA(50d) → RandomForest
│ ROC-AUC: 0.897 ± 0.014 (5-fold CV, out-of-fold 예측)
│ 47,050개 스크리닝 → 상위 1,000개 선별 (활성 확률 ≥ 0.802)
▼
[Step 3] ADMET-AI (Chemprop) — 약동학 + 물리화학적 성질 평가
│ 1,000개 → HIA/hERG/AMES/DILI 등 16속성 → 상위 300개 선별
▼
[Step 4] GNINA — 구조 기반 도킹 (Docker)
│ 수용체: PDB 1M17 (EGFR 키나제 도메인)
│ Autobox: Erlotinib 결합 부위
│ 300개 화합물 → Vina score + CNN affinity
▼
[Step 5] 종합 점수 & 화학적 다양성 필터
│ ChemBERTa 25% + ADMET 25% + CNN 30% + Vina 20%
│ Tanimoto > 0.7 중복 제거 → 최종 10개
▼
[Step 6] SELFIES 돌연변이 + 스캐폴드 호핑
│ 10개 시드 → 778개 drug-like 신규 분자 (73% novel)
```

## 핵심: 학습 데이터 vs 스크리닝 대상 분리
- **학습 데이터**: ChEMBL EGFR 바이오활성 10,005개 (실험으로 활성 확인된 것 → 모델 학습용)
- **스크리닝 대상**: ChEMBL 일반 drug-like 50,000개 (EGFR 활성 미지 → 예측 대상)
- 최종 후보는 스크리닝 대상에서 나옴 (학습 데이터 재정렬이 아님)

## 코드 위치
| 단계 | 파일 |
|------|------|
| 학습 데이터 수집 | experiment/00_data_collection/fetch_egfr_training.py |
| 스크리닝 라이브러리 수집 | experiment/00_data_collection/fetch_screening_library.py |
| 전처리 | experiment/01_preprocessing/preprocess.py |
| ChemBERTa 학습+스크리닝 | experiment/02_model_training/chemberta_screen.py |
| GNINA 도킹 준비 | experiment/03_docking/prepare_docking.py |
| GNINA 실행 | experiment/03_docking/run_gnina.sh |
| ADMET 평가 | experiment/04_admet/admet_screen.py |
| 최종 선정+분자 생성 | experiment/05_final_candidates/select_final.py |
