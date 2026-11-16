# 모델 선정 근거 및 기술 상세

## 계산 환경 제약
- CPU: AMD Ryzen 5 9600X 6-Core (GPU 미보유)
- RAM: 로컬 PC 환경
- 모든 실험은 CPU에서 수행
- 구조 기반 도킹(GNINA)은 Linux + GPU 필요 → 시퀀스 기반 대안(DeepDTA) 선택

## 모델 1: RDKit (물리화학적 필터링)
- **용도**: Lipinski Rule of Five + PAINS 필터
- **기준**: MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 (1개 위반 허용)
- **PAINS**: pan-assay interference compounds 제거
- **결과**: 10,005 → 6,608 화합물 (34% 제거)

## 모델 2: Morgan Fingerprint (분자 유사도)
- **알고리즘**: Extended-Connectivity Fingerprint (ECFP4), radius=2, 2048-bit
- **비교 대상**: Erlotinib, Gefitinib, Osimertinib, Afatinib
- **메트릭**: Tanimoto coefficient
- **용도**: 기존 약물 대비 구조적 유사성 평가

## 모델 3: DeepPurpose — MPNN_CNN (결합 친화도 예측)
- **아키텍처**: Message Passing Neural Network (drug) + CNN (target protein)
- **사전학습**: Davis kinase dataset
- **입력**: SMILES (drug) + 아미노산 시퀀스 (EGFR kinase domain)
- **출력**: 예측 Kd (nM)
- **선정 근거**: 단백질 3D 구조 없이 시퀀스만으로 결합 친화도 예측 가능, CPU 실행 가능
- **한계**: Davis 데이터셋은 키나제 편향, IC50/Kd 메트릭 차이로 실험값과의 직접 비교에 한계

## 모델 4: ChemBERTa (Transformer 기반 분자 표현)
- **모델**: seyonec/ChemBERTa-zinc-base-v1 (HuggingFace)
- **사전학습**: ZINC 데이터베이스
- **아키텍처**: RoBERTa (Transformer encoder)
- **임베딩**: CLS 토큰, 768차원 → PCA(50차원) → RandomForest
- **분류 성능**: ROC-AUC 0.831 ± 0.008 (5-fold CV)
- **활성 기준**: pChEMBL ≥ 7.0 (IC50 ≤ 100 nM)
- **선정 근거**: 대규모 사전학습으로 범용 분자 표현 학습, 다양한 downstream task 적용 가능

## 모델 5: ADMET-AI (Chemprop 앙상블)
- **패키지**: admet_ai v2.0.1
- **아키텍처**: Chemprop (Directed MPNN) 앙상블
- **예측 속성 (16종)**:
  - 흡수(A): HIA, Caco-2, Bioavailability, P-gp
  - 분포(D): BBB, VDss
  - 대사(M): CYP2D6, CYP3A4, CYP2C9
  - 배설(E): Half-life, Clearance
  - 독성(T): hERG, AMES, DILI, LD50, Skin Reaction
- **필터링**: ADMET 종합 점수(가중 평균) 상위 50% 선택
- **선정 근거**: 약물 개발 초기 단계에서 fail 원인의 ~60%가 ADMET 문제

## 모델 6: SELFIES 기반 분자 생성
- **SELFIES**: Self-Referencing Embedded Strings — 100% 유효 분자 생성 보장
- **전략**:
  1. SELFIES 토큰 돌연변이 (replace/insert/delete, 1-3회)
  2. 스캐폴드 호핑 (할로겐 교체, 작용기 변환)
- **결과**: 10개 시드 → 718개 drug-like 분자 생성 (66% novel)
- **선정 근거**: REINVENT(RNN 기반)의 경량 대안으로, GPU 없이 화학 공간 탐색 가능

## 종합 점수 가중치
| 요소 | 가중치 | 설명 |
|------|--------|------|
| pChEMBL (실험 활성) | 0.30 | IC50 기반 실험적 활성 |
| ChemBERTa 활성 확률 | 0.20 | Transformer 기반 활성 예측 |
| DeepDTA pKd | 0.15 | 시퀀스 기반 결합 친화도 |
| ADMET 점수 | 0.25 | 약동학 종합 평가 |
| 구조 신규성 | 0.10 | 기존 약물 대비 적절한 차별화 |
