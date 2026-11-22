# 모델 선정 근거

## 계산 환경
- CPU: AMD Ryzen 5 9600X 6-Core (GPU 미보유)
- Docker Desktop (WSL2 Ubuntu) — GNINA 실행용
- 모든 ML/DL 실험은 CPU에서 수행, 도킹만 Docker 컨테이너 사용

## 모델 1: RDKit (물리화학적 필터링)
- Lipinski Rule of Five: MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 (1회 위반 허용)
- PAINS 필터: 실험 가짜 양성 화합물 제거
- 결과: 50,000 → 47,050 (5.9% 제거)

## 모델 2: ChemBERTa (Transformer 기반 활성 예측)
- 모델: seyonec/ChemBERTa-zinc-base-v1 (HuggingFace)
- 사전학습: ZINC 7700만 분자, RoBERTa 아키텍처
- 임베딩: CLS 토큰 768차원 → PCA(50차원) → RandomForest 분류
- 학습 데이터: EGFR 활성 5,002개 (Active pChEMBL ≥ 7.0 / Inactive < 6.0)
- 성능: ROC-AUC 0.897 ± 0.014 (5-fold CV)
- 데이터 누수 방지: cross_val_predict로 out-of-fold 예측 사용
- 선정 근거: 대규모 사전학습 Transformer, CPU 실행 가능, SMILES 직접 입력

## 모델 3: GNINA (CNN 기반 분자 도킹)
- Docker 이미지: gnina/gnina
- 수용체: PDB 1M17 (EGFR 키나제 도메인, wild-type)
- 참조 리간드: Erlotinib (AQ4) — 결합 부위 정의용
- 출력: Vina score (kcal/mol) + CNN pose score + CNN affinity (pKd)
- 설정: exhaustiveness=8, num_modes=3, cpu=6
- 선정 근거: 구조 기반 도킹의 gold standard, CNN 스코어링으로 전통 도킹보다 정확

## 모델 4: ADMET-AI (Chemprop 앙상블)
- 패키지: admet_ai (Chemprop = Directed MPNN)
- 예측 속성 16종:
  - 흡수(A): HIA, Caco-2, Bioavailability, P-gp
  - 분포(D): BBB, VDss
  - 대사(M): CYP2D6/3A4/2C9 억제
  - 배설(E): Half-life, Clearance
  - 독성(T): hERG, AMES, DILI, LD50, Skin Reaction
- 선정 근거: 약물 개발 실패의 ~60%가 ADMET 문제

## 모델 5: Morgan Fingerprint (Tanimoto 유사도)
- ECFP4 (radius=2, 2048-bit)
- 비교 대상: Erlotinib, Gefitinib, Osimertinib, Afatinib
- 용도: 기존 약물 대비 구조적 유사성 평가 + 화학적 다양성 보장

## 모델 6: SELFIES 기반 분자 생성
- SELFIES: 100% 유효 분자 생성 보장
- 전략: 토큰 돌연변이 (replace/insert/delete) + 스캐폴드 호핑
- 결과: 10개 시드 → 778개 drug-like 분자 (73% novel)
- 선정 근거: REINVENT(RNN)의 경량 대안, GPU 불필요

## 종합 점수 가중치

| 요소 | 가중치 | 근거 |
|------|--------|------|
| GNINA CNN affinity | 0.30 | 구조 기반 결합 예측 (가장 신뢰도 높음) |
| ChemBERTa 활성 확률 | 0.25 | 리간드 기반 ML 예측 |
| ADMET 점수 | 0.25 | 약물 적합성 종합 |
| GNINA Vina score | 0.20 | 물리 기반 도킹 점수 |
