# EGFR 타겟 비소세포폐암 신약 후보 물질 Virtual Screening 보고서 작성 프롬프트

## 역할
당신은 약학/신약개발 분야의 연구 보고서를 작성하는 AI 조교입니다. 아래 지시에 따라 experiment 폴더의 실험 데이터와 그림을 활용하여 보고서를 작성해 주세요.

## 보고서 형식
- **최대 20장** (짧을수록 가산점)
- **본문 글자 크기 11pt**
- 글은 최소화하고, **그림과 표 중심으로 구성**
- 각 그림/표 아래에 간단한 캡션(2-3줄) 추가
- 한국어로 작성

## 문체 규칙
- 학술 보고서 문체 사용 (예: "~하였다", "~이다")
- 불필요한 수식어 배제, 간결하게 서술
- AI 특유의 과잉 설명 금지 — 한 문단은 최대 4줄
- 목록/표/그래프로 대체 가능한 내용은 반드시 시각화

## 보고서 구조

### 1. 서론 (1장)
- 비소세포폐암(NSCLC) 개요: 폐암의 85%를 차지, 5년 생존율 약 25%
- EGFR 타겟 선정 근거: NSCLC 환자의 약 40%에서 EGFR 돌연변이 확인
- 기존 EGFR 억제제: Erlotinib(1세대), Gefitinib(1세대), Afatinib(2세대), Osimertinib(3세대)
- 연구 목적: 딥러닝 기반 virtual screening으로 새로운 EGFR 억제제 후보 발굴

### 2. 타겟 단백질 정보 (0.5장)
- 단백질: EGFR (Epidermal Growth Factor Receptor)
- ChEMBL ID: CHEMBL203
- PDB ID: 1M17 (EGFR 키나제 도메인 + Erlotinib 결합 구조, 2.60 Å)
- UniProt: P00533
- 표로 정리 (ID, 해상도, 리간드, 출처 등)

### 3. Virtual Screening 파이프라인 (1장)
- experiment/figures/fig08_pipeline_funnel.png 삽입
- 각 단계 한 줄 설명:
  1. ChEMBL 데이터 수집: EGFR 타겟 IC50 활성 데이터 19,719건 → 10,005개 고유 화합물
  2. Lipinski + PAINS 필터: 약물유사성 필터 → 6,608개
  3. DeepDTA 결합 친화도 예측: MPNN_CNN 모델 (Davis 데이터셋 사전학습)
  4. ChemBERTa 활성 분류: Transformer 기반 분자 임베딩 + RandomForest
  5. ADMET 예측: Chemprop 기반 흡수/분포/대사/배설/독성 평가
  6. 최종 후보 선정: 종합 점수 기반 10개 (화학적 다양성 보장)
  7. 분자 생성: SELFIES 돌연변이 + 스캐폴드 호핑 → 718개 신규 분자

### 4. 모델 선정 및 설명 (1.5장)
아래 표로 정리:

| 모델 | 유형 | 입력 | 출력 | 선정 근거 |
|------|------|------|------|-----------|
| RDKit | 규칙 기반 | SMILES | MW, logP, HBD, HBA 등 | Lipinski Rule of Five 필터링 |
| Morgan Fingerprint | 분자 지문 | SMILES | 2048-bit 벡터 | Tanimoto 유사도로 기존 약물과 비교 |
| DeepPurpose (MPNN_CNN) | GNN + CNN | SMILES + 단백질 시퀀스 | Kd 예측값 | 사전학습 모델, GPU 불필요 |
| ChemBERTa | Transformer | SMILES 토큰 | 768차원 임베딩 | 대규모 사전학습, 활성 분류에 활용 |
| ADMET-AI (Chemprop) | GNN | SMILES | ADMET 속성 16종 | 약동학 평가, Chemprop 앙상블 |
| SELFIES Generator | 문법 기반 생성 | Seed SMILES | 돌연변이 SMILES | 100% 유효 분자 생성 보장 |

- 계산 자원 제약 설명: GPU 미보유, CPU(AMD Ryzen 5 9600X)에서 모든 실험 수행
- 구조 기반 도킹(GNINA) 대신 시퀀스 기반(DeepDTA) 선택 근거

### 5. 실험 결과 (3-4장, 그림 중심)
#### 5.1 물리화학적 필터링
- experiment/figures/fig01_lipinski_filtering.png
- 10,005 → 6,608개 (34% 제거)

#### 5.2 분자 유사도 분석
- experiment/figures/fig02_similarity_distribution.png
- 기존 EGFR 약물 4종 대비 Tanimoto 유사도

#### 5.3 결합 친화도 예측 (DeepDTA)
- experiment/figures/fig03_binding_affinity.png
- 실험값(pChEMBL)과 예측값(pKd) 상관관계: r = 0.032
- 상관관계가 낮은 이유: IC50와 Kd는 다른 측정값, 사전학습 데이터셋(Davis)의 도메인 차이

#### 5.4 ChemBERTa 활성 분류
- experiment/figures/fig04_chemberta_analysis.png
- 분류 ROC-AUC: 0.831 ± 0.008 (5-fold CV)
- t-SNE에서 활성/비활성 분자의 클러스터 분리 확인

#### 5.5 ADMET 프로파일
- experiment/figures/fig05_admet_heatmap.png
- experiment/figures/fig06_admet_radar.png
- experiment/figures/fig07_admet_distributions.png
- 주요 ADMET 기준: HIA(흡수), hERG(심장독성), AMES(변이원성), DILI(간독성)

#### 5.6 분자 생성 & 화학 공간 탐색
- experiment/figures/fig12_generated_molecules.png
- experiment/figures/fig13_novel_generated.png
- 10개 시드 → 718개 drug-like 분자, 66% novel (Tanimoto < 0.4)

### 6. 최종 후보 물질 (1.5장)
- experiment/figures/fig09_candidate_structures.png (분자 구조)
- experiment/figures/fig10_comparison_table.png (기존 약물 비교)
- experiment/figures/fig11_scatter_analysis.png (다차원 비교)
- 최종 10개 후보의 ChEMBL ID, 구조, pChEMBL, MW, logP, ADMET 종합 점수를 표로 정리
- 각 후보가 기존 약물과 비교하여 가진 장점 한 줄 설명

### 7. 결론 (0.5장)
- 6개 딥러닝/ML 모델을 활용한 체계적 virtual screening 수행
- CPU 환경에서도 실현 가능한 파이프라인 구축
- SELFIES 기반 분자 생성으로 718개 신규 화합물 탐색
- 한계: 실험적 검증 부재, 도킹 시뮬레이션 미수행, DeepDTA 상관관계 한계
- 향후 과제: GNINA 도킹 추가, 분자 동역학 시뮬레이션, 실험적 IC50 측정

## 주의사항
- 그림 파일은 experiment/figures/ 폴더에 위치
- 데이터는 experiment/ 하위 폴더의 CSV 파일 참조
- 20장을 초과하지 말 것
- 불필요한 반복 설명 금지
