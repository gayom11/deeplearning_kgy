# EGFR 타겟 비소세포폐암 신약 후보 물질 Virtual Screening 보고서 작성 프롬프트

## 역할
당신은 약학/신약개발 분야의 연구 보고서를 작성하는 AI 조교입니다.

## 보고서 형식
- **최대 15장** (짧을수록 가산점)
- **본문 글자 크기 11pt**
- 글은 최소화하고, **그림과 표 중심으로 구성**
- 각 그림/표 아래에 간단한 캡션(1-2줄)
- 한국어로 작성

## 문체 규칙
- 학술 보고서 문체 ("~하였다", "~이다")
- 불필요한 수식어 배제, 간결 서술
- 한 문단 최대 3줄
- 목록/표/그래프로 대체 가능한 내용은 반드시 시각화

## 보고서 구조

### 1. 서론 (0.5장)
- NSCLC 개요, EGFR 타겟 근거, 기존 EGFR 억제제 표, 연구 목적
- 표 1개 (기존 약물 비교)

### 2. 타겟 단백질 정보 (0.5장)
- EGFR 기능, ChEMBL/PDB ID, 1M17 선정 근거
- 표 1개

### 3. Virtual Screening 파이프라인 (1장)
- fig10_pipeline_funnel.png 삽입
- 핵심: 학습 데이터(모델 교육)와 스크리닝 라이브러리(예측 대상)의 분리 설명
- 각 단계 한 줄 설명

### 4. 모델 선정 (1장)
- 6종 모델 표 (모델명, 유형, 입력, 출력, 선정 근거)
- 계산 환경 제약 (CPU + Docker) 한 줄
- GNINA를 Docker로 해결한 점 명시

### 5. 실험 결과 (4-5장, 그림 중심)
#### 5.1 전처리
- fig01_preprocessing.png, 50,000 → 47,050

#### 5.2 ChemBERTa 모델 학습
- fig02_chemberta_training.png
- ROC-AUC 0.897, out-of-fold 예측 사용

#### 5.3 Virtual Screening 결과
- fig03_screening_results.png
- 47,050개 → 상위 1,000개 (활성 확률 ≥ 0.802)

#### 5.4 GNINA 도킹
- fig07_docking_results.png
- 300개 도킹, CNN affinity + Vina score

#### 5.5 ADMET 프로파일
- fig04, fig05, fig06 중 2개 선택 삽입
- 물리화학적 성질 + ADMET 결과 요약

#### 5.6 분자 생성
- fig11_generated_molecules.png
- 778개 drug-like 신규 분자, 73% novel

### 6. 최종 후보 물질 (1.5장)
- fig08_candidate_structures.png (분자 구조)
- fig09_scatter_analysis.png (다차원 비교)
- 최종 10개 표 (ID, CNN affinity, Vina score, ADMET score, 최종 점수)
- 기존 약물 대비 구조적 차이

### 7. 결론 (0.5장)
- 6개 모델 체계적 virtual screening
- 학습/스크리닝 분리의 의의
- CPU+Docker 환경에서 GNINA 도킹 실현
- 한계: 실험적 검증 부재
- 향후: in vitro IC50 측정, MD 시뮬레이션

## 주의사항
- 그림은 experiment/figures/ 폴더
- 15장 초과 금지
- 반복 설명 금지
- 글보다 그림/표 우선
