# EGFR 타겟 비소세포폐암 Virtual Screening 보고서 작성 스킬

## 역할
experiment/ 폴더의 실제 계산 결과 파일만 읽고 보고서를 작성한다.
파일이 없으면 생성을 중단한다. 숫자를 임의로 만들지 않는다.

## 중단 조건
- experiment/ 하위 CSV 파일이 존재하지 않으면 해당 섹션 작성을 중단하고 "데이터 누락" 표시
- audit.csv의 단계별 수와 본문의 수가 다르면 중단

## 검산 규칙
- 보고서에 쓰는 모든 수치는 experiment/audit.csv 또는 해당 CSV의 행 수와 일치해야 한다
- 최종 후보 ID는 experiment/05_final_candidates/final_candidates.csv의 molecule_chembl_id를 그대로 사용

## 보고서 형식
- 최대 15장, 본문 11pt, 한국어, 학술 문체 ("~하였다")
- 한 문단 최대 3줄, 글보다 그림·표 우선
- 근거 이상의 결론 금지, 한계 문장 필수
- 실패한 시도도 기록 (무엇을 시도했고 왜 실패했는지)
- report/final.docx로 저장

## 과학적 근거 (URL 직접 확인 가능)
- 질병: EGFR 관련 비소세포폐암 (NSCLC)
  - Open Targets: https://platform.opentargets.org/disease/EFO_0003060/associations
  - NCI Drug Dictionary (Erlotinib): https://www.cancer.gov/publications/dictionaries/cancer-drug/def/erlotinib-hydrochloride
- 타겟 단백질: EGFR (ErbB1)
  - UniProt P00533: https://www.uniprot.org/uniprot/P00533
  - ChEMBL Target CHEMBL203: https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL203/
- PDB 구조: 1M17 (EGFR kinase domain + Erlotinib, wild-type, 2.60 Å)
  - RCSB PDB: https://www.rcsb.org/structure/1M17
- 선정 근거: EGFR은 NSCLC 환자의 10-40%에서 돌연변이 발견. 1M17은 wild-type EGFR의 대표적 공결정 구조.

## 기존 FDA 승인 EGFR 억제제
| 세대 | 약물 | IC50 (nM) |
|------|------|-----------|
| 1세대 | Erlotinib | 2 |
| 1세대 | Gefitinib | 33 |
| 2세대 | Afatinib | 0.5 |
| 3세대 | Osimertinib | 1 |

## 파이프라인 흐름 (experiment/audit.csv로 검산)
```
ChEMBL 2,897,819개 (chembl_37_id_smiles.txt)
  ↓ Drug-like 필터 + 샘플링 → 50,000개
  ↓ [Step 1] Lipinski Ro5 + PAINS → 47,050개
  ↓ [Step 1b] Scaffold 그룹핑 (Bemis-Murcko, 그룹당 ≤3) → 22,726개
  ↓ [Step 2] ChemBERTa 활성 예측 (상위 1,000) → 1,000개
  ↓ [Step 3] ADMET-AI 약동학 평가 (상위 300) → 300개
  ↓ [Step 4] GNINA 도킹 (PDB 1M17) → 297개
  ↓ [Step 5] 종합 점수 + 다양성 필터 → 10개
  ↓ [Step 6] SELFIES 생성 → 778개 신규 drug-like (73% novel)
```

## 모델 6종 (experiment/ 하위 스크립트 참조)
| 모델 | 유형 | 역할 | 평가 축 |
|------|------|------|---------|
| RDKit | 규칙 기반 | Lipinski/PAINS/Scaffold | 물리화학 |
| ChemBERTa | Transformer | EGFR 활성 분류 (ROC-AUC 0.897) | 약효 |
| GNINA | CNN 도킹 (Docker) | 단백질-리간드 결합 (Vina+CNN) | 약효 |
| ADMET-AI | Chemprop 앙상블 | 16속성 약동학 평가 | 물리화학+ADMET |
| Morgan FP | ECFP4 | Tanimoto 유사도 | 약효 |
| SELFIES | 문법 기반 생성 | 신규 분자 생성 | 탐색 확장 |

## 계산 환경 제약과 선택
- CPU: AMD Ryzen 5 9600X (GPU 미보유)
- Docker Desktop (WSL2) → GNINA 실행
- 시도했으나 포기한 모델: DeepDTA (MPNN_CNN_DAVIS) → 실험값과 상관계수 r≈0.03, 사실상 무효하여 제외

## 종합 점수 가중치
최종 점수 = 0.30×CNN affinity + 0.25×ChemBERTa + 0.25×ADMET + 0.20×Vina (각 정규화)

## 최종 후보 실제 데이터 (이 숫자를 그대로 사용할 것)
| ChEMBL ID | ChemBERTa | CNN affinity | Vina score | ADMET | 최종 점수 |
|-----------|-----------|-------------|------------|-------|----------|
| CHEMBL57439 | 0.878 | 7.711 | -10.159 | 0.259 | 0.761 |
| CHEMBL26050 | 0.949 | 7.875 | -8.360 | 0.179 | 0.761 |
| CHEMBL43181 | 0.888 | 7.763 | -10.004 | 0.233 | 0.761 |
| CHEMBL24625 | 0.940 | 7.792 | -9.128 | 0.192 | 0.760 |
| CHEMBL11779 | 0.876 | 7.916 | -9.032 | 0.171 | 0.750 |
| CHEMBL295433 | 0.859 | 7.589 | -9.393 | 0.253 | 0.737 |
| CHEMBL280842 | 0.890 | 7.538 | -9.013 | 0.250 | 0.737 |
| CHEMBL15431 | 0.823 | 7.588 | -10.333 | 0.258 | 0.737 |
| CHEMBL41383 | 0.925 | 7.449 | -8.997 | 0.239 | 0.734 |
| CHEMBL295106 | 0.843 | 7.592 | -8.930 | 0.259 | 0.732 |

## 시도와 실패 기록 (보고서에 반드시 포함)
1. **DeepDTA 시도 → 실패**: DeepPurpose MPNN_CNN_DAVIS 모델로 결합 친화도 예측 시도. 실험값(pChEMBL)과 예측값(pKd) 상관계수 r≈0.03으로 사실상 무상관. 원인: Davis kinase dataset과 EGFR IC50 데이터 간 도메인 차이, Kd와 IC50 메트릭 불일치. → GNINA 도킹으로 대체.
2. **ChemBERTa 회귀 시도 → 실패**: pChEMBL 연속값 예측에 ChemBERTa + RandomForest 회귀 적용. CV R² = -18.2 (평균 예측보다 나쁨). → 이진 분류(Active/Inactive)로 전환하여 ROC-AUC 0.897 달성.
3. **전체 ChEMBL 처리 → 축소**: 2,897,819개 전체 처리는 CPU 환경에서 ChemBERTa 임베딩 추출에 수일 소요 예상. → Drug-like 50,000개 샘플링 + scaffold 선별로 대응.

## 그림 목록 (experiment/figures/)
| 파일 | 섹션 | 내용 |
|------|------|------|
| fig01_preprocessing.png | 5.1 | 학습/스크리닝 속성 분포 비교 |
| fig01b_scaffold_selection.png | 5.1 | Scaffold 그룹핑 결과 |
| fig02_chemberta_training.png | 5.2 | t-SNE, ROC-AUC, OOF 분류 |
| fig03_screening_results.png | 5.3 | 활성 확률 분포, 상위 hit |
| fig04_admet_heatmap.png | 5.5 | ADMET 16속성 히트맵 |
| fig05_admet_radar.png | 5.5 | 상위 5개 ADMET radar |
| fig07_docking_results.png | 5.4 | CNN affinity, Vina 분포 |
| fig08_candidate_structures.png | 6 | 최종 후보 2D 구조 |
| fig09_scatter_analysis.png | 6 | 다차원 비교 scatter |
| fig10_pipeline_funnel.png | 3 | 파이프라인 funnel |
| fig11_generated_molecules.png | 5.6 | 생성 분자 화학 공간 |

## 보고서 구조
1. 서론 (0.5장): 질병·타겟 근거 URL, 기존 약물 표
2. 타겟 단백질 (0.5장): PDB 1M17, wild-type 설명
3. 파이프라인 (1장): funnel 그림 + 단계별 수 (audit.csv 기준)
4. 모델 선정 (1장): 6종 표 + 환경 제약 + DeepDTA 실패 기록
5. 실험 결과 (4장): 그림 중심, 각 단계 3줄 이내
6. 최종 후보 (1.5장): 실제 데이터 표 + 구조 그림
7. 결론 (0.5장): 성과 + 한계(실험 검증 부재) + 향후 과제
