---
name: screening-report
description: EGFR 타겟 비소세포폐암 virtual screening 보고서 작성. experiment의 실측 파일을 읽고 후보 흐름과 근거를 재현한다.
---
저장소 루트 기준으로 작업한다. 먼저 experiment/counts.csv와 audit.csv를 직접 읽고 단계별 수, 최종 후보 수 및 최종 후보 ID를 출력한다. audit.csv의 status=FINAL_CANDIDATE 수와 counts.csv 마지막 단계가 일치하는지 확인한다. 없거나 불일치하면 보고서를 생성하지 말고 파일명을 반환한다.

최종 보고서 요청이면 아래 순서로 experiment 안의 실제 파일만 근거로 작성한다. 파일이 없으면 해당 섹션에서 중단한다.

## 과학적 근거 (Solar가 직접 URL을 열어 확인)
- 질병: EGFR 관련 비소세포폐암. Open Targets https://platform.opentargets.org/disease/EFO_0003060/associations
- 약물 근거: NCI Erlotinib https://www.cancer.gov/publications/dictionaries/cancer-drug/def/erlotinib-hydrochloride
- 타겟: UniProt P00533 https://www.uniprot.org/uniprot/P00533 / ChEMBL CHEMBL203
- 구조: PDB 1M17 (wild-type EGFR + Erlotinib, 2.60 Å) https://www.rcsb.org/structure/1M17

## 읽어야 할 파일과 용도
| 파일 | 용도 |
|------|------|
| experiment/counts.csv | 단계별 후보 수 (보고서 표·funnel 근거) |
| experiment/audit.csv | 후보별 단계 통과/탈락/이유 (ID 역추적) |
| experiment/05_final_candidates/final_candidates.csv | 최종 10개 후보 ID·점수·물성 |
| experiment/05_final_candidates/generated_molecules.csv | SELFIES 생성 분자 |
| experiment/figures/*.png | 보고서 삽입 그림 (상대경로 명시) |

## 모델 6종 (experiment/ 하위 스크립트 참조)
| 모델 | 파일 | 역할 |
|------|------|------|
| RDKit Lipinski/PAINS/Scaffold | experiment/01_preprocessing/ | 물리화학 필터 |
| ChemBERTa (seyonec/ChemBERTa-zinc-base-v1) | experiment/02_model_training/chemberta_screen.py | EGFR 활성 분류 |
| GNINA (Docker gnina/gnina) | experiment/03_docking/ | 구조 기반 도킹 |
| ADMET-AI (Chemprop) | experiment/04_admet/admet_screen.py | 약동학 16속성 |
| Morgan Fingerprint (ECFP4) | experiment/01_preprocessing/preprocess.py | Tanimoto 유사도 |
| SELFIES Generator | experiment/05_final_candidates/select_final.py | 신규 분자 생성 |

## 종합 점수 가중치
experiment/05_final_candidates/select_final.py에서 확인:
CNN affinity 30% + ChemBERTa 25% + ADMET 25% + Vina 20% (각 정규화)

## 시도와 실패 (보고서에 반드시 포함)
1. DeepDTA(DeepPurpose MPNN_CNN_DAVIS) 시도 → 실험값과 상관계수 r≈0.03, 사실상 무효하여 GNINA로 대체.
2. ChemBERTa 회귀(pChEMBL 연속값 예측) 시도 → CV R²=-18.2, 이진 분류로 전환하여 ROC-AUC 0.897 달성.
3. ChEMBL 2,897,819개 전체 처리 시도 → CPU 환경에서 임베딩 추출 수일 소요 예상, drug-like 50,000개 샘플링으로 대응.

## 보고서 구조 (최대 15쪽, 11pt, 한국어, 학술 문체)
1. 서론 (0.5쪽): 질병·타겟 근거 URL 표, 기존 EGFR 억제제 표. 한 문단 3줄 이내.
2. 타겟 단백질 (0.5쪽): PDB 1M17 wild-type. 표 1개.
3. 파이프라인 (1쪽): experiment/counts.csv 수치로 funnel 표 작성. experiment/figures/fig10_pipeline_funnel.png 삽입.
4. 모델 선정 (1쪽): 6종 표. CPU+Docker 제약. DeepDTA 실패 1줄.
5. 실험 결과 (4쪽): 그림 중심. 각 단계 3줄 이내. 그림마다 원본 상대경로 명시.
   - 5.1 전처리+Scaffold: fig01_preprocessing.png, fig01b_scaffold_selection.png
   - 5.2 ChemBERTa: fig02_chemberta_training.png
   - 5.3 Screening: fig03_screening_results.png
   - 5.4 GNINA: fig07_docking_results.png
   - 5.5 ADMET: fig04_admet_heatmap.png
   - 5.6 생성: fig11_generated_molecules.png
6. 최종 후보 (1.5쪽): experiment/05_final_candidates/final_candidates.csv에서 직접 읽은 ID·점수 표. fig08, fig09 삽입.
7. 결론 (0.5쪽): 한계(실험 검증 부재). 예측은 실험 증명이 아님. 실패를 낮은 점수로 바꾸지 않는다.

## 중단 조건
- experiment/ 하위 CSV가 없으면 중단
- counts.csv 수와 본문 수가 불일치하면 중단
- 숫자를 임의로 만들지 않는다

## 저장
report/final.docx (최대 15쪽). 저장 후 파일명·쪽수·counts.csv 마지막 수를 다시 출력한다.
문체 수정 요청은 이 지침 수정 후 새 세션에서 재생성한다.
