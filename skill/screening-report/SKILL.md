---
name: screening-report
description: EGFR 타겟 비소세포폐암 virtual screening 보고서 작성. experiment의 실측 파일을 읽고 후보 흐름과 근거를 재현한다.
---
저장소 루트 기준으로 작업한다. 먼저 experiment/counts.csv와 audit.csv를 직접 읽고 단계별 수, 최종 후보 수 및 최종 후보 ID를 출력한다. audit.csv의 status=FINAL_CANDIDATE 수와 counts.csv 마지막 단계가 일치하는지 확인한다. 없거나 불일치하면 보고서를 생성하지 말고 파일명을 반환한다.

최종 보고서 요청이면 아래 순서로 experiment 안의 실제 파일만 근거로 작성한다. 파일이 없으면 해당 섹션에서 중단한다.

## 과학적 근거 (탐색 결과, URL 직접 검증 완료)

탐색 경로: 질병·하위유형 › 단백질·역할 › 조절 방향 › 근거 URL › 구조 ID

| 항목 | 내용 | 근거 URL |
|------|------|----------|
| 질병·하위유형 | 비소세포폐암 (non-small cell lung carcinoma, MONDO_0005233) | https://platform.opentargets.org/disease/MONDO_0005233/associations |
| 단백질·역할 | EGFR (Epidermal Growth Factor Receptor) — 수용체 타이로신 키나제. Open Targets 연관성 1위 (score 0.888, 5,860개 타겟 중) | https://platform.opentargets.org/target/ENSG00000146648 |
| 조절 방향 | EGFR 억제 (저해) — Erlotinib은 ATP와 경쟁적으로 EGFR 타이로신 키나제를 억제하여 종양 성장 억제 | https://www.cancer.gov/publications/dictionaries/cancer-drug/def/erlotinib-hydrochloride |
| UniProt | P00533 — Homo sapiens, 1,210 aa, Cell membrane, Receptor tyrosine kinase | https://www.uniprot.org/uniprotkb/P00533/entry |
| ChEMBL Target | CHEMBL203 | https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL203/ |
| PDB 구조 | 1M17 — EGFR kinase domain + Erlotinib(AQ4), wild-type(변이 없음), 2.60 Å, Homo sapiens | https://www.rcsb.org/structure/1M17 |

후보 단백질 비교: EGFR(score 0.888) vs KRAS(score 0.834). EGFR은 타이로신 키나제 억제제로 직접 약물 타겟이 가능하며 기존 FDA 승인야이 4종 존재. KRAS는 이후 신호전달 단백질로 직접 결합 표적화가 어려워 EGFR을 선정함.

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
4. SELFIES 생성 분자 778개 ADMET 재평가 → ADMET 점수 최대 0.13으로 시드 분자(0.17~0.26)보다 전반적으로 낮음. SELFIES 랜덤 돌연변이가 약물 적합성을 유지하지 못하는 한계 확인. experiment/05_final_candidates/generated_admet_evaluated.csv 참조.

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
