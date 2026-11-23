---
name: screening-report
description: BRAF V600E 타겟 피부 흑색종 virtual screening 보고서 작성. experiment의 실측 파일을 읽고 후보 흐름과 근거를 재현한다.
---
저장소 루트 기준으로 작업한다. 먼저 experiment/counts.csv와 audit.csv를 직접 읽고 단계별 수, 최종 후보 수 및 최종 후보 ID를 출력한다. audit.csv의 status=FINAL_CANDIDATE 수와 counts.csv 마지막 단계가 일치하는지 확인한다. 없거나 불일치하면 보고서를 생성하지 말고 파일명을 반환한다.

최종 보고서 요청이면 아래 순서로 experiment 안의 실제 파일만 근거로 작성한다. 파일이 없으면 해당 섹션에서 중단한다.

## 과학적 근거 (탐색 완료, URL 검증 완료)

탐색 경로: 질병·하위유형 › 단백질·역할 › 조절 방향 › 근거 URL › 구조 ID

| 항목 | 내용 | 근거 URL |
|------|------|----------|
| 질병 | 흑색종 (melanoma, MONDO_0005105) | https://platform.opentargets.org/disease/MONDO_0005105/associations |
| 하위유형 | 피부 흑색종 (cutaneous melanoma, MONDO_0005012). BRAF V600E 돌연변이가 환자의 ~50%에서 발견 | https://platform.opentargets.org/disease/MONDO_0005012/associations |
| 단백질·역할 | BRAF (B-Raf proto-oncogene, serine/threonine kinase). cutaneous melanoma 연관 2위 (score 0.701, 3,674개 타겟 중). RAS/MAPK 신호전달 경로에서 MAP2K1을 인산화하여 세포 증식 촉진 | https://platform.opentargets.org/target/ENSG00000157764 |
| 후보 비교 | BRAF(0.701) vs CDKN2A(0.787, 1위). CDKN2A는 종양 억제 유전자로 소분자 억제제 개발 부적합. BRAF는 과활성 키나제로 직접 억제 가능하며 FDA 승인약 3종 존재 → BRAF 선정 | Open Targets 동일 |
| 조절 방향 | BRAF V600E 키나제 억제 — Vemurafenib(PLX4032)은 BRAF V600E를 선택적으로 억제하여 MAPK 신호전달 차단 → 종양 성장 억제 | https://www.cancer.gov/publications/dictionaries/cancer-drug/def/vemurafenib |
| UniProt | P15056 — Homo sapiens, 766 aa, Nucleus/Cell membrane, Serine/threonine-protein kinase | https://www.uniprot.org/uniprotkb/P15056/entry |
| ChEMBL Target | CHEMBL5145 (IC50 데이터 11,017건) | https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL5145/ |
| PDB 구조 | 3OG7 — BRAF V600E kinase domain + PLX4032(Vemurafenib), **V600E mutant**, 2.45 Å, Homo sapiens | https://www.rcsb.org/structure/3OG7 |

참고: EGFR(1M17)은 wild-type이지만, BRAF(3OG7)는 V600E mutant 구조를 사용한다. BRAF V600E 돌연변이 자체가 약물 타겟이므로 mutant 구조가 도킹에 적합하다.

기존 FDA 승인 BRAF 억제제:
| 약물 | 세대 | 작용기전 | 승인 |
|------|------|---------|------|
| Vemurafenib (Zelboraf) | 1세대 | BRAF V600E 선택적 억제 | 2011 |
| Dabrafenib (Tafinlar) | 1세대 | BRAF V600E/K 선택적 억제 | 2013 |
| Encorafenib (Braftovi) | 2세대 | BRAF V600E 선택적 억제 | 2018 |

## 읽어야 할 파일과 용도
| 파일 | 용도 |
|------|------|
| experiment/counts.csv | 단계별 후보 수 (보고서 표·funnel 근거) |
| experiment/audit.csv | 후보별 단계 통과/탈락/이유 (ID 역추적) |
| experiment/05_final_candidates/final_candidates.csv | 최종 후보 ID·점수·물성 |
| experiment/05_final_candidates/generated_molecules.csv | SELFIES 생성 분자 |
| experiment/05_final_candidates/generated_admet_evaluated.csv | 생성 분자 ADMET 재평가 |
| experiment/figures/*.png | 보고서 삽입 그림 (상대경로 명시) |

## 모델 6종
| 모델 | 파일 | 역할 |
|------|------|------|
| RDKit Lipinski/PAINS/Scaffold | experiment/01_preprocessing/ | 물리화학 필터 + Bemis-Murcko scaffold 그룹핑 |
| ChemBERTa (seyonec/ChemBERTa-zinc-base-v1) | experiment/02_model_training/ | BRAF 활성 분류 (Transformer) |
| GNINA (Docker gnina/gnina) | experiment/03_docking/ | PDB 3OG7 기반 구조 도킹 (CNN scoring) |
| ADMET-AI (Chemprop) | experiment/04_admet/ | 약동학 16속성 + 물리화학적 성질 |
| Morgan Fingerprint (ECFP4) | experiment/01_preprocessing/ | Tanimoto 유사도 |
| SELFIES Generator | experiment/05_final_candidates/ | 신규 분자 생성 |

## 종합 점수 가중치
experiment/05_final_candidates/select_final.py에서 확인:
CNN affinity 30% + ChemBERTa 25% + ADMET 25% + Vina 20% (각 정규화)

## 시도와 실패 (보고서에 반드시 포함)
1. DeepDTA(DeepPurpose MPNN_CNN_DAVIS) 시도 → 실험값과 상관계수 r≈0.03, 사실상 무효하여 GNINA로 대체.
2. ChemBERTa 회귀(pChEMBL 연속값 예측) 시도 → CV R²=-18.2, 이진 분류로 전환하여 ROC-AUC 달성.
3. ChEMBL 2,897,819개 전체 처리 시도 → CPU 환경에서 임베딩 추출 수일 소요 예상, drug-like 50,000개 샘플링으로 대응.
4. SELFIES 생성 분자 ADMET 재평가 → ADMET 점수가 시드 대비 전반적으로 낮음. 랜덤 돌연변이의 약물 적합성 한계 확인. experiment/05_final_candidates/generated_admet_evaluated.csv 참조.

## 보고서 구조 (최대 15쪽, 11pt, 한국어, 학술 문체)
1. 서론 (0.5쪽): 질병·타겟 근거 URL 표, 기존 BRAF 억제제 표.
2. 타겟 단백질 (0.5쪽): PDB 3OG7, V600E mutant vs wild-type 설명.
3. 파이프라인 (1쪽): experiment/counts.csv 수치로 funnel 표 작성. 그림 삽입.
4. 모델 선정 (1쪽): 6종 표. CPU+Docker 제약. 실패 기록.
5. 실험 결과 (4쪽): 그림 중심. 각 단계 3줄 이내. 그림마다 원본 상대경로 명시.
6. 최종 후보 (1.5쪽): final_candidates.csv에서 직접 읽은 ID·점수 표. 구조 그림.
7. 결론 (0.5쪽): 한계(실험 검증 부재). 예측은 실험 증명이 아님.

## 중단 조건
- experiment/ 하위 CSV가 없으면 중단
- counts.csv 수와 본문 수가 불일치하면 중단
- 숫자를 임의로 만들지 않는다

## 저장
report/final.docx (최대 15쪽). 저장 후 파일명·쪽수·counts.csv 마지막 수를 다시 출력한다.
문체 수정 요청은 이 지침 수정 후 새 세션에서 재생성한다.
