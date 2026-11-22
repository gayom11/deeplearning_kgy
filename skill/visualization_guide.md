# 시각화 자료 가이드

모든 그림은 `experiment/figures/` 폴더에 위치.

| 파일명 | 보고서 섹션 | 설명 |
|--------|-----------|------|
| fig01_preprocessing.png | 5.1 | 학습 데이터 vs 스크리닝 라이브러리 — MW, logP, HBD, HBA, TPSA, 유사도 분포 비교 |
| fig02_chemberta_training.png | 5.2 | ChemBERTa 학습 결과 — t-SNE, CV ROC-AUC(0.897), OOF 활성 분류, pChEMBL vs 예측 |
| fig03_screening_results.png | 5.3 | Virtual Screening — 활성 확률 분포, Top 1000 drug-likeness, 구조 신규성 |
| fig04_admet_heatmap.png | 5.5 | 상위 20개 × 16 ADMET 속성 히트맵 |
| fig05_admet_radar.png | 5.5 | 상위 5개 후보 ADMET radar chart |
| fig06_admet_distributions.png | 5.5 | HIA, Bioavailability, Caco-2, hERG, AMES, DILI 분포 |
| fig07_docking_results.png | 5.4 | GNINA — CNN affinity 분포, Vina score 분포, Vina vs CNN scatter |
| fig08_candidate_structures.png | 6 | 최종 10개 후보 + 기존 4개 약물 2D 분자 구조 |
| fig09_scatter_analysis.png | 6 | 다차원 비교 — MW-logP, ChemBERTa-GNINA, ADMET-GNINA |
| fig10_pipeline_funnel.png | 3 | 파이프라인 funnel (50,000 → 10) |
| fig11_generated_molecules.png | 5.6 | 생성 분자 화학 공간, 유사도 분포, novelty 파이차트 |
| fig12_novel_generated.png | 5.6 | 가장 novel한 생성 분자 10개 2D 구조 |

## 핵심 해석 포인트

### fig02 (ChemBERTa)
- ROC-AUC 0.897: 양호한 분류 성능
- t-SNE에서 Active(초록)/Inactive(빨강)이 분리됨 → 임베딩이 활성 패턴 포착

### fig07 (GNINA)
- CNN affinity는 pKd 단위 (높을수록 강한 결합)
- Vina score는 kcal/mol (낮을수록 강한 결합)
- 두 점수 간 양의 상관관계 확인

### fig10 (Funnel)
- 50,000 → 47,050 → 1,000 → 300 → 297 → 10 단계적 축소
- 각 단계에서 사용된 모델과 필터 기준 명시
