# 시각화 자료 가이드

보고서에 삽입할 그림 목록과 해석 가이드. 모든 그림은 `experiment/figures/` 폴더에 위치.

## 그림 목록

| 파일명 | 보고서 섹션 | 설명 |
|--------|-----------|------|
| fig01_lipinski_filtering.png | 5.1 | Lipinski 필터링 전후 6개 속성 분포 비교 (Before: 빨강, After: 초록) |
| fig02_similarity_distribution.png | 5.2 | 기존 EGFR 약물 대비 Tanimoto 유사도 분포 (빨간 점선: 0.4 기준) |
| fig03_binding_affinity.png | 5.3 | DeepDTA 예측 pKd 분포, 실험값 대비 scatter, Top 50 bar chart |
| fig04_chemberta_analysis.png | 5.4 | t-SNE 임베딩, CV ROC-AUC, 회귀 scatter, 활성 확률 분포 |
| fig05_admet_heatmap.png | 5.5 | 상위 20개 화합물 × 16 ADMET 속성 히트맵 |
| fig06_admet_radar.png | 5.5 | 상위 5개 후보의 ADMET radar chart |
| fig07_admet_distributions.png | 5.5 | 6개 주요 ADMET 속성 분포 |
| fig08_pipeline_funnel.png | 3 | 파이프라인 필터링 funnel (10,005 → 10) |
| fig09_candidate_structures.png | 6 | 최종 10개 후보 + 기존 4개 약물의 2D 분자 구조 |
| fig10_comparison_table.png | 6 | 후보 vs 기존 약물 물리화학적 성질 비교 표 |
| fig11_scatter_analysis.png | 6 | 다차원 비교 scatter (MW-logP, 활성-ADMET, 활성-유사도) |
| fig12_generated_molecules.png | 5.6 | 생성 분자 화학 공간, 유사도 분포, novelty 파이차트 |
| fig13_novel_generated.png | 5.6 | 가장 novel한 생성 분자 10개의 2D 구조 |

## 그림 해석 핵심 포인트

### fig01 (Lipinski 필터링)
- MW 500 이상, logP 5 이상 화합물이 필터에 의해 제거됨
- 필터링 후 drug-like 화합물 공간에 집중

### fig04 (ChemBERTa)
- t-SNE에서 pChEMBL 값에 따른 그라데이션 확인 → 임베딩이 활성 정보를 포착
- ROC-AUC 0.831 → 양호한 분류 성능
- 회귀 R²가 음수 → 연속값 예측보다 이진 분류가 적합

### fig06 (ADMET Radar)
- 각 꼭짓점 점수가 1에 가까울수록 약물 적합성 높음
- BBB는 항암제에 불필요 (뇌 투과 불필요)
- hERG Safe, Non-Mutagen이 높은 후보가 안전

### fig08 (Pipeline Funnel)
- 10,005개에서 10개로 단계적 축소
- 각 단계에서 사용된 모델과 필터 기준 명시

### fig12 (분자 생성)
- 66% novel 분자 → 기존 화학 공간 외 새로운 영역 탐색 성공
- 시드 분자(빨간 별)를 중심으로 화학 공간이 확장됨
