# 타겟 정보: EGFR × NSCLC

## 질병: 비소세포폐암 (Non-Small Cell Lung Cancer, NSCLC)
- 전체 폐암의 약 85% 차지
- 5년 생존율: 약 25%
- EGFR 돌연변이는 NSCLC 환자의 약 10-40%에서 발견 (동아시아 인구에서 더 높은 빈도)

## 타겟 단백질: EGFR (Epidermal Growth Factor Receptor)
- **UniProt ID**: P00533
- **ChEMBL Target ID**: CHEMBL203
- **기능**: 수용체 타이로신 키나제 → RAS/MAPK, PI3K/AKT 경로 활성화 → 세포 증식
- **돌연변이 시**: 비정상적 세포 증식 → 종양 형성

## PDB 구조

| PDB ID | 설명 | 해상도 | 리간드 |
|--------|------|--------|--------|
| 1M17 | EGFR 키나제 도메인 + Erlotinib | 2.60 Å | AQ4 (Erlotinib) |

**선정 근거**: 1M17은 wild-type EGFR의 대표적 공결정 구조로, 다수 virtual screening 연구에서 레퍼런스로 사용됨. GNINA 도킹 시 Erlotinib 결합 부위(autobox)로 사용.

## 기존 FDA 승인 EGFR 억제제

| 세대 | 약물 | 작용기전 | IC50 (nM) |
|------|------|---------|-----------|
| 1세대 | Erlotinib | 가역적 TKI | 2 |
| 1세대 | Gefitinib | 가역적 TKI | 33 |
| 2세대 | Afatinib | 비가역적 TKI | 0.5 |
| 3세대 | Osimertinib | T790M 선택적 | 1 |
