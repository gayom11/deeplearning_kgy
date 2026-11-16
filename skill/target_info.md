# 타겟 정보: EGFR × NSCLC

## 질병: 비소세포폐암 (Non-Small Cell Lung Cancer, NSCLC)
- 전체 폐암의 약 85% 차지
- 주요 아형: 선암(adenocarcinoma), 편평세포암(squamous cell carcinoma), 대세포암(large cell carcinoma)
- EGFR 돌연변이는 NSCLC 환자의 약 10-40%에서 발견 (동아시아 인구에서 더 높은 빈도)
- 5년 생존율: 약 25% (전체), EGFR TKI 치료 시 60% 이상 반응률

## 타겟 단백질: EGFR (Epidermal Growth Factor Receptor)
- **UniProt ID**: P00533
- **ChEMBL Target ID**: CHEMBL203
- **유전자**: ERBB1 (ErbB family)
- **기능**: 수용체 타이로신 키나제 (Receptor Tyrosine Kinase)
- **위치**: 세포막 관통 단백질
- **신호전달**: EGFR → RAS/MAPK, PI3K/AKT 경로 활성화 → 세포 증식, 생존, 이동 촉진
- **돌연변이 시**: 비정상적 세포 증식 → 종양 형성

## PDB 구조
| PDB ID | 설명 | 해상도 | 리간드 | 돌연변이 |
|--------|------|--------|--------|---------|
| 1M17 | EGFR 키나제 도메인 + Erlotinib | 2.60 Å | AQ4 (Erlotinib) | 없음 (wild-type) |
| 4HJO | EGFR 키나제 도메인 (inactive) + Erlotinib | 2.75 Å | AQ4 | 있음 |

**선정 근거**: 1M17은 wild-type EGFR의 가장 대표적인 공결정 구조로, 다수의 virtual screening 연구에서 레퍼런스로 사용됨.

## 기존 FDA 승인 EGFR 억제제
| 세대 | 약물 | 작용기전 | 승인 | IC50 (nM) |
|------|------|---------|------|-----------|
| 1세대 | Erlotinib (Tarceva) | 가역적 TKI | 2004 | 2 |
| 1세대 | Gefitinib (Iressa) | 가역적 TKI | 2003 | 33 |
| 2세대 | Afatinib (Gilotrif) | 비가역적 TKI | 2013 | 0.5 |
| 3세대 | Osimertinib (Tagrisso) | T790M 선택적 | 2015 | 1 |

## EGFR 키나제 도메인 시퀀스 (UniProt P00533, 696-1022)
```
FKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDN
PHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDR
RLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIY
THQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCW
MIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVV
DADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDR
```
