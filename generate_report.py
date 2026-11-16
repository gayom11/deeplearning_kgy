"""
Solar Pro4 API를 사용하여 skill 폴더 기반 보고서 생성 + python-docx로 최종 docx 작성
"""
import os
import sys
import json
import time
import requests
from docx import Document
from docx.shared import Inches, Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn

API_KEY = "up_r2f9Y6uBUMZTg7S885tCGDdnGZUxO"
API_URL = "https://api.upstage.ai/v1/chat/completions"
MODEL = "solar-pro4"

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SKILL_DIR = os.path.join(PROJECT_DIR, "skill")
FIGURES_DIR = os.path.join(PROJECT_DIR, "experiment", "figures")
REPORT_DIR = os.path.join(PROJECT_DIR, "report")
os.makedirs(REPORT_DIR, exist_ok=True)


def read_skill_files():
    """skill 폴더의 모든 md 파일 읽기"""
    content = {}
    for fname in sorted(os.listdir(SKILL_DIR)):
        if fname.endswith(".md"):
            with open(os.path.join(SKILL_DIR, fname), "r", encoding="utf-8") as f:
                content[fname] = f.read()
    return content


def call_solar(system_prompt, user_prompt, max_tokens=4000):
    """Solar Pro4 API 호출"""
    headers = {"Authorization": f"Bearer {API_KEY}", "Content-Type": "application/json"}
    payload = {
        "model": MODEL,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        "max_tokens": max_tokens,
        "temperature": 0.3,
    }

    for attempt in range(3):
        try:
            resp = requests.post(API_URL, headers=headers, json=payload, timeout=120)
            resp.raise_for_status()
            data = resp.json()
            return data["choices"][0]["message"]["content"].strip()
        except Exception as e:
            print(f"  API error (attempt {attempt+1}): {e}")
            time.sleep(5)

    return None


def generate_sections(skills):
    """각 섹션별로 Solar Pro4 호출하여 보고서 텍스트 생성"""
    # 전체 skill 내용을 하나의 컨텍스트로
    full_skill = "\n\n---\n\n".join(f"## {k}\n{v}" for k, v in skills.items())

    system_prompt = """당신은 약학/신약개발 분야의 학술 보고서를 작성하는 전문가입니다.
규칙:
- 한국어로 작성
- 학술 보고서 문체 ("~하였다", "~이다")
- 간결하게 서술, 한 문단 최대 4줄
- 불필요한 수식어, 반복 설명 금지
- 표와 그림 참조를 적극 활용
- AI 특유의 과잉 설명 하지 말 것"""

    sections = {}

    # 섹션 1: 서론
    print("Generating: 1. 서론...")
    sections["intro"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 보고서의 "1. 서론" 섹션을 작성하세요.

내용:
- 비소세포폐암(NSCLC) 개요 (폐암의 85%, 5년 생존율 약 25%)
- EGFR 타겟 선정 근거 (NSCLC 환자 10-40%에서 EGFR 돌연변이)
- 기존 EGFR 억제제 (Erlotinib, Gefitinib, Afatinib, Osimertinib)
- 연구 목적: 딥러닝 기반 virtual screening으로 새로운 EGFR 억제제 후보 발굴

분량: 10-15줄 이내. 표 1개 포함 (기존 EGFR 억제제 비교).

{full_skill}""", max_tokens=1500)

    # 섹션 2: 타겟 단백질
    print("Generating: 2. 타겟 단백질...")
    sections["target"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "2. 타겟 단백질 정보" 섹션을 작성하세요.

내용:
- EGFR 단백질 기능, ChEMBL ID, PDB ID 정보
- PDB 1M17 선정 근거
- 표로 정리 (PDB ID, 해상도, 리간드, 돌연변이 여부)

분량: 8-10줄 + 표 1개.

{full_skill}""", max_tokens=1200)

    # 섹션 3: 파이프라인
    print("Generating: 3. Virtual Screening 파이프라인...")
    sections["pipeline"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "3. Virtual Screening 파이프라인" 섹션을 작성하세요.

내용:
- 전체 파이프라인 개요 (그림 fig08_pipeline_funnel.png 참조 언급)
- 각 단계 한두 줄 설명:
  1. ChEMBL 데이터 수집 (19,719건 → 10,005개)
  2. Lipinski + PAINS 필터 (→ 6,608개)
  3. DeepDTA 결합 친화도 예측
  4. ChemBERTa 활성 분류
  5. ADMET 예측
  6. 종합 랭킹 → 최종 10개
  7. SELFIES 분자 생성 → 718개

분량: 15줄 이내.

{full_skill}""", max_tokens=1500)

    # 섹션 4: 모델 선정
    print("Generating: 4. 모델 선정...")
    sections["models"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "4. 모델 선정 및 설명" 섹션을 작성하세요.

내용:
- 사용 모델 6종을 표로 정리 (모델명, 유형, 입력, 출력, 선정 근거)
- 계산 자원 제약 (GPU 미보유, CPU만 사용) 설명
- GNINA 대신 DeepDTA 선택 근거

분량: 표 1개 + 설명 10줄 이내.

{full_skill}""", max_tokens=1500)

    # 섹션 5.1-5.2: 물리화학 + 유사도
    print("Generating: 5.1-5.2 필터링 & 유사도...")
    sections["results_filter"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "5.1 물리화학적 필터링" 및 "5.2 분자 유사도 분석" 섹션을 작성하세요.

5.1:
- Lipinski Rule of Five 기준 (MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10)
- PAINS 필터 적용
- 10,005 → 6,608개 (34% 제거)
- 그림 참조: fig01_lipinski_filtering.png

5.2:
- Morgan Fingerprint (ECFP4, radius=2, 2048-bit) 사용
- Erlotinib, Gefitinib, Osimertinib, Afatinib 대비 Tanimoto 유사도
- 그림 참조: fig02_similarity_distribution.png

분량: 각 5-7줄.

{full_skill}""", max_tokens=1200)

    # 섹션 5.3: DeepDTA
    print("Generating: 5.3 결합 친화도...")
    sections["results_binding"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "5.3 결합 친화도 예측 (DeepDTA)" 섹션을 작성하세요.

내용:
- DeepPurpose MPNN_CNN 모델 (Davis kinase dataset 사전학습)
- EGFR kinase domain 시퀀스 + SMILES → pKd 예측
- 실험값(pChEMBL)과 예측값(pKd) 상관관계: r = 0.032
- 상관관계가 낮은 이유: IC50와 Kd는 다른 측정, Davis 데이터셋 도메인 차이
- 그림 참조: fig03_binding_affinity.png

분량: 7-10줄.

{full_skill}""", max_tokens=1200)

    # 섹션 5.4: ChemBERTa
    print("Generating: 5.4 ChemBERTa...")
    sections["results_chemberta"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "5.4 ChemBERTa 활성 분류" 섹션을 작성하세요.

내용:
- seyonec/ChemBERTa-zinc-base-v1 (RoBERTa 아키텍처, ZINC 사전학습)
- 768차원 CLS 임베딩 → PCA(50차원) → RandomForest
- 분류 ROC-AUC: 0.831 ± 0.008 (5-fold CV)
- 활성 기준: pChEMBL ≥ 7.0 (IC50 ≤ 100 nM)
- t-SNE에서 활성/비활성 분자 클러스터 확인
- 그림 참조: fig04_chemberta_analysis.png

분량: 7-10줄.

{full_skill}""", max_tokens=1200)

    # 섹션 5.5: ADMET
    print("Generating: 5.5 ADMET...")
    sections["results_admet"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "5.5 ADMET 프로파일" 섹션을 작성하세요.

내용:
- ADMET-AI (Chemprop 앙상블) 사용
- 16종 ADMET 속성 예측 (HIA, Caco-2, BBB, hERG, AMES, DILI 등)
- 흡수/분포/대사/배설/독성 각 카테고리 설명
- ADMET 종합 점수로 상위 50% 선별
- 그림 참조: fig05, fig06, fig07

분량: 7-10줄.

{full_skill}""", max_tokens=1200)

    # 섹션 5.6: 분자 생성
    print("Generating: 5.6 분자 생성...")
    sections["results_generation"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "5.6 분자 구조 생성 및 화학 공간 탐색" 섹션을 작성하세요.

내용:
- SELFIES 기반 돌연변이 (replace/insert/delete, 1-3회)
- 스캐폴드 호핑 (할로겐 교체, 작용기 변환)
- Murcko 스캐폴드 분해
- 10개 시드 → 754개 생성 → 718개 drug-like (66% novel, Tanimoto < 0.4)
- REINVENT의 경량 대안으로서 GPU 없이 화학 공간 확장
- 그림 참조: fig12_generated_molecules.png, fig13_novel_generated.png

분량: 7-10줄.

{full_skill}""", max_tokens=1200)

    # 섹션 6: 최종 후보
    print("Generating: 6. 최종 후보 물질...")
    sections["candidates"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "6. 최종 후보 물질" 섹션을 작성하세요.

내용:
- 종합 점수 가중치: 활성 0.30, ChemBERTa 0.20, DeepDTA 0.15, ADMET 0.25, 신규성 0.10
- 화학적 다양성 보장 (Tanimoto > 0.8 중복 제거)
- 최종 10개 후보 선정
- 그림 참조: fig09 (분자 구조), fig10 (비교 표), fig11 (scatter)
- 각 후보의 특징 간략히 (quinazoline 계열 scaffold 등)
- 기존 약물 대비 구조적 차이점

분량: 10-12줄 + 표 1개 (상위 5개 후보 ID, pChEMBL, MW, logP, 최종점수).

{full_skill}""", max_tokens=1500)

    # 섹션 7: 결론
    print("Generating: 7. 결론...")
    sections["conclusion"] = call_solar(system_prompt, f"""아래 skill 자료를 참고하여 "7. 결론" 섹션을 작성하세요.

내용:
- 6개 딥러닝/ML 모델을 활용한 체계적 virtual screening 수행
- CPU 환경에서 실현 가능한 파이프라인 구축
- SELFIES 기반 718개 신규 분자 탐색
- 한계: 실험적 검증 부재, 도킹 시뮬레이션 미수행, DeepDTA 상관관계
- 향후 과제: GNINA 도킹, MD 시뮬레이션, 실험적 IC50 측정

분량: 8줄 이내.

{full_skill}""", max_tokens=1000)

    return sections


def build_docx(sections):
    """python-docx로 보고서 docx 생성"""
    doc = Document()

    # 기본 스타일 설정
    style = doc.styles["Normal"]
    font = style.font
    font.name = "맑은 고딕"
    font.size = Pt(11)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), "맑은 고딕")

    # 페이지 여백
    for section in doc.sections:
        section.top_margin = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin = Cm(2.5)
        section.right_margin = Cm(2.5)

    # 제목
    title = doc.add_heading("딥러닝 기반 EGFR 타겟 비소세포폐암\n신약 후보 물질 Virtual Screening", level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in title.runs:
        run.font.size = Pt(18)

    # 저자 정보
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run("Deep Learning Final Project")
    run.font.size = Pt(11)
    run.font.color.rgb = RGBColor(100, 100, 100)

    doc.add_page_break()

    # 섹션 매핑
    section_config = [
        ("1. 서론", "intro", []),
        ("2. 타겟 단백질 정보", "target", []),
        ("3. Virtual Screening 파이프라인", "pipeline", ["fig08_pipeline_funnel.png"]),
        ("4. 모델 선정 및 설명", "models", []),
        ("5. 실험 결과", None, []),
        ("5.1 물리화학적 필터링 및 분자 유사도", "results_filter",
         ["fig01_lipinski_filtering.png", "fig02_similarity_distribution.png"]),
        ("5.2 결합 친화도 예측 (DeepDTA)", "results_binding", ["fig03_binding_affinity.png"]),
        ("5.3 ChemBERTa 활성 분류", "results_chemberta", ["fig04_chemberta_analysis.png"]),
        ("5.4 ADMET 프로파일", "results_admet",
         ["fig05_admet_heatmap.png", "fig06_admet_radar.png"]),
        ("5.5 분자 구조 생성 및 화학 공간 탐색", "results_generation",
         ["fig12_generated_molecules.png", "fig13_novel_generated.png"]),
        ("6. 최종 후보 물질", "candidates",
         ["fig09_candidate_structures.png", "fig10_comparison_table.png", "fig11_scatter_analysis.png"]),
        ("7. 결론", "conclusion", []),
    ]

    for heading_text, key, figures in section_config:
        # 헤딩
        level = 1 if not heading_text.startswith("5.") else 2
        if heading_text == "5. 실험 결과":
            level = 1
        h = doc.add_heading(heading_text, level=level)
        for run in h.runs:
            run.font.size = Pt(13) if level == 1 else Pt(12)

        # 본문 텍스트
        if key and key in sections and sections[key]:
            text = sections[key]
            # 마크다운 테이블을 일반 텍스트로 변환 (docx에서는 별도 처리)
            for line in text.split("\n"):
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    # 하위 헤딩은 bold 문단으로
                    clean = line.lstrip("#").strip()
                    p = doc.add_paragraph()
                    run = p.add_run(clean)
                    run.bold = True
                    run.font.size = Pt(11)
                elif line.startswith("|") and "---" not in line:
                    # 테이블 행 → 그냥 텍스트로
                    p = doc.add_paragraph(line)
                    p.style.font.size = Pt(10)
                elif line.startswith("---"):
                    continue
                else:
                    p = doc.add_paragraph(line)
                    p.style.font.size = Pt(11)

        # 그림 삽입
        for fig_name in figures:
            fig_path = os.path.join(FIGURES_DIR, fig_name)
            if os.path.exists(fig_path):
                doc.add_paragraph()  # 여백
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                run = p.add_run()
                run.add_picture(fig_path, width=Inches(5.5))

                # 캡션
                caption = doc.add_paragraph()
                caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
                cap_run = caption.add_run(f"Figure: {fig_name.replace('.png', '').replace('_', ' ')}")
                cap_run.font.size = Pt(9)
                cap_run.font.color.rgb = RGBColor(100, 100, 100)
                cap_run.italic = True

    # 저장
    output_path = os.path.join(REPORT_DIR, "EGFR_Virtual_Screening_Report.docx")
    doc.save(output_path)
    print(f"\nReport saved: {output_path}")
    return output_path


def main():
    print("=" * 60)
    print("Solar Pro4 보고서 생성 시작")
    print("=" * 60)

    # 1) Skill 파일 로드
    print("\n[1] Loading skill files...")
    skills = read_skill_files()
    for name, content in skills.items():
        print(f"  {name}: {len(content)} chars")

    # 2) Solar Pro4 호출하여 각 섹션 생성
    print("\n[2] Generating report sections via Solar Pro4...")
    sections = generate_sections(skills)

    # 결과 확인
    for key, text in sections.items():
        if text:
            print(f"  {key}: {len(text)} chars OK")
        else:
            print(f"  {key}: FAILED")

    # 원본 텍스트 저장 (디버깅용)
    raw_path = os.path.join(REPORT_DIR, "raw_sections.json")
    with open(raw_path, "w", encoding="utf-8") as f:
        json.dump(sections, f, ensure_ascii=False, indent=2)
    print(f"\nRaw sections saved: {raw_path}")

    # 3) DOCX 빌드
    print("\n[3] Building DOCX report...")
    output = build_docx(sections)

    print("\n" + "=" * 60)
    print("보고서 생성 완료!")
    print(f"파일: {output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
