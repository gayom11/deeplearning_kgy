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
from docx.oxml.ns import qn

API_KEY = os.environ.get("SOLAR_API_KEY", "")
if not API_KEY:
    print("ERROR: SOLAR_API_KEY 환경변수를 설정하세요. (예: set SOLAR_API_KEY=your_key)")
    sys.exit(1)

API_URL = "https://api.upstage.ai/v1/chat/completions"
MODEL = "solar-pro4"

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SKILL_DIR = os.path.join(PROJECT_DIR, "skill")
FIGURES_DIR = os.path.join(PROJECT_DIR, "experiment", "figures")
REPORT_DIR = os.path.join(PROJECT_DIR, "report")
os.makedirs(REPORT_DIR, exist_ok=True)


def read_skills():
    content = {}
    for f in sorted(os.listdir(SKILL_DIR)):
        if f.endswith(".md"):
            with open(os.path.join(SKILL_DIR, f), "r", encoding="utf-8") as fh:
                content[f] = fh.read()
    return content


def call_solar(system_prompt, user_prompt, max_tokens=4000):
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
            return resp.json()["choices"][0]["message"]["content"].strip()
        except Exception as e:
            print(f"  API error (attempt {attempt+1}): {e}")
            time.sleep(5)
    return None


SYSTEM = """당신은 약학/신약개발 분야의 학술 보고서를 작성하는 전문가입니다.
규칙:
- 한국어로 작성
- 학술 보고서 문체 ("~하였다", "~이다")
- 간결하게 서술, 한 문단 최대 3줄
- 불필요한 수식어, 반복 설명 금지
- 표와 그림 참조를 적극 활용
- AI 특유의 과잉 설명 하지 말 것"""


def generate_sections(skills):
    full_skill = "\n\n---\n\n".join(f"## {k}\n{v}" for k, v in skills.items())
    sections = {}

    prompts = {
        "intro": ("1. 서론", "NSCLC 개요, EGFR 타겟 근거, 기존 EGFR 억제제 표(세대/약물/IC50), 연구 목적. 분량: 8줄 + 표 1개.", 1200),
        "target": ("2. 타겟 단백질 정보", "EGFR 단백질 기능, ChEMBL ID, PDB 1M17 선정 근거. 분량: 6줄 + 표 1개.", 1000),
        "pipeline": ("3. Virtual Screening 파이프라인", "fig10_pipeline_funnel.png 참조. 학습 데이터와 스크리닝 라이브러리 분리 설명. 각 단계 한 줄 설명: ChEMBL 50K → Lipinski 47K → ChemBERTa 1K → ADMET 300 → GNINA → 최종 10. 분량: 12줄.", 1200),
        "models": ("4. 모델 선정", "6종 모델 표(모델명/유형/입력/출력/선정 근거). CPU+Docker 환경 제약 1줄. 분량: 표 1개 + 5줄.", 1500),
        "results_preprocess": ("5.1 전처리", "Lipinski+PAINS 필터. fig01 참조. 50,000 → 47,050 (5.9% 제거). 분량: 4줄.", 800),
        "results_chemberta": ("5.2 ChemBERTa 모델 학습", "fig02 참조. Active 3,704 / Inactive 1,298. ROC-AUC 0.897. OOF 예측 사용. 분량: 5줄.", 800),
        "results_screening": ("5.3 Virtual Screening 결과", "fig03 참조. 47,050개 → 상위 1,000개 (확률 ≥ 0.802). 분량: 4줄.", 800),
        "results_docking": ("5.4 GNINA 도킹", "fig07 참조. PDB 1M17 수용체, Erlotinib autobox, 300개 도킹. CNN affinity + Vina score. 분량: 5줄.", 800),
        "results_admet": ("5.5 ADMET 프로파일", "fig04/fig05 참조. 16속성 평가. HIA 평균 0.958, AMES 0.333. 물리화학적 성질 포함. 분량: 5줄.", 800),
        "results_generation": ("5.6 분자 생성", "fig11 참조. SELFIES 돌연변이 + 스캐폴드 호핑. 778개 drug-like, 73% novel. 분량: 5줄.", 800),
        "candidates": ("6. 최종 후보 물질", "fig08/fig09 참조. 종합 점수 가중치(CNN 30%, ChemBERTa 25%, ADMET 25%, Vina 20%). Tanimoto > 0.7 중복 제거. 최종 10개 표(ID, CNN affinity, Vina, ADMET, 최종점수). 분량: 6줄 + 표 1개.", 1500),
        "conclusion": ("7. 결론", "6개 모델 활용, 학습/스크리닝 분리 의의, CPU+Docker로 GNINA 실현, 한계(실험 검증 부재), 향후 과제. 분량: 6줄.", 800),
    }

    for key, (title, instruction, tokens) in prompts.items():
        print(f"Generating: {title}...")
        sections[key] = call_solar(SYSTEM,
            f'아래 skill 자료를 참고하여 보고서의 "{title}" 섹션을 작성하세요.\n\n내용:\n{instruction}\n\n{full_skill}',
            max_tokens=tokens)

    return sections


def parse_md_table(text):
    rows = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line.startswith("|") or "---" in line:
            continue
        cells = [c.strip() for c in line.split("|")[1:-1]]
        if cells:
            rows.append(cells)
    return rows


def add_table(doc, table_data):
    if len(table_data) < 2:
        return
    n_cols = len(table_data[0])
    table = doc.add_table(rows=len(table_data), cols=n_cols)
    table.style = "Light Grid Accent 1"
    table.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for i, row_data in enumerate(table_data):
        for j, cell_text in enumerate(row_data):
            if j < n_cols:
                cell = table.cell(i, j)
                cell.text = cell_text
                for p in cell.paragraphs:
                    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    for run in p.runs:
                        run.font.size = Pt(9)
                        run.font.name = "맑은 고딕"
                        if i == 0:
                            run.bold = True


def add_section_text(doc, text):
    if not text:
        return
    table_buf = []
    in_table = False
    for line in text.split("\n"):
        s = line.strip()
        if not s:
            if in_table and table_buf:
                add_table(doc, parse_md_table("\n".join(table_buf)))
                table_buf = []; in_table = False
            continue
        if s.startswith("|"):
            in_table = True; table_buf.append(s); continue
        if in_table and table_buf:
            add_table(doc, parse_md_table("\n".join(table_buf)))
            table_buf = []; in_table = False
        if s.startswith("#"):
            continue
        p = doc.add_paragraph()
        parts = s.split("**")
        for idx, part in enumerate(parts):
            if not part:
                continue
            run = p.add_run(part)
            run.font.size = Pt(11)
            run.font.name = "맑은 고딕"
            if idx % 2 == 1:
                run.bold = True
    if table_buf:
        add_table(doc, parse_md_table("\n".join(table_buf)))


def build_docx(sections):
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "맑은 고딕"
    style.font.size = Pt(11)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), "맑은 고딕")
    for sec in doc.sections:
        sec.top_margin = Cm(2.5); sec.bottom_margin = Cm(2.5)
        sec.left_margin = Cm(2.5); sec.right_margin = Cm(2.5)

    # 표지
    for _ in range(6):
        doc.add_paragraph()
    t = doc.add_heading("딥러닝 기반 EGFR 타겟 비소세포폐암\n신약 후보 물질 Virtual Screening", level=0)
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for r in t.runs:
        r.font.size = Pt(20)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("\nDeep Learning Final Project\n")
    r.font.size = Pt(13)
    r.font.color.rgb = RGBColor(80, 80, 80)
    doc.add_page_break()

    config = [
        ("1. 서론", "intro", []),
        ("2. 타겟 단백질 정보", "target", []),
        ("3. Virtual Screening 파이프라인", "pipeline", ["fig10_pipeline_funnel.png"]),
        ("4. 모델 선정 및 설명", "models", []),
        ("5.1 전처리", "results_preprocess", ["fig01_preprocessing.png"]),
        ("5.2 ChemBERTa 모델 학습", "results_chemberta", ["fig02_chemberta_training.png"]),
        ("5.3 Virtual Screening 결과", "results_screening", ["fig03_screening_results.png"]),
        ("5.4 GNINA 도킹", "results_docking", ["fig07_docking_results.png"]),
        ("5.5 ADMET 프로파일", "results_admet", ["fig04_admet_heatmap.png", "fig05_admet_radar.png"]),
        ("5.6 분자 생성", "results_generation", ["fig11_generated_molecules.png"]),
        ("6. 최종 후보 물질", "candidates", ["fig08_candidate_structures.png", "fig09_scatter_analysis.png"]),
        ("7. 결론", "conclusion", []),
    ]

    for heading, key, figs in config:
        level = 2 if heading.startswith("5.") else 1
        h = doc.add_heading(heading, level=level)
        for r in h.runs:
            r.font.size = Pt(14) if level == 1 else Pt(12)

        if key in sections and sections[key]:
            add_section_text(doc, sections[key])

        for fig in figs:
            fp = os.path.join(FIGURES_DIR, fig)
            if os.path.exists(fp):
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.add_run().add_picture(fp, width=Inches(5.5))
                cap = doc.add_paragraph()
                cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
                r = cap.add_run(fig.replace(".png","").replace("_"," ").replace("fig","Figure "))
                r.font.size = Pt(9); r.font.color.rgb = RGBColor(100,100,100); r.italic = True

    out = os.path.join(REPORT_DIR, "EGFR_Virtual_Screening_Report.docx")
    doc.save(out)
    print(f"\nReport saved: {out}")
    return out


def main():
    print("="*50)
    print("Solar Pro4 보고서 생성")
    print("="*50)

    skills = read_skills()
    for name, content in skills.items():
        print(f"  {name}: {len(content)} chars")

    sections = generate_sections(skills)

    raw = os.path.join(REPORT_DIR, "raw_sections.json")
    with open(raw, "w", encoding="utf-8") as f:
        json.dump(sections, f, ensure_ascii=False, indent=2)

    for k, v in sections.items():
        print(f"  {k}: {'OK' if v else 'FAILED'} ({len(v) if v else 0} chars)")

    build_docx(sections)
    print("완료!")


if __name__ == "__main__":
    main()
