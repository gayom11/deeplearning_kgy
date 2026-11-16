"""수정된 raw_sections.json으로 DOCX 재빌드"""
import os
import json
from docx import Document
from docx.shared import Inches, Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(PROJECT_DIR, "experiment", "figures")
REPORT_DIR = os.path.join(PROJECT_DIR, "report")

with open(os.path.join(REPORT_DIR, "raw_sections.json"), "r", encoding="utf-8") as f:
    sections = json.load(f)


def parse_markdown_table(text):
    """마크다운 테이블 텍스트를 행/열 리스트로 파싱"""
    rows = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line.startswith("|") or "---" in line:
            continue
        cells = [c.strip() for c in line.split("|")[1:-1]]
        if cells:
            rows.append(cells)
    return rows


def add_table_to_doc(doc, table_data):
    """파싱된 테이블을 docx 테이블로 추가"""
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
                for paragraph in cell.paragraphs:
                    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    for run in paragraph.runs:
                        run.font.size = Pt(9)
                        run.font.name = "맑은 고딕"
                        if i == 0:
                            run.bold = True


def add_section_text(doc, text):
    """마크다운 텍스트를 docx 문단으로 변환"""
    lines = text.split("\n")
    table_buffer = []
    in_table = False

    for line in lines:
        stripped = line.strip()

        # 빈 줄
        if not stripped:
            if in_table and table_buffer:
                # 테이블 종료 → 삽입
                add_table_to_doc(doc, parse_markdown_table("\n".join(table_buffer)))
                table_buffer = []
                in_table = False
            continue

        # 테이블 행
        if stripped.startswith("|"):
            in_table = True
            table_buffer.append(stripped)
            continue

        # 테이블 중이었으면 먼저 플러시
        if in_table and table_buffer:
            add_table_to_doc(doc, parse_markdown_table("\n".join(table_buffer)))
            table_buffer = []
            in_table = False

        # 헤딩 (## → 볼드)
        if stripped.startswith("#"):
            clean = stripped.lstrip("#").strip()
            # 최상위 섹션 헤딩은 스킵 (이미 doc.add_heading으로 추가됨)
            if stripped.startswith("# ") or stripped.startswith("## "):
                continue
            p = doc.add_paragraph()
            run = p.add_run(clean)
            run.bold = True
            run.font.size = Pt(11)
            run.font.name = "맑은 고딕"
            continue

        # **볼드 텍스트**로 시작하면 (표 캡션 등)
        if stripped.startswith("**") and stripped.endswith("**"):
            p = doc.add_paragraph()
            run = p.add_run(stripped.strip("*"))
            run.bold = True
            run.font.size = Pt(10)
            run.font.name = "맑은 고딕"
            continue

        # 일반 텍스트
        p = doc.add_paragraph()
        # **bold** 처리
        parts = stripped.split("**")
        for idx, part in enumerate(parts):
            if not part:
                continue
            run = p.add_run(part)
            run.font.size = Pt(11)
            run.font.name = "맑은 고딕"
            if idx % 2 == 1:
                run.bold = True

    # 남은 테이블 플러시
    if table_buffer:
        add_table_to_doc(doc, parse_markdown_table("\n".join(table_buffer)))


def build():
    doc = Document()

    # 기본 스타일
    style = doc.styles["Normal"]
    font = style.font
    font.name = "맑은 고딕"
    font.size = Pt(11)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), "맑은 고딕")

    for section in doc.sections:
        section.top_margin = Cm(2.5)
        section.bottom_margin = Cm(2.5)
        section.left_margin = Cm(2.5)
        section.right_margin = Cm(2.5)

    # 표지
    for _ in range(6):
        doc.add_paragraph()

    title = doc.add_heading("딥러닝 기반 EGFR 타겟 비소세포폐암\n신약 후보 물질 Virtual Screening", level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in title.runs:
        run.font.size = Pt(20)
        run.font.name = "맑은 고딕"

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run("\nDeep Learning Final Project\n")
    run.font.size = Pt(13)
    run.font.color.rgb = RGBColor(80, 80, 80)

    doc.add_page_break()

    # 섹션 구성
    section_config = [
        ("1. 서론", "intro", []),
        ("2. 타겟 단백질 정보", "target", []),
        ("3. Virtual Screening 파이프라인", "pipeline", ["fig08_pipeline_funnel.png"]),
        ("4. 모델 선정 및 설명", "models", []),
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
        # 헤딩 레벨 결정
        if heading_text.startswith("5."):
            level = 2
        else:
            level = 1

        h = doc.add_heading(heading_text, level=level)
        for run in h.runs:
            run.font.size = Pt(14) if level == 1 else Pt(12)

        # 본문
        if key and key in sections and sections[key]:
            add_section_text(doc, sections[key])

        # 그림 삽입
        for fig_name in figures:
            fig_path = os.path.join(FIGURES_DIR, fig_name)
            if os.path.exists(fig_path):
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                run = p.add_run()
                # 넓은 그림 vs 일반 그림
                width = Inches(6.0) if "heatmap" in fig_name or "structure" in fig_name or "table" in fig_name else Inches(5.5)
                run.add_picture(fig_path, width=width)

                # 캡션
                caption_text = fig_name.replace(".png", "").replace("_", " ").replace("fig", "Figure ")
                cap = doc.add_paragraph()
                cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
                cap_run = cap.add_run(caption_text)
                cap_run.font.size = Pt(9)
                cap_run.font.color.rgb = RGBColor(100, 100, 100)
                cap_run.italic = True

    output_path = os.path.join(REPORT_DIR, "EGFR_Virtual_Screening_Report.docx")
    doc.save(output_path)
    print(f"Report saved: {output_path}")
    return output_path


if __name__ == "__main__":
    build()
