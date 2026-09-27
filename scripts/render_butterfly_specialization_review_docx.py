#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from PIL import Image


RUNNING_TITLE = "Dimensions of butterfly specialization"
TICK = chr(96)
CODE_FENCE = TICK * 3

FIGURES = {
    1: "Figure1_taxonomic_vs_geographic_specialization.png",
    2: "Figure2_anthropogenic_resource_expansion.png",
    3: "Figure3_host_contribution_architecture.png",
    4: "Figure4_within_family_specialization_hierarchy.png",
    5: "Figure5_independent_climate_filtering.png",
}

INLINE = re.compile(
    r"(\*\*[^*]+\*\*|\*[^*]+\*|"
    + re.escape(TICK)
    + r"[^"
    + re.escape(TICK)
    + r"]+"
    + re.escape(TICK)
    + r")"
)


def add_line_numbering(section) -> None:
    sect_pr = section._sectPr
    old = sect_pr.find(qn("w:lnNumType"))
    if old is not None:
        sect_pr.remove(old)
    node = OxmlElement("w:lnNumType")
    node.set(qn("w:countBy"), "1")
    node.set(qn("w:distance"), "360")
    node.set(qn("w:restart"), "continuous")
    sect_pr.append(node)


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    separate = OxmlElement("w:fldChar")
    separate.set(qn("w:fldCharType"), "separate")
    display = OxmlElement("w:t")
    display.text = "1"
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    for node in (begin, instr, separate, display, end):
        run._r.append(node)


def configure_document(doc: Document) -> None:
    section = doc.sections[0]
    section.top_margin = Inches(1.0)
    section.bottom_margin = Inches(1.0)
    section.left_margin = Inches(1.0)
    section.right_margin = Inches(1.0)
    add_line_numbering(section)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    normal.paragraph_format.line_spacing = 2.0
    normal.paragraph_format.space_after = Pt(0)

    for name in ("Title", "Heading 1", "Heading 2", "Heading 3"):
        style = doc.styles[name]
        style.font.name = "Times New Roman"
        style.font.color.rgb = RGBColor(0, 0, 0)

    doc.styles["Title"].font.size = Pt(16)
    doc.styles["Title"].font.bold = True
    doc.styles["Heading 1"].font.size = Pt(13)
    doc.styles["Heading 1"].font.bold = True
    doc.styles["Heading 2"].font.size = Pt(12)
    doc.styles["Heading 2"].font.bold = True
    doc.styles["Heading 3"].font.size = Pt(12)
    doc.styles["Heading 3"].font.bold = True

    header = section.header.paragraphs[0]
    header.text = RUNNING_TITLE
    header.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in header.runs:
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)

    add_page_number(section.footer.paragraphs[0])

    cp = doc.core_properties
    cp.author = ""
    cp.last_modified_by = ""
    cp.comments = ""
    cp.keywords = ""
    cp.subject = ""
    cp.category = ""


def add_inline(paragraph, text: str) -> None:
    pos = 0
    for match in INLINE.finditer(text):
        if match.start() > pos:
            paragraph.add_run(text[pos:match.start()])
        token = match.group(0)
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        elif token.startswith("*"):
            run = paragraph.add_run(token[1:-1])
            run.italic = True
        else:
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Courier New"
            run.font.size = Pt(10.5)
        pos = match.end()
    if pos < len(text):
        paragraph.add_run(text[pos:])
    for run in paragraph.runs:
        if run.font.name is None:
            run.font.name = "Times New Roman"
        if run.font.size is None:
            run.font.size = Pt(12)


def add_figure(doc: Document, image_path: Path) -> None:
    with Image.open(image_path) as img:
        w, h = img.size
    max_w = 6.25
    max_h = 7.25
    ratio = w / h
    width = min(max_w, max_h * ratio)
    para = doc.add_paragraph()
    para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = para.add_run()
    picture = run.add_picture(str(image_path), width=Inches(width))
    picture._inline.docPr.set("descr", image_path.stem.replace("_", " "))
    para.paragraph_format.line_spacing = 1.0
    para.paragraph_format.space_before = Pt(6)
    para.paragraph_format.space_after = Pt(6)


def add_body_paragraph(doc: Document, text: str, *, first_line: bool = True):
    p = doc.add_paragraph()
    if first_line:
        p.paragraph_format.first_line_indent = Inches(0.25)
    add_inline(p, text)
    return p


def parse_markdown(source: str, figures_dir: Path) -> Document:
    doc = Document()
    configure_document(doc)

    lines = source.splitlines()
    in_code = False
    code_lines: list[str] = []
    in_figure_legends = False
    in_references = False

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()

        if line.startswith(CODE_FENCE):
            if not in_code:
                in_code = True
                code_lines = []
            else:
                p = doc.add_paragraph()
                p.paragraph_format.left_indent = Inches(0.35)
                p.paragraph_format.right_indent = Inches(0.35)
                p.paragraph_format.line_spacing = 1.0
                run = p.add_run("\n".join(code_lines))
                run.font.name = "Courier New"
                run.font.size = Pt(9.5)
                in_code = False
            i += 1
            continue

        if in_code:
            code_lines.append(line)
            i += 1
            continue

        if not line.strip():
            i += 1
            continue

        if line.strip() == "---":
            i += 1
            continue

        if line.startswith("# "):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.first_line_indent = Inches(0)
            p.paragraph_format.space_after = Pt(12)
            run = p.add_run(line[2:].strip())
            run.bold = True
            run.font.name = "Times New Roman"
            run.font.size = Pt(16)
            run.font.color.rgb = RGBColor(0, 0, 0)
            i += 1
            continue

        if line.startswith("## "):
            heading = line[3:].strip()
            p = doc.add_paragraph(heading, style="Heading 1")
            p.paragraph_format.keep_with_next = True
            in_figure_legends = heading.lower().startswith("figure legends")
            in_references = heading.lower() == "references"
            i += 1
            continue

        if line.startswith("### "):
            p = doc.add_paragraph(line[4:].strip(), style="Heading 2")
            p.paragraph_format.keep_with_next = True
            i += 1
            continue

        if line.startswith("#### "):
            p = doc.add_paragraph(line[5:].strip(), style="Heading 3")
            p.paragraph_format.keep_with_next = True
            i += 1
            continue

        if line.startswith("> "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.35)
            p.paragraph_format.right_indent = Inches(0.35)
            p.paragraph_format.first_line_indent = Inches(0)
            run = p.add_run(line[2:].strip())
            run.italic = True
            run.font.name = "Times New Roman"
            run.font.size = Pt(12)
            i += 1
            continue

        fig_match = re.match(r"^\*\*Figure\s+(\d+)\.", line)
        if in_figure_legends and fig_match:
            num = int(fig_match.group(1))
            image = figures_dir / FIGURES[num]
            if not image.exists():
                raise FileNotFoundError(image)
            doc.add_page_break()
            add_figure(doc, image)
            p = doc.add_paragraph()
            p.paragraph_format.first_line_indent = Inches(0)
            add_inline(p, line)
            i += 1
            continue

        if line.startswith("- "):
            if in_references:
                p = doc.add_paragraph()
                p.paragraph_format.left_indent = Inches(0.3)
                p.paragraph_format.first_line_indent = Inches(-0.3)
                p.paragraph_format.line_spacing = 2.0
                add_inline(p, line[2:].strip())
            else:
                p = doc.add_paragraph(style="List Bullet")
                p.paragraph_format.line_spacing = 2.0
                add_inline(p, line[2:].strip())
            i += 1
            continue

        if re.match(r"^\d+\.\s+", line):
            text = re.sub(r"^\d+\.\s+", "", line)
            p = doc.add_paragraph(style="List Number")
            p.paragraph_format.line_spacing = 2.0
            add_inline(p, text)
            i += 1
            continue

        if line.startswith("**Running title:**"):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.first_line_indent = Inches(0)
            add_inline(p, line)
            i += 1
            continue

        para_lines = [line]
        j = i + 1
        while j < len(lines):
            nxt = lines[j].rstrip()
            if not nxt.strip():
                break
            if (
                nxt.startswith("#")
                or nxt.startswith("- ")
                or nxt.startswith("> ")
                or nxt.startswith(CODE_FENCE)
                or nxt.strip() == "---"
                or re.match(r"^\d+\.\s+", nxt)
                or (in_figure_legends and re.match(r"^\*\*Figure\s+\d+\.", nxt))
            ):
                break
            para_lines.append(nxt)
            j += 1

        text = " ".join(part.strip() for part in para_lines)
        p = add_body_paragraph(doc, text, first_line=True)
        if any(
            text.startswith(prefix)
            for prefix in (
                "**Aim:**",
                "**Location:**",
                "**Time period:**",
                "**Major taxa studied:**",
                "**Methods:**",
                "**Results:**",
                "**Main conclusions:**",
                "**Keywords:**",
            )
        ):
            p.paragraph_format.first_line_indent = Inches(0)
        if text.startswith("**Keywords:**"):
            p.paragraph_format.space_after = Pt(6)
        i = j if j > i + 1 else i + 1

    return doc


def scrub_core_properties(path: Path) -> None:
    tmp = path.with_suffix(".scrubbed.docx")
    ns = {
        "dc": "http://purl.org/dc/elements/1.1/",
        "cp": "http://schemas.openxmlformats.org/package/2006/metadata/core-properties",
    }
    with zipfile.ZipFile(path, "r") as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "docProps/core.xml":
                root = ET.fromstring(data)
                for xpath in ("dc:creator", "cp:lastModifiedBy"):
                    node = root.find(xpath, ns)
                    if node is not None:
                        node.text = ""
                data = ET.tostring(root, encoding="utf-8", xml_declaration=True)
            zout.writestr(item, data)
    tmp.replace(path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-markdown", type=Path, required=True)
    ap.add_argument("--figures-dir", type=Path, required=True)
    ap.add_argument("--output-docx", type=Path, required=True)
    args = ap.parse_args()

    text = args.input_markdown.read_text(encoding="utf-8")
    doc = parse_markdown(text, args.figures_dir)
    args.output_docx.parent.mkdir(parents=True, exist_ok=True)
    doc.save(args.output_docx)
    scrub_core_properties(args.output_docx)
    print(args.output_docx)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
