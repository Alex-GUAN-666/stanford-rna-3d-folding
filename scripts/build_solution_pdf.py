"""Build the public case study and its text companion. Requires reportlab.

Run from any directory: python scripts/build_solution_pdf.py
This documentation build does not run models or read competition data.
"""
from pathlib import Path
import html
import re
import reportlab

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak,
    Flowable,
)

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://github.com/Alex-GUAN-666/stanford-rna-3d-folding"
INK = colors.HexColor("#172D40")
TEAL = colors.HexColor("#007F83")
MUTED = colors.HexColor("#526777")
PALE = colors.HexColor("#EFF5F7")
WIDTH = A4[0] - 100
FONT_DIR = Path(reportlab.__file__).resolve().parent / "fonts"
pdfmetrics.registerFont(TTFont("CaseSans", str(FONT_DIR / "Vera.ttf")))
pdfmetrics.registerFont(TTFont("CaseSans-Bold", str(FONT_DIR / "VeraBd.ttf")))
pdfmetrics.registerFont(TTFont("CaseSans-Italic", str(FONT_DIR / "VeraIt.ttf")))
pdfmetrics.registerFontFamily("CaseSans", normal="CaseSans", bold="CaseSans-Bold", italic="CaseSans-Italic", boldItalic="CaseSans-Bold")
STYLES = {
    "title": ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=29, leading=33, textColor=INK, spaceAfter=12),
    "subtitle": ParagraphStyle("subtitle", fontName="Helvetica", fontSize=13, leading=18, textColor=TEAL, spaceAfter=13),
    "h2": ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=12, leading=16, textColor=INK, spaceBefore=11, spaceAfter=6),
    "body": ParagraphStyle("body", fontName="Helvetica", fontSize=10.2, leading=13.5, textColor=INK, spaceAfter=7),
    "small": ParagraphStyle("small", fontName="Helvetica", fontSize=8.5, leading=11.6, textColor=MUTED, spaceAfter=6),
    "cell": ParagraphStyle("cell", fontName="Helvetica", fontSize=9, leading=12, textColor=INK),
    "head": ParagraphStyle("head", fontName="Helvetica-Bold", fontSize=9, leading=12, textColor=colors.white),
    "code": ParagraphStyle("code", fontName="Courier", fontSize=8.4, leading=12, textColor=INK, backColor=PALE, borderPadding=9, spaceAfter=10),
}
for style in STYLES.values():
    if style.fontName == "Helvetica":
        style.fontName = "CaseSans"
    elif style.fontName == "Helvetica-Bold":
        style.fontName = "CaseSans-Bold"
STYLES["title"].fontSize = 26
STYLES["title"].leading = 31
STYLES["body"].fontSize = 9.7
STYLES["cell"].fontSize = 8.5
STYLES["small"].fontSize = 8


def p(text, style="body"):
    return Paragraph(text, STYLES[style])


class Workflow(Flowable):
    """Compare the two separate historical notebook workflows."""
    def __init__(self):
        super().__init__()
        self.width = WIDTH
        self.height = 138

    def draw(self):
        c = self.canv
        gap = 17
        w = (WIDTH - gap) / 2
        for x, title, lines in [
            (0, "Protenix + trRNA backup", ["Five Protenix candidates", "PyRosetta ordering", "For <=300 nt: replace slot 5", "with trRNA candidate 1"]),
            (w + gap, "DRfold2 backup", ["Length-dependent checkpoints", "480 nt windows for long inputs", "Order by distance to mean score", "Keep five candidates"]),
        ]:
            c.setFillColor(PALE)
            c.roundRect(x, 6, w, 127, 8, fill=1, stroke=0)
            c.setFillColor(TEAL)
            c.roundRect(x, 102, w, 31, 8, fill=1, stroke=0)
            c.rect(x, 102, w, 10, fill=1, stroke=0)
            c.setFillColor(colors.white)
            c.setFont("CaseSans-Bold", 9)
            c.drawString(x + 12, 114, title)
            c.setFont("CaseSans", 8.6)
            c.setFillColor(INK)
            for i, line in enumerate(lines):
                c.drawString(x + 12, 85 - i * 18, line)


def table(headers, rows, widths):
    data = [[p(html.escape(x), "head") for x in headers]]
    data += [[p(x, "cell") for x in row] for row in rows]
    t = Table(data, colWidths=[WIDTH * x for x in widths], hAlign="LEFT", repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), TEAL),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [PALE, colors.white]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 9),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
        ("LINEBELOW", (0, 0), (-1, 0), 0.5, TEAL),
    ]))
    return t


PAGES = [{'title': 'Stanford RNA\n3D Folding',
  'subtitle': 'Model integration and RNA candidate processing',
  'blocks': [('small', 'Yuzhen Guan (Alex) | Team GZSL | September 2026'),
             ('table',
              ['COMPETITION RESULT', 'TASK', 'REPOSITORY'],
              [['<b>0.42295 private</b><br/>76 / 1,516 teams<br/>Bronze medal [S1,S2]',
                "Five candidate structures per RNA; one C1' coordinate per nucleotide",
                'Competition notebooks, technical notes, and CPU utilities']],
              [0.31, 0.35, 0.34]),
             ('h2', 'About the project'),
             ('body',
              'I competed in Stanford RNA 3D Folding as part of Team GZSL. My work focused on coordinating '
              'experiments, analyzing results by sequence length, and integrating predictions from pretrained '
              'models. We finished 76th out of 1,516 teams with a final private leaderboard score of 0.42295.'),
             ('body',
              'Our approach combined Protenix, DRfold2 and a pipeline recorded in my notes as trRNA. The practical '
              'questions were how to allocate inference work across sequence lengths, how to choose five candidate '
              'structures, and how to keep every predicted coordinate matched to the right nucleotide.'),
             ('h2', 'Results'),
             ('table',
              ['Evaluation', 'Score', 'Context'],
              [['Final private leaderboard', '<b>0.42295</b>', 'Team GZSL, rank 76; bronze.'],
               ['Later public test', '0.409', 'Three-model ensemble result recorded in my solution notes.'],
               ['Development comparison',
                '0.334 to 0.368',
                'Retained report; the exact validation split is no longer clear.']],
              [0.3, 0.23, 0.47]),
             ('small',
              'I keep these results separate because they come from different evaluation contexts. The '
              'public/private difference is not a same-split improvement. Detailed experiment records are in '
              'docs/RESULTS.md.'),
             ('h2', 'Why I maintain this repository'),
             ('body',
              'I brought together the inference backups I still have and a small package for working with '
              'predicted coordinates. I added the CPU utilities after the competition so that the data checks, '
              'candidate handling and segment alignment can be run independently of a model installation.'),
             ('body',
              'The original folding models remain upstream dependencies. This repository runs the '
              'coordinate-processing workflow and synthetic examples; recovering the exact scored GPU run still '
              'requires the matching model assets and final submission CSV.')]},
 {'title': 'Model integration',
  'subtitle': 'What is in the retained inference notebooks',
  'blocks': [('body',
              'My solution notes describe length-dependent preferences for trRNA, DRfold2 and Protenix. The two '
              'retained notebooks capture separate experiments from that work: Protenix with trRNA candidate replacement, '
              'and DRfold2 inference. I keep their actual settings visible instead of treating them as one final '
              'submission script.'),
             ('workflow',),
             ('h2', 'Protenix with trRNA candidate replacement'),
             ('body',
              'The notebook uses Protenix 0.4.6, model_v0.2.0 weights, 10 cycles, 200 diffusion steps and seed 42, '
              'with MSA disabled. It sorts five Protenix candidates by PyRosetta score. For targets up to 300 nt, '
              'it retains the first four and replaces slot 5 with the first trRNA candidate. The structures stay '
              'separate; this step does not average their coordinates.'),
             ('h2', 'DRfold2 candidate generation'),
             ('body',
              'The DRfold2 backup uses fp32 and cfg_97. Its checkpoint counts are 20, 10, 5 and 5 for lengths '
              '<=100, 101-200, 201-480 and >480 nt. With GET_CENTER=True, candidates are ordered by distance from '
              'the mean geometric heuristic score, rather than minimum energy.'),
             ('h2', 'Long sequences'),
             ('body',
              'The Protenix backup predicts up to 960 nt and pads the remaining coordinates with zeros. DRfold2 '
              'uses 480 nt windows with a one-residue overlap and truncates beyond 2400 nt. These shortcuts leave '
              'incomplete structural predictions. In the current utilities I require full residue coverage and '
              'enough overlap to determine a rigid alignment.'),
             ('h2', 'Status of the backups'),
             ('body',
              'The Protenix/trRNA notebook retained outputs for 12 Protenix targets and nine trRNA replacements '
              'before a final concatenation error. The DRfold2 notebook has no saved outputs. I documented '
              'dependencies and known issues in docs/HISTORICAL_NOTES.md. The two Python backups use a different '
              'Protenix configuration: 400 diffusion steps and seed 0.'),
             ('small',
              'The retained files contain inference code. Training records for the fine-tuning mentioned in my '
              'earlier notes are missing, as are the exact historical trRNA source snapshot and the final scored '
              'submission.')]},
 {'title': 'Coordinate processing',
  'subtitle': 'The CPU package and its design choices',
  'blocks': [('body',
              'I use the current Python package to inspect, assemble and export coordinates produced by external '
              'predictors. Separating this work from neural inference makes the data contracts easier to check and '
              'keeps the basic workflow runnable on a laptop.'),
             ('table',
              ['Component', 'Implementation', 'Purpose'],
              [['Input validation',
                'Check target IDs, residue order, finite values and array shapes.',
                'Keep each coordinate tied to the correct nucleotide.'],
               ['Window assembly',
                'Kabsch alignment on corresponding overlap residues; equal blending after alignment.',
                'Bring independently predicted windows into a common frame.'],
               ['Candidate selection',
                'Round-robin across model pools; energy ordering within declared comparable groups.',
                'Retain multiple sources without mixing incompatible energy scales.'],
               ['Submission export',
                'Five distinct candidate IDs and the 18-column CSV schema.',
                'Catch incomplete candidate sets and malformed submissions.']],
              [0.22, 0.42, 0.36]),
             ('h2', 'Routing and window boundaries'),
             ('body',
              'The current planner uses the model labels trrosettarna for 1-100 nt, drfold2 for 101-480 nt, and '
              'protenix above 480 nt. These defaults follow the written strategy; the notebooks use their own '
              'policies. The trrosettarna label is a working name for the trRNA route, whose exact upstream '
              'snapshot is unresolved.'),
             ('body',
              'Default windows are 300 nt long with 50 nt overlap. The first window anchors the coordinate frame. '
              'Each additional window needs at least three non-collinear corresponding points; gaps and ambiguous '
              'fits are rejected. Averaging the aligned overlaps does not resolve long-range contacts between '
              'distant windows.'),
             ('h2', 'Selection and evaluation'),
             ('body',
              'Different candidate IDs or model sources do not guarantee different folds. Round-robin selection is '
              'a transparent baseline, and energy values are only compared within explicitly compatible groups. '
              'Structural diversity and best-of-five accuracy need to be measured on real targets.'),
             ('body',
              'For evaluation, the organizer reference compares candidate structures with native conformers and '
              'aggregates the best comparisons by target. The optional local USalign helper evaluates one PDB '
              'pair. It does not implement the complete competition scoring procedure. [S3]')]},
 {'title': 'Running the project',
  'subtitle': 'Reproducibility, limitations and next experiments',
  'blocks': [('body',
              'After cloning the repository, I recommend a fresh Python environment and a normal package '
              'installation. The default workflow needs CPU and NumPy; it does not download model weights or '
              'competition data.'),
             ('code',
              'python -m pip install .<br/>python scripts/verify_repo.py --output outputs/verification.json'),
             ('body',
              'The verifier runs 45 tests and 11 command checks, including the installed package outside the '
              'checkout, input export, submission import and numeric roundtrip checks. The synthetic demo contains '
              'three targets and produces 752 rows, 18 columns and five candidates per target. Source hashes check '
              'that the installed package matches the checkout.'),
             ('body',
              'The GitHub workflow checks Windows and Ubuntu with Python 3.10 and 3.12. Its Ubuntu / Python 3.12 '
              'job also executes the seven-cell walkthrough in a Jupyter kernel. I link each run and its '
              'downloadable report so readers can inspect the results for a specific commit.'),
             ('small',
              '<link href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/docs/VERIFY.md" '
              'color="#007F83">Installation and verification guide</link> | <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions/workflows/tests.yml" '
              'color="#007F83">GitHub Actions</link> | <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/docs/REPRODUCIBILITY.md" '
              'color="#007F83">Input formats and commands</link>'),
             ('h2', 'What I would evaluate next'),
             ('body',
              'With matching model assets and a fixed allowed-data split, I would compare a common-model baseline '
              'with length routing, then compare selection policies on the same saved candidate pool. For long '
              'sequences I would report both fold quality and boundary geometry, alongside runtime, GPU memory and '
              'complete-target coverage. These are proposed experiments, not completed benchmarks.'),
             ('h2', 'References and project records'),
             ('small',
              '[S1] <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/evidence/certificate.jpg" '
              'color="#007F83">Kaggle award certificate</link>: Alex GUAN, 76 / 1,516 teams, bronze, September 24, '
              '2025.'),
             ('small',
              '[S2] <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/evidence/private_leaderboard.png" '
              'color="#007F83">Final leaderboard screenshot</link>: Team GZSL, 0.42295. My <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/docs/RESULTS.md" '
              'color="#007F83">results notes</link> separate this private result from the later public score of '
              '0.409.'),
             ('small',
              '[S3] <link href="https://github.com/DasLab/k1_tools" color="#007F83">DasLab competition '
              'tools</link> and <link href="https://github.com/DasLab/k1_tools/blob/main/ribonanza_tmscore.py" '
              'color="#007F83">TM-score reference</link>; <link href="https://github.com/pylelab/USalign" '
              'color="#007F83">USalign</link>.'),
             ('small',
              '[S4] Upstream models: <link href="https://github.com/bytedance/Protenix" '
              'color="#007F83">Protenix</link>, <link href="https://github.com/leeyang/DRfold2" '
              'color="#007F83">DRfold2</link>. See <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/docs/REFERENCES.md" '
              'color="#007F83">references</link> for trRNA-related components and attribution.'),
             ('small',
              '[S5] <link href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/tree/main/historical" '
              'color="#007F83">Retained inference notebooks</link>; <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/docs/HISTORICAL_NOTES.md" '
              'color="#007F83">dependency and implementation notes</link>; <link '
              'href="https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/blob/main/evidence/source_manifest.json" '
              'color="#007F83">source file inventory</link>.'),
             ('small',
              'This September 2026 write-up brings together my competition notes and the current repository. The '
              'original upstream model authors retain credit for their architectures and implementations.')]}]


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(TEAL)
    canvas.setLineWidth(2)
    canvas.line(50, A4[1] - 37, A4[0] - 50, A4[1] - 37)
    canvas.setFillColor(MUTED)
    canvas.setFont("CaseSans", 7)
    canvas.drawString(50, 29, "ALEX GUAN  /  RNA 3D FOLDING  /  TECHNICAL CASE STUDY")
    canvas.drawRightString(A4[0] - 50, 29, f"{doc.page} / 4")
    canvas.restoreState()


def markdown(text):
    text = re.sub(r'<link href="([^"]+)"[^>]*>(.*?)</link>', r'[\2](\1)', text)
    text = text.replace("<br/>", " ").replace("<b>", "**").replace("</b>", "**")
    return html.unescape(text)


def build():
    story = []
    md = ["# Stanford RNA 3D Folding: technical case study", "", "September 2026. This is the text companion to [the PDF](../evidence/solution_notes.pdf).", "", "To rebuild both: install `reportlab`, then run `python scripts/build_solution_pdf.py` from the repository root. Edit the `PAGES` content in that script.", ""]
    for i, page in enumerate(PAGES):
        if i:
            story.append(PageBreak())
        story += [p(page["title"].replace("\n", "<br/>"), "title"), p(page["subtitle"], "subtitle")]
        md += ["## " + page["title"].replace("\n", " "), "", page["subtitle"], ""]
        for block in page["blocks"]:
            kind = block[0]
            if kind == "table":
                story += [table(*block[1:]), Spacer(1, 5)]
                md += ["| " + " | ".join(block[1]) + " |", "| " + " | ".join(["---"] * len(block[1])) + " |"]
                md += ["| " + " | ".join(markdown(v) for v in row) + " |" for row in block[2]]
            elif kind == "workflow":
                story.append(Workflow())
                md += ["| Protenix + trRNA backup | DRfold2 backup |", "| --- | --- |", "| Five Protenix candidates; PyRosetta ordering; for <=300 nt replace slot 5 with trRNA candidate 1. | Length-dependent checkpoints; 480 nt windows for long inputs; order by distance to mean score; keep five. |"]
            else:
                story.append(p(block[1], kind))
                if kind == "h2":
                    md.append("### " + markdown(block[1]))
                elif kind == "code":
                    md += ["```bash", block[1].replace("<br/>", "\n"), "```"]
                else:
                    md.append(markdown(block[1]))
            md.append("")
    out = ROOT / "evidence/solution_notes.pdf"
    doc = SimpleDocTemplate(str(out), pagesize=A4, leftMargin=50, rightMargin=50,
                            topMargin=56, bottomMargin=51, title="Stanford RNA 3D Folding - Technical Case Study",
                            author="Guan Yuzhen (Alex Guan)", subject="RNA model integration and coordinate processing; September 2026",
                            creator="ReportLab; technical documentation", pageCompression=1)
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    (ROOT / "docs/SOLUTION.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    build()
