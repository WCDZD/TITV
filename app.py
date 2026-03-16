from __future__ import annotations

import csv
import html
import io
import sqlite3
from datetime import datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs
from uuid import uuid4

BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "data" / "tv_mock_data.csv"
DB_PATH = BASE_DIR / "ti_tv.db"
REPORT_DIR = BASE_DIR / "reports"
CSS_FILE = BASE_DIR / "static" / "style.css"


def init_storage() -> None:
    REPORT_DIR.mkdir(exist_ok=True)
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS ti_requests (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                disease_background TEXT,
                target_mechanism TEXT,
                data_requirements TEXT,
                priority TEXT,
                output_format TEXT,
                status TEXT NOT NULL DEFAULT 'pending',
                created_at TEXT NOT NULL
            )
            """
        )


def load_tv_data() -> list[dict[str, str]]:
    with DATA_FILE.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def disease_categories(rows: list[dict[str, str]]) -> list[str]:
    return sorted({row["disease_category"] for row in rows})


def summarize_tv(rows: list[dict[str, str]], targets: list[str], disease: str) -> list[dict[str, float | int | str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        if row["disease_category"] == disease and row["gene"].upper() in targets:
            grouped.setdefault(row["gene"].upper(), []).append(row)

    summary: list[dict[str, float | int | str]] = []
    for gene, items in sorted(grouped.items()):
        expr = sum(float(i["expression"]) for i in items) / len(items)
        mut = 100 * (sum(int(i["mutation_status"]) for i in items) / len(items))
        summary.append({"gene": gene, "mean_expression": expr, "mutation_rate": mut, "sample_count": len(items)})
    return summary


class TinyPDF:
    def __init__(self) -> None:
        self.stream = io.StringIO()

    @staticmethod
    def _esc(text: str) -> str:
        return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")

    def text(self, x: float, y: float, size: int, content: str) -> None:
        self.stream.write(f"BT /F1 {size} Tf {x:.1f} {y:.1f} Td ({self._esc(content)}) Tj ET\n")

    def rect(self, x: float, y: float, w: float, h: float) -> None:
        self.stream.write(f"{x:.1f} {y:.1f} {w:.1f} {h:.1f} re S\n")

    def filled_rect(self, x: float, y: float, w: float, h: float) -> None:
        self.stream.write(f"{x:.1f} {y:.1f} {w:.1f} {h:.1f} re f\n")

    def rgb(self, r: float, g: float, b: float) -> None:
        self.stream.write(f"{r:.2f} {g:.2f} {b:.2f} rg\n")

    def line(self, x1: float, y1: float, x2: float, y2: float) -> None:
        self.stream.write(f"{x1:.1f} {y1:.1f} m {x2:.1f} {y2:.1f} l S\n")

    def to_pdf(self) -> bytes:
        content = self.stream.getvalue().encode("latin-1", errors="replace")
        objects = []
        objects.append(b"1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\n")
        objects.append(b"2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\n")
        objects.append(b"3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >> endobj\n")
        objects.append(b"4 0 obj << /Type /Font /Subtype /Type1 /BaseFont /Helvetica >> endobj\n")
        objects.append(f"5 0 obj << /Length {len(content)} >> stream\n".encode("latin-1") + content + b"endstream endobj\n")

        out = io.BytesIO()
        out.write(b"%PDF-1.4\n")
        offsets = [0]
        for obj in objects:
            offsets.append(out.tell())
            out.write(obj)
        xref_pos = out.tell()
        out.write(f"xref\n0 {len(objects)+1}\n".encode())
        out.write(b"0000000000 65535 f \n")
        for off in offsets[1:]:
            out.write(f"{off:010d} 00000 n \n".encode())
        out.write(f"trailer << /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF".encode())
        return out.getvalue()


def build_pdf(summary: list[dict[str, float | int | str]], disease: str, targets: list[str]) -> bytes:
    pdf = TinyPDF()
    y = 810
    pdf.text(40, y, 16, "Target Validation (TV) Report")
    y -= 22
    pdf.text(40, y, 11, f"Disease Category: {disease}")
    y -= 16
    pdf.text(40, y, 11, f"Targets: {', '.join(targets)}")
    y -= 16
    pdf.text(40, y, 10, f"Generated(UTC): {datetime.utcnow().isoformat(timespec='seconds')}")

    y -= 26
    pdf.text(40, y, 12, "Summary Table")
    y -= 16
    pdf.text(40, y, 10, "Gene")
    pdf.text(120, y, 10, "Mean Expression")
    pdf.text(260, y, 10, "Mutation Rate(%)")
    pdf.text(410, y, 10, "N")
    y -= 8
    pdf.line(40, y, 560, y)
    for row in summary:
        y -= 18
        pdf.text(40, y, 10, str(row["gene"]))
        pdf.text(120, y, 10, f"{row['mean_expression']:.2f}")
        pdf.text(260, y, 10, f"{row['mutation_rate']:.1f}")
        pdf.text(410, y, 10, str(row["sample_count"]))

    base_y = 220
    chart_h = 140
    bar_w = 60
    spacing = 90

    # Expression chart
    pdf.text(40, 390, 11, "Expression by Target")
    pdf.line(40, base_y, 260, base_y)
    pdf.line(40, base_y, 40, base_y + chart_h)
    max_expr = max(float(r["mean_expression"]) for r in summary) if summary else 1
    for idx, row in enumerate(summary):
        h = (float(row["mean_expression"]) / max_expr) * (chart_h - 10)
        x = 55 + idx * spacing
        pdf.rgb(0.15, 0.38, 0.72)
        pdf.filled_rect(x, base_y, bar_w, h)
        pdf.rgb(0, 0, 0)
        pdf.text(x, base_y - 14, 9, str(row["gene"]))

    # Mutation chart
    pdf.text(320, 390, 11, "Mutation Rate by Target")
    pdf.line(320, base_y, 560, base_y)
    pdf.line(320, base_y, 320, base_y + chart_h)
    for idx, row in enumerate(summary):
        h = (float(row["mutation_rate"]) / 100.0) * (chart_h - 10)
        x = 335 + idx * spacing
        pdf.rgb(0.76, 0.19, 0.17)
        pdf.filled_rect(x, base_y, bar_w, h)
        pdf.rgb(0, 0, 0)
        pdf.text(x, base_y - 14, 9, str(row["gene"]))

    return pdf.to_pdf()


def render_page(message: str = "", message_type: str = "success") -> bytes:
    rows = load_tv_data()
    categories_html = "".join(f'<option value="{html.escape(c)}">{html.escape(c)}</option>' for c in disease_categories(rows))
    message_html = ""
    if message:
        message_html = f'<section class="messages"><div class="message {message_type}">{html.escape(message)}</div></section>'

    body = f"""<!doctype html>
<html lang='zh-CN'>
<head>
<meta charset='utf-8'>
<title>TI/TV Portal</title>
<link rel='stylesheet' href='/static/style.css'>
</head>
<body>
<header><h1>Target Identification & Target Validation Portal</h1><p>用于 TI 请求提交与 TV 报告自动生成。</p></header>
{message_html}
<main class='grid'>
<section class='card'>
<h2>1) TI：Target Identification Request</h2>
<form method='post' action='/submit_ti'>
<label>疾病背景<textarea name='disease_background' required></textarea></label>
<label>目标机制<textarea name='target_mechanism' required></textarea></label>
<label>数据需求<textarea name='data_requirements' required></textarea></label>
<label>优先级<select name='priority' required>
<option value='P0'>P0 (紧急)</option><option value='P1'>P1 (高)</option><option value='P2'>P2 (中)</option><option value='P3'>P3 (低)</option></select></label>
<label>输出形式<input type='text' name='output_format' required></label>
<button type='submit'>提交 TI Request</button></form>
<p class='hint'>提交后数据会存放在后台数据库，等待人工处理。</p>
</section>
<section class='card'>
<h2>2) TV：Target Validation</h2>
<form method='post' action='/submit_tv'>
<label>Select Target（单基因/双基因）<input type='text' name='target' placeholder='EGFR 或 KRAS+TP53' required></label>
<label>Select Disease Category<select name='disease_category' required>{categories_html}</select></label>
<button type='submit'>生成 TV PDF 报告</button></form>
<p class='hint'>后台使用本地 WES/RNA/scRNA 数据，返回表达与突变结果图。</p>
</section></main></body></html>"""
    return body.encode("utf-8")


class Handler(BaseHTTPRequestHandler):
    def do_GET(self) -> None:  # noqa: N802
        if self.path == "/":
            body = render_page()
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        if self.path == "/static/style.css":
            css = CSS_FILE.read_bytes()
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/css; charset=utf-8")
            self.send_header("Content-Length", str(len(css)))
            self.end_headers()
            self.wfile.write(css)
            return

        self.send_error(HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:  # noqa: N802
        content_length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(content_length).decode("utf-8")
        form = {k: v[0] for k, v in parse_qs(raw).items()}

        if self.path == "/submit_ti":
            self._handle_submit_ti(form)
            return
        if self.path == "/submit_tv":
            self._handle_submit_tv(form)
            return

        self.send_error(HTTPStatus.NOT_FOUND)

    def _handle_submit_ti(self, form: dict[str, str]) -> None:
        with sqlite3.connect(DB_PATH) as conn:
            conn.execute(
                """
                INSERT INTO ti_requests (
                    disease_background, target_mechanism, data_requirements,
                    priority, output_format, created_at
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    form.get("disease_background", ""),
                    form.get("target_mechanism", ""),
                    form.get("data_requirements", ""),
                    form.get("priority", ""),
                    form.get("output_format", ""),
                    datetime.utcnow().isoformat(timespec="seconds"),
                ),
            )
        body = render_page("TI Request 已提交，数据已存入后台等待人工处理。")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _handle_submit_tv(self, form: dict[str, str]) -> None:
        raw_target = form.get("target", "").upper().replace(" ", "")
        disease = form.get("disease_category", "")
        targets = [t for t in raw_target.replace("+", ",").split(",") if t]
        if len(targets) not in (1, 2):
            body = render_page("Select Target 仅支持单基因或双基因组合。", "error")
            self.send_response(HTTPStatus.BAD_REQUEST)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        summary = summarize_tv(load_tv_data(), targets, disease)
        if not summary:
            body = render_page("当前本地数据中没有匹配的 target/disease 结果。", "error")
            self.send_response(HTTPStatus.BAD_REQUEST)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        pdf_data = build_pdf(summary, disease, targets)
        filename = f"TV_Report_{disease}_{'_'.join(targets)}_{uuid4().hex[:8]}.pdf"
        (REPORT_DIR / filename).write_bytes(pdf_data)

        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/pdf")
        self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.send_header("Content-Length", str(len(pdf_data)))
        self.end_headers()
        self.wfile.write(pdf_data)


def run() -> None:
    init_storage()
    server = ThreadingHTTPServer(("0.0.0.0", 8000), Handler)
    print("Server started at http://0.0.0.0:8000")
    server.serve_forever()


if __name__ == "__main__":
    run()
