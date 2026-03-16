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


def percentile(sorted_values: list[float], p: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    k = (len(sorted_values) - 1) * p
    f = int(k)
    c = min(f + 1, len(sorted_values) - 1)
    if f == c:
        return sorted_values[f]
    return sorted_values[f] + (sorted_values[c] - sorted_values[f]) * (k - f)


def summarize_tv(rows: list[dict[str, str]], targets: list[str], disease: str) -> dict[str, object]:
    selected = [
        row
        for row in rows
        if row["disease_category"] == disease and row["gene"].upper() in targets
    ]

    grouped: dict[str, list[dict[str, str]]] = {}
    for row in selected:
        grouped.setdefault(row["gene"].upper(), []).append(row)

    summary: list[dict[str, float | int | str]] = []
    boxplot_data: dict[str, dict[str, float]] = {}
    oncoplot_samples: list[str] = sorted(
        {
            row["sample_id"]
            for row in selected
            if "WES" in row["data_type"].upper()
        }
    )
    oncoplot_matrix: dict[tuple[str, str], int] = {}

    for gene, items in sorted(grouped.items()):
        expr_values = [float(i["expression"]) for i in items if "RNA" in i["data_type"].upper()]
        mut_values = [int(i["mutation_status"]) for i in items if "WES" in i["data_type"].upper()]

        if expr_values:
            sorted_expr = sorted(expr_values)
            boxplot_data[gene] = {
                "min": sorted_expr[0],
                "q1": percentile(sorted_expr, 0.25),
                "median": percentile(sorted_expr, 0.50),
                "q3": percentile(sorted_expr, 0.75),
                "max": sorted_expr[-1],
            }

        for item in items:
            if "WES" in item["data_type"].upper():
                oncoplot_matrix[(item["sample_id"], gene)] = int(item["mutation_status"])

        mean_expression = sum(expr_values) / len(expr_values) if expr_values else 0.0
        mutation_rate = (100 * sum(mut_values) / len(mut_values)) if mut_values else 0.0

        summary.append(
            {
                "gene": gene,
                "mean_expression": mean_expression,
                "mutation_rate": mutation_rate,
                "sample_count": len(items),
            }
        )

    return {
        "summary": summary,
        "boxplot_data": boxplot_data,
        "oncoplot_samples": oncoplot_samples,
        "oncoplot_matrix": oncoplot_matrix,
    }


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
        objects = [
            b"1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\n",
            b"2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\n",
            b"3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >> endobj\n",
            b"4 0 obj << /Type /Font /Subtype /Type1 /BaseFont /Helvetica >> endobj\n",
            f"5 0 obj << /Length {len(content)} >> stream\n".encode("latin-1") + content + b"endstream endobj\n",
        ]

        out = io.BytesIO()
        out.write(b"%PDF-1.4\n")
        offsets = [0]
        for obj in objects:
            offsets.append(out.tell())
            out.write(obj)

        xref_pos = out.tell()
        out.write(f"xref\n0 {len(objects) + 1}\n".encode())
        out.write(b"0000000000 65535 f \n")
        for off in offsets[1:]:
            out.write(f"{off:010d} 00000 n \n".encode())
        out.write(f"trailer << /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF".encode())
        return out.getvalue()


def draw_rna_boxplot(pdf: TinyPDF, boxplot_data: dict[str, dict[str, float]]) -> None:
    pdf.text(40, 390, 11, "RNA Expression Boxplot")
    axis_x0, axis_y0 = 40, 220
    axis_x1, axis_y1 = 270, 360
    pdf.line(axis_x0, axis_y0, axis_x1, axis_y0)
    pdf.line(axis_x0, axis_y0, axis_x0, axis_y1)

    if not boxplot_data:
        pdf.text(50, 330, 9, "No RNA data available for selected target/disease.")
        return

    all_values = [v for gene in boxplot_data.values() for v in gene.values()]
    min_v = min(all_values)
    max_v = max(all_values)
    scale = (axis_y1 - axis_y0 - 10) / (max_v - min_v if max_v > min_v else 1)

    for idx, (gene, stat) in enumerate(sorted(boxplot_data.items())):
        cx = 75 + idx * 95

        def y(v: float) -> float:
            return axis_y0 + 5 + (v - min_v) * scale

        y_min, y_q1, y_med, y_q3, y_max = y(stat["min"]), y(stat["q1"]), y(stat["median"]), y(stat["q3"]), y(stat["max"])
        pdf.line(cx, y_min, cx, y_q1)
        pdf.line(cx, y_q3, cx, y_max)
        pdf.rect(cx - 18, y_q1, 36, max(y_q3 - y_q1, 1))
        pdf.line(cx - 18, y_med, cx + 18, y_med)
        pdf.line(cx - 10, y_min, cx + 10, y_min)
        pdf.line(cx - 10, y_max, cx + 10, y_max)
        pdf.text(cx - 12, axis_y0 - 14, 9, gene)


def draw_wes_oncoplot(
    pdf: TinyPDF,
    targets: list[str],
    samples: list[str],
    matrix: dict[tuple[str, str], int],
) -> None:
    pdf.text(320, 390, 11, "WES Mutation Oncoplot")
    x0, y0 = 320, 220
    cell_w, cell_h = 24, 16

    if not samples:
        pdf.text(330, 330, 9, "No WES data available for selected target/disease.")
        return

    shown_samples = samples[:8]
    shown_targets = targets[:2]

    for cidx, gene in enumerate(shown_targets):
        pdf.text(x0 + 45 + cidx * cell_w, y0 + (len(shown_samples) + 1) * cell_h + 8, 8, gene)

    for ridx, sample_id in enumerate(shown_samples):
        y = y0 + (len(shown_samples) - ridx) * cell_h
        pdf.text(x0, y + 4, 7, sample_id)
        for cidx, gene in enumerate(shown_targets):
            x = x0 + 42 + cidx * cell_w
            status = matrix.get((sample_id, gene), 0)
            if status == 1:
                pdf.rgb(0.82, 0.18, 0.18)
            else:
                pdf.rgb(0.90, 0.90, 0.90)
            pdf.filled_rect(x, y, cell_w - 2, cell_h - 2)
            pdf.rgb(0, 0, 0)
            pdf.rect(x, y, cell_w - 2, cell_h - 2)

    ly = y0 - 22
    pdf.rgb(0.82, 0.18, 0.18)
    pdf.filled_rect(x0 + 4, ly, 10, 10)
    pdf.rgb(0, 0, 0)
    pdf.rect(x0 + 4, ly, 10, 10)
    pdf.text(x0 + 18, ly + 2, 8, "Mutated")
    pdf.rgb(0.90, 0.90, 0.90)
    pdf.filled_rect(x0 + 95, ly, 10, 10)
    pdf.rgb(0, 0, 0)
    pdf.rect(x0 + 95, ly, 10, 10)
    pdf.text(x0 + 109, ly + 2, 8, "Wild-type / missing")


def build_pdf(aggregated: dict[str, object], disease: str, targets: list[str]) -> bytes:
    summary = aggregated["summary"]
    boxplot_data = aggregated["boxplot_data"]
    oncoplot_samples = aggregated["oncoplot_samples"]
    oncoplot_matrix = aggregated["oncoplot_matrix"]

    pdf = TinyPDF()
    y = 810
    pdf.text(40, y, 16, "Target Validation (TV) Report")
    y -= 22
    pdf.text(40, y, 11, f"Disease Category: {disease}")
    y -= 16
    pdf.text(40, y, 11, f"Targets: {', '.join(targets)}")
    y -= 16
    pdf.text(40, y, 10, f"Generated(UTC): {datetime.utcnow().isoformat(timespec='seconds')}")

    y -= 22
    pdf.text(40, y, 9, "scRNA featureplot: not implemented yet (planned for future release).")

    y -= 22
    pdf.text(40, y, 12, "Summary Table")
    y -= 14
    pdf.text(40, y, 9, "Gene")
    pdf.text(120, y, 9, "Mean RNA Exp")
    pdf.text(240, y, 9, "WES Mutation Rate(%)")
    pdf.text(390, y, 9, "N Samples")
    y -= 8
    pdf.line(40, y, 560, y)

    for row in summary:  # type: ignore[assignment]
        y -= 16
        row = row  # type: ignore[no-redef]
        pdf.text(40, y, 9, str(row["gene"]))
        pdf.text(120, y, 9, f"{row['mean_expression']:.2f}")
        pdf.text(240, y, 9, f"{row['mutation_rate']:.1f}")
        pdf.text(390, y, 9, str(row["sample_count"]))

    draw_rna_boxplot(pdf, boxplot_data)  # type: ignore[arg-type]
    draw_wes_oncoplot(pdf, targets, oncoplot_samples, oncoplot_matrix)  # type: ignore[arg-type]

    return pdf.to_pdf()


def render_page(message: str = "", message_type: str = "success") -> bytes:
    rows = load_tv_data()
    categories_html = "".join(
        f'<option value="{html.escape(c)}">{html.escape(c)}</option>'
        for c in disease_categories(rows)
    )
    message_html = ""
    if message:
        message_html = (
            f'<section class="messages"><div class="message {message_type}">'
            f"{html.escape(message)}</div></section>"
        )

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
<p class='hint'>RNA 将输出 boxplot；WES 将输出 mutation oncoplot；scRNA featureplot 暂不支持。</p>
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

        aggregated = summarize_tv(load_tv_data(), targets, disease)
        if not aggregated["summary"]:
            body = render_page("当前本地数据中没有匹配的 target/disease 结果。", "error")
            self.send_response(HTTPStatus.BAD_REQUEST)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        pdf_data = build_pdf(aggregated, disease, targets)
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
