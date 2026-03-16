# TI/TV Web Portal

一个用于 **Target Identification (TI)** 和 **Target Validation (TV)** 的轻量网页示例（仅使用 Python 标准库实现）。

## 功能

- TI 请求提交：
  - 前端表单收集疾病背景、目标机制、数据需求、优先级、输出形式。
  - 后端将请求写入 SQLite 数据库并标记为 `pending`，等待人工处理。
- TV 验证报告：
  - 基于本地存储的 WES/RNA/scRNA 示例数据。
  - 支持单基因或双基因组合（例如 `EGFR` 或 `KRAS+TP53`）。
  - 按疾病类别过滤后生成 PDF 报告，包含统计表、RNA expression boxplot、WES mutation oncoplot。
  - scRNA featureplot 功能暂时搁置，后续版本再增加。

## 启动

```bash
python app.py
```

访问 `http://localhost:8000`。

## 目录

- `app.py`：HTTP 服务、数据库写入、TV 统计与 PDF 生成逻辑
- `static/style.css`：页面样式
- `data/tv_mock_data.csv`：本地示例数据
- `ti_tv.db`：运行后自动生成的 SQLite 数据库
- `reports/`：运行后生成的 PDF 报告
