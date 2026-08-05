from pypdf import PdfReader
from pathlib import Path
import sys

def pdf_to_txt(pdf_path: str, txt_path: str | None = None):
    pdf_file = Path(pdf_path)
    if not pdf_file.exists():
        raise FileNotFoundError(f"找不到文件: {pdf_file}")

    if txt_path is None:
        txt_file = pdf_file.with_suffix(".txt")
    else:
        txt_file = Path(txt_path)

    reader = PdfReader(str(pdf_file))
    pages_text = []
    for i, page in enumerate(reader.pages, start=1):
        text = page.extract_text() or ""
        pages_text.append(f"--- Page {i} ---\n{text}\n")

    txt_file.write_text("\n".join(pages_text), encoding="utf-8")
    print(f"已生成: {txt_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python pdf_to_txt.py <pdf路径> [txt路径]")
        sys.exit(1)
    pdf_to_txt(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)