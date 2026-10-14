import argparse
import csv
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


PROJECT = next(parent for parent in Path(__file__).resolve().parents
               if (parent / "R" / "project_paths.R").is_file())


@dataclass
class Figure:
    document: Path
    number: str
    label: str
    image: Path
    graphic: str
    panel_macro: str
    panels: list[str]


def strip_comments(source):
    return re.sub(r"(?m)(?<!\\)%.*$", "", source)


def read_group(source, start, opening="{", closing="}"):
    if start >= len(source) or source[start] != opening:
        raise ValueError(f"Expected {opening!r} in LaTeX source.")
    depth = 1
    position = start + 1
    while position < len(source):
        if source[position] == "\\":
            position += 2
            continue
        if source[position] == opening:
            depth += 1
        elif source[position] == closing:
            depth -= 1
            if depth == 0:
                return source[start + 1:position], position + 1
        position += 1
    raise ValueError("Unbalanced LaTeX group.")


def image_argument(graphic):
    command = re.search(r"\\(?:begin\{overpic\}|includegraphics\*?)", graphic)
    if command is None:
        raise ValueError("No supported image command found.")
    position = command.end()
    while position < len(graphic) and graphic[position].isspace():
        position += 1
    if graphic[position] == "[":
        _, position = read_group(graphic, position, "[", "]")
    while position < len(graphic) and graphic[position].isspace():
        position += 1
    filename, end = read_group(graphic, position)
    return filename.strip(), end


def panel_definition(source):
    definition = re.search(
        r"\\(?:newcommand|renewcommand)\s*\{\\jrsscpanellabel\}\s*\[1\]\s*", source
    )
    if definition is None:
        return ""
    body, _ = read_group(source, definition.end())
    return "\\newcommand{\\jrsscpanellabel}[1]{" + body + "}"


def collect_figures(directory, documents):
    figures = []
    counts = {"main": 0, "supplementary": 0}
    for filename in documents:
        document = (directory / filename).resolve()
        if not document.is_file():
            raise FileNotFoundError(f"Missing manuscript source: {document}")
        source = strip_comments(document.read_text(encoding="utf-8-sig"))
        blocks = list(re.finditer(
            r"\\begin\{(figure\*?)\}(.*?)\\end\{\1\}", source, re.S
        ))
        if not blocks:
            raise ValueError(f"No figure environments found in {document.name}.")
        if len(blocks) != len(re.findall(r"\\begin\{figure\*?\}", source)):
            raise ValueError(f"Unbalanced figure environments in {document.name}.")
        macro = panel_definition(source)
        category = "main" if document.stem == "main" else "supplementary"
        for block in blocks:
            body = block.group(2)
            labels = re.findall(r"\\label\{(fig:[^}]+)\}", body)
            if len(labels) != 1:
                raise ValueError(f"Expected one fig: label in a figure in {document.name}.")
            overpics = re.findall(r"\\begin\{overpic\}.*?\\end\{overpic\}", body, re.S)
            if overpics:
                if len(overpics) != 1:
                    raise ValueError(f"Unsupported multiple overpic blocks: {labels[0]}")
                graphic = overpics[0]
            else:
                images = list(re.finditer(r"\\includegraphics\*?", body))
                if len(images) != 1:
                    raise ValueError(f"Expected one complete image in {labels[0]}.")
                graphic = body[images[0].start():]
                _, end = image_argument(graphic)
                graphic = graphic[:end]
            image_name, _ = image_argument(graphic)
            if "\\" in image_name or "{" in image_name:
                raise ValueError(f"Use a literal image filename in {labels[0]}.")
            image = document.parent / image_name
            if not image.is_file():
                image = document.parent / "fig" / image_name
            if not image.is_file():
                raise FileNotFoundError(f"Missing image for {labels[0]}: {image_name}")
            panels = re.findall(r"\\jrsscpanellabel\{([^}]+)\}", graphic)
            if panels and not macro:
                raise ValueError(f"Missing jrsscpanellabel definition in {document.name}.")
            counts[category] += 1
            number = str(counts[category])
            if category == "supplementary":
                number = "S" + number
            figures.append(Figure(document, number, labels[0], image.resolve(),
                                  graphic, macro, panels))
    stems = [figure.image.stem.casefold() for figure in figures]
    if len(stems) != len(set(stems)):
        raise ValueError("Repeated image stems would collide; give the source figures unique filenames.")
    return figures


def standalone_source(figure, width_mm):
    graphic = figure.graphic
    image_name, _ = image_argument(graphic)
    graphic = graphic.replace("{" + image_name + "}",
                              "{\\detokenize{" + figure.image.as_posix() + "}}", 1)
    return (
        "\\documentclass[border={2mm 2mm 2mm 6mm}]{standalone}\n"
        "\\pdfminorversion=7\n"
        "\\usepackage[T1]{fontenc}\n"
        "\\usepackage{graphicx}\n"
        "\\usepackage[percent]{overpic}\n"
        "\\renewcommand{\\normalsize}{\\sffamily\\fontsize{9bp}{11.5pt}\\selectfont}\n"
        f"\\setlength{{\\textwidth}}{{{width_mm:g}mm}}\n"
        + figure.panel_macro + "\n\\begin{document}\n"
        + graphic + "\n\\end{document}\n"
    )


def export_figures(figures, directory, executable, width_mm):
    targets = [directory / (figure.image.stem + ".pdf") for figure in figures]
    if any(target.resolve() == figure.image for target, figure in zip(targets, figures)):
        raise ValueError("Output would overwrite a source image; use a separate export directory.")
    directory.mkdir(parents=True, exist_ok=True)
    build = directory / "_build"
    build.mkdir(exist_ok=True)
    records = []
    for figure, target in zip(figures, targets):
        stem = figure.image.stem
        tex = build / (stem + ".tex")
        tex.write_text(standalone_source(figure, width_mm), encoding="utf-8")
        print(f"Exporting figure {figure.number}: {target.name}", flush=True)
        result = subprocess.run(
            [executable, "-no-shell-escape", "-interaction=nonstopmode", "-halt-on-error",
             "-output-directory=" + str(build), str(tex)],
            cwd=figure.document.parent, capture_output=True, text=True, errors="replace",
        )
        (build / (stem + "_console.txt")).write_text(result.stdout + result.stderr,
                                                    encoding="utf-8")
        pdf = build / (stem + ".pdf")
        log = build / (stem + ".log")
        if result.returncode != 0 or not pdf.is_file() or not log.is_file():
            raise RuntimeError(f"Compilation failed for {figure.label}; inspect {build / (stem + '_console.txt')}")
        compiled = re.search(r"Output written on.*?\((\d+)\s+page",
                             log.read_text(encoding="utf-8", errors="replace"), re.S)
        if compiled is None or compiled.group(1) != "1":
            raise RuntimeError(f"Expected one PDF page for {figure.label}; inspect {log}")
        records.append({
            "document": figure.document.name, "figure_number": figure.number,
            "latex_label": figure.label, "source_image": figure.image.name,
            "panel_labels": ",".join(figure.panels), "exported_pdf": target.name,
            "width_mm": width_mm, "file_bytes": pdf.stat().st_size,
        })
    for target, figure in zip(targets, figures):
        shutil.copy2(build / (figure.image.stem + ".pdf"), target)
    manifest = directory / "figure_export_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f"Created {len(figures)} standalone PDFs in {directory}")
    print(f"Figure numbering and source mapping: {manifest}")


def main():
    parser = argparse.ArgumentParser(description="Export manuscript figures with their existing LaTeX panel labels.")
    parser.add_argument("--manuscript-dir", type=Path, default=PROJECT / "Luc_i_a_JRSSC")
    parser.add_argument("--documents", nargs="+", default=["main.tex", "supplementary.tex"])
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--pdflatex", default="pdflatex")
    parser.add_argument("--width-mm", type=float, default=157.0)
    parser.add_argument("--dry-run", action="store_true")
    arguments = parser.parse_args()
    if not 0 < arguments.width_mm < 1000:
        parser.error("Figure width must be positive and below 1,000 mm.")
    manuscript = arguments.manuscript_dir.resolve()
    output = arguments.output_dir.resolve() if arguments.output_dir else manuscript / "fig_exported"
    figures = collect_figures(manuscript, arguments.documents)
    if arguments.dry_run:
        for figure in figures:
            print(f"{figure.document.name}: figure {figure.number}, {figure.label}, "
                  f"panels={','.join(figure.panels) or 'none'} -> {figure.image.stem}.pdf")
            standalone_source(figure, arguments.width_mm)
        print(f"Validated {len(figures)} figures. Planned output: {output}. No files written or compiled.")
        return
    executable = shutil.which(arguments.pdflatex)
    if executable is None:
        parser.error("pdflatex was not found; supply its full path with --pdflatex.")
    export_figures(figures, output, executable, arguments.width_mm)


if __name__ == "__main__":
    main()
