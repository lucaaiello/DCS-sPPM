import argparse
import csv
import shutil
import subprocess
from pathlib import Path

from pypdf import PageObject, PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject, ContentStream, NameObject


PROJECT = next(parent for parent in Path(__file__).resolve().parents
               if (parent / "R" / "project_paths.R").is_file())
DEFAULT_SOURCE = PROJECT / "Luc_i_a" / "fig"
DEFAULT_OUTPUT = PROJECT / "Luc_i_a_JRSSC" / "fig"


def read_page(path):
    reader = PdfReader(path)
    if len(reader.pages) != 1:
        raise ValueError(f"Expected a single-page figure: {path}")
    return reader.pages[0]


def write_pdf(page, path, labels):
    writer = PdfWriter()
    writer.add_page(page)
    writer.compress_identical_objects(remove_duplicates=True, remove_unreferenced=True)
    writer.pages[0].compress_content_streams()
    writer.add_metadata({"/Title": path.stem, "/Creator": "DCS-sPPM figure packaging"})
    with path.open("wb") as stream:
        writer.write(stream)


def combine_panels(paths, output):
    pages = [read_page(path) for path in paths]
    panel_width = max(float(page.mediabox.width) for page in pages)
    panel_height = max(float(page.mediabox.height) for page in pages)
    gutter = 12
    label_band = 22
    combined = PageObject.create_blank_page(
        width=len(pages) * panel_width + (len(pages) - 1) * gutter,
        height=panel_height + label_band,
    )
    labels = []
    for index, page in enumerate(pages):
        scale = min(
            panel_width / float(page.mediabox.width),
            panel_height / float(page.mediabox.height),
        )
        horizontal = index * (panel_width + gutter)
        vertical = (panel_height - float(page.mediabox.height) * scale) / 2
        combined.merge_transformed_page(
            page, Transformation().scale(scale).translate(horizontal, vertical)
        )
        labels.append((chr(ord("A") + index), horizontal + 8, panel_height + 6))
    write_pdf(combined, output, labels)


def render_png(pdf, executable):
    subprocess.run(
        [executable, "-q", "-singlefile", "-r", "350", "-png", str(pdf), str(pdf.with_suffix(""))],
        check=True,
    )


def combine_response_rows(paths, output, upper_row):
    pages = [read_page(path) for path in paths]
    width = float(pages[0].mediabox.width)
    height = float(pages[0].mediabox.height)
    combined = PageObject.create_blank_page(width=width, height=height)
    labels = []
    for index, page in enumerate(pages):
        if (float(page.mediabox.width), float(page.mediabox.height)) != (width, height):
            raise ValueError("Response figures must have identical dimensions")
        source_bottom = height / 2 if upper_row else 0
        page.add_transformation(Transformation().translate(0, -source_bottom))
        bounds = RectangleObject([0, 0, width, height / 2])
        page.mediabox = bounds
        page.cropbox = bounds
        page.trimbox = bounds
        destination_bottom = height / 2 if index == 0 else 0
        combined.merge_transformed_page(
            page, Transformation().translate(0, destination_bottom)
        )
        for column in range(2):
            labels.append((chr(ord("A") + index * 2 + column),
                           column * width / 2 + 10, destination_bottom + height / 2 - 15))
    write_pdf(combined, output, labels)


def package_response_rows(source, output, executable):
    paths = [source / f"y_{date}.pdf" for date in ("2024-04-09", "2024-12-15")]
    for stem, upper_row in (("europe_response_mean_sd_jrssc", True),
                            ("europe_response_limits_jrssc", False)):
        destination = output / f"{stem}.pdf"
        combine_response_rows(paths, destination, upper_row)
        render_png(destination, executable)


def label_comparisons(directory, executable):
    specifications = {
        "comparative_simulation_partition_recovery": (0.055, 0.29, 0.525),
        "comparative_real_data_structural_sensitivity": (0.015, 0.37),
    }
    for stem, positions in specifications.items():
        path = directory / f"{stem}.pdf"
        if not path.exists():
            raise FileNotFoundError(f"Run generate_comparative_figures.R before packaging comparisons: {path}")
        reader = PdfReader(path)
        if reader.metadata.get("/Creator") == "DCS-sPPM figure packaging":
            page = reader.pages[0]
            label_fonts = {name for name, font in page["/Resources"]["/Font"].items()
                           if "BitstreamVeraSans-Bold" in str(font.get_object().get("/BaseFont", ""))}
            content = ContentStream(page.get_contents(), reader)
            active_font = None
            operations = []
            for operands, operator in content.operations:
                if operator == b"Tf":
                    active_font = operands[0]
                if operator in (b"Tj", b"TJ", b"'", b'"') and active_font in label_fonts:
                    continue
                operations.append((operands, operator))
            content.operations = operations
            page[NameObject("/Contents")] = content
            write_pdf(page, path, [])
        if reader.metadata.get("/Creator") != "DCS-sPPM figure packaging":
            page = reader.pages[0]
            width = float(page.mediabox.width)
            height = float(page.mediabox.height)
            labelled = PageObject.create_blank_page(width=width, height=height + 22)
            labelled.merge_page(page)
            labels = [(chr(ord("A") + index), width * position, height + 6)
                      for index, position in enumerate(positions)]
            write_pdf(labelled, path, labels)
        render_png(path, executable)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--pdftoppm", default="pdftoppm")
    parser.add_argument("--comparisons-only", action="store_true")
    parser.add_argument("--responses-only", action="store_true")
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=True)
    if arguments.responses_only:
        package_response_rows(arguments.source, arguments.output, arguments.pdftoppm)
        print(f"Updated response mean/SD and predictive-limit figures in {arguments.output}")
        return
    if arguments.comparisons_only:
        label_comparisons(arguments.output, arguments.pdftoppm)
        print(f"Updated comparison panel labels in {arguments.output}")
        return
    specifications = {
        "europe_pm10_exploratory_jrssc": ["europe_pm10_mean", "europe_pm10_lag1_acf"],
        "vi_partition_jrssc": ["vi_partition_graph_colours", "vi_partition_stability"],
        "vi_regime_parameters_jrssc": ["vi_partition_phi", "vi_partition_tau2"],
        "europe_predictive_temporal_regime_jrssc": [
            "europe_predictive_phi", "europe_predictive_tau", "europe_predictive_regime_uncertainty"
        ],
        "europe_predictive_regime_map_jrssc": [
            "europe_predictive_regime_map", "europe_predictive_regime_entropy"
        ],
    }
    records = []
    package_response_rows(arguments.source, arguments.output, arguments.pdftoppm)
    for stem in ("europe_response_mean_sd_jrssc", "europe_response_limits_jrssc"):
        records.append((stem, " | ".join(str(arguments.source / f"y_{date}.pdf")
                                       for date in ("2024-04-09", "2024-12-15")),
                        "selected response rows assembled as vector PDF"))
    for stem, components in specifications.items():
        paths = [arguments.source / f"{component}.pdf" for component in components]
        output = arguments.output / f"{stem}.pdf"
        combine_panels(paths, output)
        render_png(output, arguments.pdftoppm)
        records.append((stem, " | ".join(str(path) for path in paths), "assembled vector PDF"))
    for stem in ("dcs_sppm_conceptual_scheme", "europe_pm10_temporal_summary"):
        for suffix in (".pdf", ".png"):
            source = arguments.source / f"{stem}{suffix}"
            destination = arguments.output / source.name
            if source.resolve() != destination.resolve():
                shutil.copy2(source, destination)
        records.append((stem, str(arguments.source / f"{stem}.pdf"), "existing single figure"))
    for stem in ("comparative_simulation_partition_recovery", "comparative_real_data_structural_sensitivity"):
        records.append((stem, "Manuscript_figures/scripts/generate_comparative_figures.R",
                        "comparison figure; panel labels added by LaTeX"))
    with (arguments.output / "jrssc_figure_manifest.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("stem", "source", "operation"))
        for stem, sources, operation in records:
            portable_sources = []
            for source in sources.split(" | "):
                path = Path(source)
                try:
                    portable_sources.append(path.relative_to(PROJECT).as_posix())
                except ValueError:
                    portable_sources.append(source)
            writer.writerow((stem, " | ".join(portable_sources), operation))
    label_comparisons(arguments.output, arguments.pdftoppm)
    print(f"Packaged {len(records)} figures in {arguments.output}")


if __name__ == "__main__":
    main()
