import argparse
import os
import subprocess
from pathlib import Path

from pypdf import PdfReader


def main():
    parser = argparse.ArgumentParser(description="Optimise a compiled supplement without downsampling; preserve the original.")
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--ghostscript", type=Path, default=Path("C:/texlive/2025/tlpkg/tlgs/bin/gswin64c.exe"))
    arguments = parser.parse_args()
    source = arguments.input.resolve()
    destination = (arguments.output or source.with_name(source.stem + "_optimised.pdf")).resolve()
    if source == destination:
        parser.error("The output must differ from the original PDF.")
    reader = PdfReader(source)
    environment = os.environ.copy()
    installation = arguments.ghostscript.resolve().parent.parent
    if (installation / "Resource/Init/gs_init.ps").is_file():
        environment["GS_LIB"] = os.pathsep.join(str(installation / folder) for folder in ["Resource/Init", "lib", "Resource", "kanji"])
    subprocess.run([
        str(arguments.ghostscript), "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7",
        "-dNOPAUSE", "-dBATCH", "-dSAFER", "-dDetectDuplicateImages=true",
        "-dCompressFonts=true", "-dSubsetFonts=true", "-dDownsampleColorImages=false",
        "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false",
        f"-sOutputFile={destination}", str(source),
    ], check=True, env=environment)
    result = PdfReader(destination)
    if len(result.pages) != len(reader.pages):
        raise RuntimeError("Page count changed.")
    for original_page, optimised_page in zip(reader.pages, result.pages):
        if original_page.mediabox != optimised_page.mediabox:
            raise RuntimeError("Page dimensions changed.")
    original_size = source.stat().st_size
    result_size = destination.stat().st_size
    print(f"Original: {original_size:,} bytes")
    print(f"Optimised: {result_size:,} bytes ({100 * (1 - result_size / original_size):.1f}% reduction)")
    print(f"Verified {len(result.pages)} pages and unchanged page dimensions. Visually check the output before submission.")
    print(destination)


if __name__ == "__main__":
    main()
