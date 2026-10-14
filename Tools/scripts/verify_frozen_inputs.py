import argparse
import re
from pathlib import Path

from package_frozen_inputs import INPUTS, sha256_file


def main():
    project = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Verify the three downloaded analysis assets against the release SHA256SUMS.")
    parser.add_argument("--input-dir", type=Path, default=project / "Data")
    parser.add_argument("--checksums", type=Path, required=True)
    arguments = parser.parse_args()
    expected = {}
    for line in arguments.checksums.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{64})\s+\*?([^/\\]+)", line)
        if not match or match[2] not in INPUTS or match[2] in expected:
            parser.error(f"Invalid checksum entry: {line}")
        expected[match[2]] = match[1].lower()
    if set(expected) != set(INPUTS):
        parser.error("SHA256SUMS must list exactly the three production input files.")
    for name in INPUTS:
        path = arguments.input_dir / name
        if not path.is_file() or sha256_file(path) != expected[name]:
            raise SystemExit(f"Missing or mismatched input: {path}")
        print(f"OK: {name}")
    print("All release checksums match. Run the R input validator next; no fitting was performed.")


if __name__ == "__main__":
    main()
