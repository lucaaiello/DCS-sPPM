import argparse
import json
import shutil
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Copy Git-visible reproduction sources to an existing repository, without committing or pushing.")
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    arguments = parser.parse_args()
    source = Path(__file__).resolve().parents[2]
    destination = arguments.destination.resolve()
    if destination == source or source in destination.parents or destination in source.parents:
        parser.error("The destination must be a separate checkout.")
    if not (destination / ".git").exists():
        parser.error("The destination must already be a Git repository.")
    output = subprocess.check_output([
        "git", "-C", str(source.parent), "ls-files", "--cached", "--others",
        "--exclude-standard", "-z", "--", source.name,
    ])
    relative_paths = sorted(set(Path(path.decode("utf-8")).relative_to(source.name)
                                for path in output.split(b"\0") if path))
    forbidden = {".rds", ".rda", ".rdata", ".pdf", ".png", ".parquet", ".nc", ".zip", ".dll", ".o"}
    for relative in relative_paths:
        if relative.suffix.lower() in forbidden or any(part in ("results", ".git", "Luc_i_a", "Luc_i_a_JRSSC", "raw") for part in relative.parts):
            raise RuntimeError(f"Non-code artifact in packaging list: {relative}")
        if (source / relative).stat().st_size > 10_000_000:
            raise RuntimeError(f"Unexpected large input: {relative}")
        if (source / relative).is_symlink():
            raise RuntimeError(f"Refusing symbolic link: {relative}")
        target = (destination / relative).resolve()
        if destination not in target.parents:
            raise RuntimeError(f"Destination escapes checkout: {relative}")
    print(json.dumps({"source": str(source), "destination": str(destination),
                      "files": len(relative_paths), "bytes": sum((source / path).stat().st_size for path in relative_paths),
                      "dry_run": arguments.dry_run}, indent=2))
    if arguments.dry_run:
        return
    status = subprocess.check_output(["git", "-C", str(destination), "status", "--porcelain"])
    if status.strip():
        raise RuntimeError("Destination has local changes. Review them before synchronising; nothing was copied.")
    for relative in relative_paths:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / relative, target)
    print("Source files synchronised. Existing extra files, .git, remote and history are untouched. No commit or push performed.")


if __name__ == "__main__":
    main()
