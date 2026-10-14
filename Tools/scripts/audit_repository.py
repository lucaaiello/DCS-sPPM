import argparse
import json
import re
from pathlib import Path


PROJECT = next(parent for parent in Path(__file__).resolve().parents
               if (parent / "R" / "project_paths.R").is_file())


def main():
    parser = argparse.ArgumentParser(description="Read-only manuscript dependency audit")
    parser.add_argument("--require-results", action="store_true")
    parser.add_argument("--require-manuscripts", action="store_true")
    arguments = parser.parse_args()
    errors = []
    migration = json.loads((PROJECT / "docs/layout_migration_2026-10-13.json").read_text(encoding="utf-8"))
    for category in ("scripts", "results"):
        for before, after in migration[category].items():
            if not (PROJECT / after).exists() and (category == "scripts" or arguments.require_results):
                errors.append(f"Missing migrated {category} path: {after}")
            if (PROJECT / before).exists():
                errors.append(f"Obsolete duplicate {category} path: {before}")
    documents = [PROJECT / "README.md", PROJECT / "Data/README.md"]
    documents.extend((PROJECT / "Data/preparation").rglob("*.md"))
    documents.extend((PROJECT / "docs").glob("*.md"))
    for document in documents:
        source = document.read_text(encoding="utf-8-sig")
        for reference in re.findall(r"\[[^\]]*\]\(([^)]+)\)", source):
            if reference.startswith(("http:", "https:", "#")):
                continue
            if not (document.parent / reference.split("#", 1)[0]).exists():
                errors.append(f"Missing documentation link: {document.relative_to(PROJECT)}: {reference}")
    figures = {}
    missing_manuscripts = []
    for relative in ("Luc_i_a/main_v2.tex", "Luc_i_a_JRSSC/main.tex",
                     "Luc_i_a_JRSSC/supplementary.tex"):
        document = PROJECT / relative
        if not document.is_file():
            missing_manuscripts.append(relative)
            if arguments.require_manuscripts:
                errors.append(f"Missing manuscript: {relative}")
            continue
        source = re.sub(r"(?m)(?<!\\)%.*$", "", document.read_text(encoding="utf-8-sig"))
        references = re.findall(
            r"\\(?:includegraphics\*?|begin\{overpic\})\s*(?:\[[^\]]*\]\s*)?\{([^}]+)\}",
            source,
        )
        figures[relative] = references
        for reference in references:
            path = document.parent / reference
            if not path.suffix:
                path = path.with_suffix(".pdf")
            if not path.exists():
                errors.append(f"Missing figure: {relative}: {reference}")
        for reference in re.findall(r"\\(?:input|include)\{([^}]+)\}", source):
            if not (document.parent / reference).with_suffix(".tex").exists():
                errors.append(f"Missing included TeX file: {relative}: {reference}")
        for reference in re.findall(r"\\bibliography\{([^}]+)\}", source):
            for bibliography in reference.split(","):
                if not (document.parent / bibliography.strip()).with_suffix(".bib").exists():
                    errors.append(f"Missing bibliography: {relative}: {bibliography}")

    code_directories = [PROJECT / "R", PROJECT / "src"]
    code_directories.extend([PROJECT / "Data/preparation", PROJECT / "Data/preparation/manifests"])
    code_directories.extend(PROJECT.glob("*/scripts"))
    code_directories.extend(PROJECT.glob("*/*/scripts"))
    code_files = [path for directory in code_directories
                  for path in directory.glob("*") if path.is_file()
                  and path.suffix in (".R", ".cpp", ".py", ".ps1")]
    known_names = {path.name for path in code_files}
    for path in code_files:
        if path.suffix not in (".R", ".ps1"):
            continue
        source = path.read_text(encoding="utf-8-sig")
        for reference in re.findall(r'["\']([^"\'\n]+\.(?:R|cpp))["\']', source):
            name = reference.replace("\\", "/").split("/")[-1]
            if re.fullmatch(r"[A-Za-z0-9_.-]+\.(?:R|cpp)", name) and name not in known_names:
                errors.append(f"Missing code dependency: {path.relative_to(PROJECT)}: {reference}")

    result_root = PROJECT / "Real_data_analysis/model_comparison/results/production_30000_unique_locations"
    inputs = [PROJECT / "Data" / name for name in (
        "pm10_combined_2024.rds", "era5_combined_2024.rds", "pm10_unique_locations_2024.rds",
        "production_30000_online_grid_preparation.rds", "production_10000_single_chain_preparation.rds")]
    inputs.extend(result_root / model / "pilot_fit.rds" for model in (
        "dcs_sppm", "ppm_common_shock", "k1_common_shock", "sppm_no_common_shock"))
    inputs.extend([
        result_root / "dcs_sppm/maps/clustering/vi_partition.rds",
        result_root / "dcs_sppm/predictive_temporal_regime/predictive_temporal_regime_summary.rds",
        PROJECT / "Simulations/clustering_recovery/results/clustering_recovery_100rep/parameter_summary.csv",
        PROJECT / "Simulations/spatial_prediction/results/prediction_experiment_100rep/prediction_summary.csv",
        PROJECT / "Simulations/model_comparison/results/production/simulation_comparison_table.csv",
        PROJECT / "Simulations/gp_nngp_benchmark/results/summary.csv",
    ])
    missing_inputs = [str(path.relative_to(PROJECT)) for path in inputs if not path.exists()]
    if arguments.require_results:
        errors.extend(f"Missing frozen input/result: {path}" for path in missing_inputs)
    print(json.dumps({
        "project": str(PROJECT), "manuscript_figures": figures,
        "retained_code_files": len(code_files), "missing_local_inputs": missing_inputs,
        "missing_local_manuscripts": missing_manuscripts,
        "errors": sorted(set(errors)),
        "external_data": "Raw EEA/ERA5 downloads are external; the frozen-cohort route uses Data/preparation manifests without parent caches or historical fits. Historical selection audits need parent caches.",
        "scientific_validation": "Not performed: paths are checked, not numerical results or convergence.",
    }, indent=2))
    raise SystemExit(1 if errors else 0)


if __name__ == "__main__":
    main()
