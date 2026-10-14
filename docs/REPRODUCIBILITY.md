# Reproduction workflow

Run commands from the project directory. Use `Rscript`/`python` on PATH, or
replace them with installed executable paths. No workflow below was launched
by the repository cleanup. Avoid inherited experimental environment variables;
these examples deliberately identify the final result names.

## 1. Read-only audit

For a fresh clone, first download the three frozen input assets and SHA256SUMS
from [release v1.0.0](https://github.com/lucaaiello/DCS-sPPM/releases/tag/v1.0.0). Follow [input installation](../Data/README.md)
and verify checksums before refitting. The release link must identify the input
version matching the code, not a moving latest download. Raw preprocessing is
not part of this reproduction route. Historical diagnostic inputs are optional.

```powershell
python Tools/scripts/audit_repository.py
Rscript Tools/scripts/audit_source_syntax.R
```

This checks file connections, not scientific correctness. A fresh clone will
report missing local inputs and manuscripts without failing the code-only audit.
After installing the release assets, use `--require-inputs` and run
`Rscript Tools/scripts/validate_preparation_inputs.R` to check the final
data/design/grid without fitting. `--require-results --require-manuscripts`
additionally requires saved results and local manuscripts. Results, figures,
data and both LaTeX folders are excluded from Git.

## 2. Recreate numerical summaries from checkpoints

These run post-processing, not MCMC. The final real-data summary also performs
VI searches; it is not a cheap plotting step. Do not rerun it just to change
map styling: use existing summary CSVs and cached VI instead.

```powershell
Rscript Simulations/clustering_recovery/scripts/summarize_clustering_recovery.R
$env:DCS_PREDICTION_RESULTS_NAME = "prediction_experiment_100rep"
Rscript Simulations/spatial_prediction/scripts/summarize_prediction_experiment.R
$env:DCS_COMPARISON_RESULTS_NAME = "production"
Rscript Simulations/model_comparison/scripts/summarize_model_comparison.R
$env:DCS_REAL_FULL_RESULTS_NAME = "production_30000_unique_locations"
Rscript Real_data_analysis/model_comparison/scripts/summarize_full_analysis.R
```

Consult each script for its environment settings before changing an experiment.
The production DCS simulation checkpoints are reused by `summarize_model_comparison.R` to avoid
refitting the same data solely for model comparison.

## 3. Refresh manuscript graphics using cached output

```powershell
./Manuscript_figures/scripts/refresh_assets.ps1 -RscriptPath "C:/Program Files/R/R-4.4.1/bin/Rscript.exe" -PythonPath "python"
```

This requires the final fit, saved VI and saved predictive-regime summary before
starting, as well as the separately supplied manuscript projects/source assets.
It sets final-cohort paths, skips the full daily PNG sweep, reuses the
cached partition and regime summaries, then composes only the JRSS-C figure
variants currently needed. Regime post-processing expects the retained
stride-10 settings; incompatible caches stop this wrapper rather than trigger
recomputation. Never use it as a substitute for validating a new fit.
SERRA/source panels remain in `Luc_i_a/fig`; composed JRSS-C panels are in
`Luc_i_a_JRSSC/fig`. Letter overlays remain in LaTeX.

To export separate submission figures with the same LaTeX panel letters:

```powershell
python Manuscript_figures/scripts/export_jrssc_figures.py --pdflatex "C:/texlive/2025/bin/windows/pdflatex.exe"
```

This author-only packaging step reads the current main/supplement figure blocks
and existing images, not fitted models. It writes standalone, caption-free PDFs
and a numbering manifest to `Luc_i_a_JRSSC/fig_exported/`; TeX sources and logs
remain in its `_build/` subfolder. Original figures and manuscript sources are
not changed. The 157 mm width and 9 bp sans-serif labels match the current
modern/medium OUP setup; change `--width-mm` if manuscript dimensions change.
Top padding prevents the conceptual figure's above-image labels being clipped.
`--dry-run` validates extraction without compiling or writing output.

## 4. Refit the four real-data specifications (expensive)

Only launch deliberately after the frozen release assets pass checksum and
input validation and are installed under Data/. This runs the matched analyses,
not a lightweight repository check.

```powershell
$env:DCS_RSCRIPT = "C:/Program Files/R/R-4.4.1/bin/Rscript.exe"
./Real_data_analysis/model_comparison/scripts/run_four_models.ps1
```

The wrapper sets the unique-location input, all four models, 30,000/20,000
iterations, no thinning, online prediction at every retained iteration and
eight OpenMP threads. Reuse is enabled; existing completed fits are retained.
It then summarises the comparison, which can be costly. Do not change the
result name or discard caches without recording that decision.

## 5. Refit simulations or the spatial benchmark (expensive)

The run/summarise scripts under `Simulations/clustering_recovery`,
`Simulations/spatial_prediction` and `Simulations/model_comparison` control
those experiments. `Simulations/gp_nngp_benchmark/scripts/benchmark_gp_nngp.R`
is the standalone spatial benchmark. Preserve scenario settings and seeds.
`run_prediction_experiment.R` supports smoke and production runs; set
`DCS_PREDICTION_SMOKE=false` and use its explicit production result setting.
Do not assume a smoke run reproduces the manuscript's 100-replication design.

## 6. Compile the manuscripts

This local author workflow requires the manuscript projects, which are not
included in the code-only Git repository.

```powershell
Push-Location Luc_i_a_JRSSC
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex supplementary.tex
Pop-Location
Push-Location Luc_i_a
latexmk -pdf -interaction=nonstopmode -halt-on-error main_v2.tex
Pop-Location
```

Final submission checks remain separate: verify referenced numbers, panel-letter
exports, data/code availability URLs, author disclosures and supplementary-file
size. Cleanup does not establish convergence or real-data predictive calibration.
