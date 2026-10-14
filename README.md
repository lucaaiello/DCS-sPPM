# DCS-sPPM: European PM10 and simulation reproducibility

This is the reproduction-code project for a Bayesian decomposition of shared spatial variation and heterogeneous local temporal dynamics. Git contains source code, documentation and text/CSV input manifests, not raw observations, fitted results or manuscripts. The [download-to-analysis workflow](Data/preparation/README.md) rebuilds the final application inputs from PM10 and ERA5 downloads and checks them against the frozen manuscript sample before use.

## Start here

| Folder | Purpose |
| --- | --- |
| `Luc_i_a_JRSSC/` | Local-only JRSS-C manuscript and figures; ignored by Git |
| `Luc_i_a/` | Local-only SERRA manuscript and source figures; ignored by Git |
| `src/` | Production/prediction C++ kernels, spatial benchmark and regime post-processing |
| `R/` | Shared preparation, simulation, sampling, summary and plotting helpers |
| `Simulations/` | Experiment folders, each with its own `scripts/` and `results/` |
| `Real_data_analysis/` | Preparation, model comparison, maps, mixing diagnostics and reference inputs |
| `Manuscript_figures/` | Shared schematic/comparative/composition scripts and supporting assets |
| `Tools/` | Read-only path and source-syntax audits |
| `Data/` | Local inputs; `Data/preparation/` contains versioned download/preparation code and manuscript manifests |
| `docs/` | Reproduction instructions, dependencies, manuscript-to-code map and layout migration; local audit reports/editorial notes are ignored |

See [the manuscript map](docs/MANUSCRIPT_MAP.md), [reproduction instructions](docs/REPRODUCIBILITY.md) and [dependencies](docs/DEPENDENCIES.md). The [layout migration map](docs/layout_migration_2026-10-13.json) records former numbered script names and result locations. Script prefixes were removed; meaningful years, iteration counts and replication counts in scientific identifiers were preserved.

## Frozen manuscript analysis

The application uses **2,989 unique monitored locations**, **366 days in 2024** and **5,297 prediction locations**. All four full-data fits use 30,000 MCMC iterations, 20,000 burn-in iterations and 10,000 retained iterations. Final results are under:

```
Real_data_analysis/model_comparison/results/production_30000_unique_locations/
```

The specifications are `dcs_sppm`, `ppm_common_shock`, `k1_common_shock` and `sppm_no_common_shock`. These are matched structural comparisons, not held-out real-data accuracy scores. Cached VI and predictive-regime summaries are retained and should be reused.

Simulation sources are `Simulations/clustering_recovery/results/clustering_recovery_100rep` and `Simulations/spatial_prediction/results/prediction_experiment_100rep`. Competitor checkpoints and four-model summaries are in `Simulations/model_comparison/results/production`. The isolated spatial benchmark is in `Simulations/gp_nngp_benchmark/results`.

## Production kernels

- `src/sampler_global_beta_single_w.cpp`: clustering recovery.
- `src/sampler_global_beta_single_w_prediction.cpp`: online simulation prediction.
- `src/sampler_global_beta_real_data_online_prediction.cpp`: final real-data fits.
- `src/predictive_regime_postprocess.cpp`: cached-partition predictive summaries.
- `src/benchmark_gp_nngp.cpp`: isolated spatial computation.

The other real-data kernels and `Real_data_analysis/mixing_diagnostics/scripts/benchmark_mixing.R` remain as historical diagnostic code, not production alternatives. Final analyses start from fixed values and do not require historical posterior fits. A fresh clone recreates the compact preparation metadata from versioned station/grid/scaling manifests; the historical diagnostic warm start is optional.

## Safe first check

From this directory:

```powershell
python Tools/scripts/audit_repository.py
```

The audit reads paths and file metadata only. It does not fit models, recompute VI, render figures or load posterior arrays. Missing local results and manuscript projects are reported but do not fail this code-only check. On a complete local installation, add `--require-results --require-manuscripts` to require both.

## GitHub and data release

Code, preparation notebooks/scripts, input manifests and useful documentation are versioned. All results, raw downloads, generated RDS inputs, graphics, caches and both LaTeX projects remain ignored. Obtain PM10 and ERA5 from their original providers and follow [the preparation instructions](Data/preparation/README.md). The final-cohort reconstruction does not require a separate historical-results repository. Saved results are needed only to regenerate figures without refitting. Manuscript projects and source images are local author-only resources.

The supplied preparation differs from current manuscript wording on daily/hourly precedence, verification/hour-count filtering and meteorological definitions. Portable entry points preserve the supplied source rules and require agreement with frozen numeric fingerprints instead of silently imposing different scientific choices. Raw-download regeneration must pass that check before the workflow is described as independently validated. No fitting is needed for this check.

This cleanup does not publish, push, commit or license the project. Data redistribution terms and the authors' choice of code license must be checked before release. Manuscript availability placeholders require author action.

## Superseded material

Unused drafts, pilot/holdout outputs and previous crossover work were moved to the sibling folder `DCS_sPPM superseded experimental result`, preserving relative paths. They were not permanently deleted. The inventory is in `docs/cleanup_2026-10-13.csv`. No active workflow should use archived files without restoring and documenting that dependency.
