# DCS-sPPM: European PM10 and simulation reproducibility

This is the reproduction-code project for a Bayesian decomposition of shared spatial variation and heterogeneous local temporal dynamics. Git contains source code and reproduction documentation only. Scientific inputs, results, figures and both manuscript projects remain available locally but are excluded from Git. Cleanup does not change the statistical model or fitted results.

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
| `Data/` | Local frozen inputs; see `Data/README.md` before sharing |
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

The two other real-data kernels and `Real_data_analysis/mixing_diagnostics/scripts/benchmark_mixing.R` remain because SERRA still discusses their mixing diagnostics. They are not production alternatives. Preparation now reads compact frozen metadata from `Data/*_preparation.rds`, not historical posterior fits. The historical diagnostic warm start contains only one set of starting values; final unique-location analyses still start from the same fixed values as before.

## Safe first check

From this directory:

```powershell
python Tools/scripts/audit_repository.py
```

The audit reads paths and file metadata only. It does not fit models, recompute VI, render figures or load posterior arrays. Missing local results and manuscript projects are reported but do not fail this code-only check. On a complete local installation, add `--require-results --require-manuscripts` to require both.

## GitHub and data release

Only reproduction code, configuration/path helpers and useful reproduction documentation are versioned. Every experiment's `results/` directory is ignored, including compact CSV summaries and run metadata. Both LaTeX projects, input data, exported graphics, compiled objects, caches and local maintenance/editorial reports are also ignored; none is deleted from disk. Scientific inputs can be distributed as GitHub Release assets, separately from the code history, once redistribution permissions and source attribution have been checked. A clone requires those inputs to refit, and saved current results to regenerate graphics without refitting. See `Data/README.md` for the input inventory. Figure composition and manuscript compilation additionally require the separately supplied manuscript projects and source assets. No input release has been uploaded yet.

This cleanup does not publish, push, commit or license the project. Data redistribution terms and the authors' choice of code license must be checked before release. Manuscript availability placeholders require author action.

## Superseded material

Unused drafts, pilot/holdout outputs and previous crossover work were moved to the sibling folder `DCS_sPPM superseded experimental result`, preserving relative paths. They were not permanently deleted. The inventory is in `docs/cleanup_2026-10-13.csv`. No active workflow should use archived files without restoring and documenting that dependency.
