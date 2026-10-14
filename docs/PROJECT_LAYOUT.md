# Experiment-based layout

```
Simulations/
  clustering_recovery/{scripts,results}/
  spatial_prediction/{scripts,results}/
  model_comparison/{scripts,results}/
  gp_nngp_benchmark/{scripts,results}/
  validation/{scripts,results}/
Real_data_analysis/
  preparation/{scripts,results}/
  model_comparison/{scripts,results}/
  maps/scripts/
  mixing_diagnostics/{scripts,results}/
  reference_inputs/results/
Manuscript_figures/{scripts,results}/
Tools/scripts/
Data/preparation/{scripts,manifests}/
```

Real-data map outputs remain beside their generating final fits, under each
model's `maps` or `predictive_temporal_regime` folder. This avoids duplicating
large prediction outputs. The manuscript figure folders remain beside their
TeX sources. Shared R helpers, C++ sources and scientific inputs stay in `R`,
`src` and `Data`.

Git retains scripts, shared helpers, C++ sources and reproduction documentation.
The diagram above describes the complete local installation: all `results/`
directories, scientific inputs, exported graphics, both LaTeX folders and local
maintenance/editorial reports are ignored, not deleted. A code-only clone must
obtain raw inputs and follow Data/preparation/README.md before running analyses.
Preparation code and station/grid/scaling manifests are versioned; raw
downloads and generated input arrays are not. The final-cohort builder is
independent of the parent repository and historical posterior fits.

Ordinal script prefixes were removed. Result names containing iteration counts,
replication counts and years are scientific configuration identifiers and were
not removed. Environment-variable names and supported model identifiers remain
unchanged.

`layout_migration_2026-10-13.json` maps old paths to new paths;
`layout_moves_2026-10-13.csv` records directory moves. The earlier cleanup audit
and the migration map deliberately preserve historical source paths. Frozen RDS
objects also retain original provenance strings; their binary contents were
not rewritten. Active code and current path manifests use the new locations.

Validation consists of dependency/documentation-link checks, R/Python/PowerShell
syntax checks, a dry-run of the production launcher and unchanged checksums for
16 scientific inputs/kernels. No MCMC, VI search, result summarisation or
figure regeneration was needed for the reorganisation.
