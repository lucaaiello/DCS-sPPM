# Manuscript-to-code and result map

JRSS-C uses `Luc_i_a_JRSSC/main.tex` and `supplementary.tex`. SERRA uses
`Luc_i_a/main_v2.tex` with an embedded supplement. Labels, not changing figure
or table numbers, identify results.

| Result | Numerical source | Scripts |
| --- | --- | --- |
| Clustering and parameter recovery | `Simulations/clustering_recovery/results/clustering_recovery_100rep/*summary.csv` | [Run](../Simulations/clustering_recovery/scripts/run_clustering_recovery.R), [summarise](../Simulations/clustering_recovery/scripts/summarize_clustering_recovery.R) |
| Process/regime prediction | `Simulations/spatial_prediction/results/prediction_experiment_100rep/*summary.csv` | [Run](../Simulations/spatial_prediction/scripts/run_prediction_experiment.R), [summarise](../Simulations/spatial_prediction/scripts/summarize_prediction_experiment.R) |
| Four-model simulation comparison | `Simulations/model_comparison/results/production` | [Run](../Simulations/model_comparison/scripts/run_model_comparison.R), [summarise](../Simulations/model_comparison/scripts/summarize_model_comparison.R) |
| GP/NNGP timing table | `Simulations/gp_nngp_benchmark/results` | [Benchmark](../Simulations/gp_nngp_benchmark/scripts/benchmark_gp_nngp.R) |
| Final sample and provenance | Frozen input-release assets; historical audit files remain local | [Install](../Data/README.md), [validate](../Tools/scripts/validate_preparation_inputs.R) |
| Final application and comparators | `Real_data_analysis/model_comparison/results/production_30000_unique_locations` | [Four-model launcher](../Real_data_analysis/model_comparison/scripts/run_four_models.ps1), [summarise](../Real_data_analysis/model_comparison/scripts/summarize_full_analysis.R) |
| VI, stability and representatives | Final `dcs_sppm/maps/clustering/vi_partition.rds` and CSVs | [Render](../Real_data_analysis/maps/scripts/render_prediction_and_clustering.R); reuse VI |
| Predictive regime summaries | Final `dcs_sppm/predictive_temporal_regime` | [Render](../Real_data_analysis/maps/scripts/generate_predictive_regime_maps.R); reuse summary |
| April/December response maps | Final stored online prediction moments | [Render](../Real_data_analysis/maps/scripts/render_prediction_and_clustering.R), [compose](../Manuscript_figures/scripts/package_jrssc_figures.py) |
| Empirical maps and temporal summary | Final monitored data | [Generate](../Real_data_analysis/maps/scripts/generate_exploratory_figures.R) |
| Comparative figures | Final simulation/application CSV summaries | [Generate](../Manuscript_figures/scripts/generate_comparative_figures.R) |
| Conceptual diagram | Code-native schematic | [Generate](../Manuscript_figures/scripts/generate_conceptual_scheme.R) |
| JRSS-C composed figures | Required source panel PDFs | [Compose](../Manuscript_figures/scripts/package_jrssc_figures.py), [refresh](../Manuscript_figures/scripts/refresh_assets.ps1) |

Manuscript tables are manually typeset. This map is provenance, not a claim
that cleanup rewrote TeX tables or revalidated every numerical entry.

## Frozen preparation inputs and historical diagnostics

`Data/production_30000_online_grid_preparation.rds` supplies original station
metadata, weather matching, scaling, geographic settings and the prediction
grid to `run_full_analysis.R`. It contains
no posterior chains. Historical fits are not required to refit the
final application. The comparison summary requires the current DCS-sPPM fit
in the same result directory as its competitors and never substitutes an older fit.

Historical mixing checks use `benchmark_mixing.R`, its two real-data
kernels, `Data/production_10000_single_chain_preparation.rds` and the
`mixing_benchmark_500*` results remain for that reason. They do not establish
convergence of the final unique-location fit.

Raw data preparation and selection code was archived in the local
`DCS_sPPM superseded experimental result/frozen_input_workflow_2026-10-14`
folder. It is not a dependency of the public frozen-input workflow.

`R/pilot_simulation_utils.R` and `R/full_simulation_utils.R` are dependencies
of final simulation workflows despite their historical names.

## Interpretation and submission

- Prediction limits are pointwise moment-matched Gaussian approximations.
- Agreement with DCS-sPPM is not observed-truth accuracy.
- Certainty/entropy concern VI-projected allocations, not PM10 uncertainty.
- Panel letters are added in LaTeX. Submitted individual figures need a composed
  export, not the unlabelled source PDF alone.
- The isolated spatial benchmark does not establish approximation accuracy.
