param([string]$RscriptPath = "Rscript", [string]$PythonPath = "python")

$ErrorActionPreference = "Stop"
$project = $PSScriptRoot
while (-not (Test-Path -LiteralPath (Join-Path $project "R/project_paths.R"))) {
    $parent = Split-Path $project -Parent
    if (-not $parent -or $parent -eq $project) { throw "Could not locate the manuscript project." }
    $project = $parent
}
$resultDir = Join-Path $project "Real_data_analysis/model_comparison/results/production_30000_unique_locations/dcs_sppm"
foreach ($relative in @("pilot_fit.rds", "maps/clustering/vi_partition.rds", "predictive_temporal_regime/predictive_temporal_regime_summary.rds")) {
    if (-not (Test-Path -LiteralPath (Join-Path $resultDir $relative))) {
        throw "Missing cached manuscript input: $relative. This wrapper will not fit models or search for VI."
    }
}
$settings = @{
    DCS_GRID_RESULT_DIR = $resultDir
    DCS_REAL_FULL_RESULTS_NAME = "production_30000_unique_locations"
    DCS_RENDER_PREDICTIVE_MAPS = "false"
    DCS_REUSE_VI = "true"
    DCS_REUSE_REGIME_SUMMARY = "true"
    DCS_REQUIRE_CACHED_MANUSCRIPT = "true"
    DCS_REGIME_DRAW_STRIDE = "10"
    DCS_COMPARISON_RESULTS_NAME = "production"
    DCS_MANUSCRIPT_FIGURE_DIR = (Join-Path $project "Luc_i_a_JRSSC/fig")
    DCS_MANUSCRIPT_ASSET_DIR = (Join-Path $project "Manuscript_figures/results/jrssc")
}
$previous = @{}
foreach ($name in $settings.Keys) {
    $previous[$name] = [Environment]::GetEnvironmentVariable($name, "Process")
    [Environment]::SetEnvironmentVariable($name, $settings[$name], "Process")
}
Push-Location -LiteralPath $project
try {
    foreach ($script in @("Real_data_analysis/maps/scripts/render_prediction_and_clustering.R", "Manuscript_figures/scripts/generate_conceptual_scheme.R", "Real_data_analysis/maps/scripts/generate_exploratory_figures.R", "Real_data_analysis/maps/scripts/generate_predictive_regime_maps.R", "Manuscript_figures/scripts/generate_comparative_figures.R")) {
        & $RscriptPath (Join-Path $project $script)
        if ($LASTEXITCODE -ne 0) { throw "Figure refresh failed: $script" }
    }
    & $PythonPath (Join-Path $PSScriptRoot "package_jrssc_figures.py")
    if ($LASTEXITCODE -ne 0) { throw "JRSS-C figure composition failed." }
} finally {
    Pop-Location
    foreach ($name in $settings.Keys) {
        [Environment]::SetEnvironmentVariable($name, $previous[$name], "Process")
    }
}
