$ErrorActionPreference = "Stop"
$project = $PSScriptRoot
while (-not (Test-Path -LiteralPath (Join-Path $project "R/project_paths.R"))) {
    $parent = Split-Path $project -Parent
    if (-not $parent -or $parent -eq $project) { throw "Could not locate the manuscript project." }
    $project = $parent
}
Set-Location -LiteralPath $project
$rscript = if ($env:DCS_RSCRIPT) { $env:DCS_RSCRIPT } else { "Rscript" }
$scripts = $PSScriptRoot
$env:DCS_REAL_FULL_DATA_FILE = Join-Path $project "Data\pm10_unique_locations_2024.rds"
$env:DCS_REAL_FULL_RESULTS_NAME = "production_30000_unique_locations"
$env:DCS_REAL_FULL_MODELS = "dcs_sppm,ppm_common_shock,k1_common_shock,sppm_no_common_shock"
$env:DCS_REAL_FULL_RUNS = "30000"
$env:DCS_REAL_FULL_BURN = "20000"
$env:DCS_REAL_FULL_THIN = "1"
$env:DCS_REAL_FULL_PREDICTION_THIN = "1"
$env:DCS_REAL_FULL_THREADS = "8"
$env:DCS_REAL_FULL_NEIGHBORS = "10"
$env:DCS_REAL_FULL_REUSE = "true"
$env:DCS_REAL_FULL_RENDER = "false"
$env:DCS_REAL_FULL_RENDER_DAILY_MAPS = "false"
$env:DCS_REAL_FULL_SUMMARIZE = "true"
Remove-Item Env:DCS_REAL_FULL_DRY_RUN -ErrorAction SilentlyContinue
Remove-Item Env:DCS_REAL_FULL_PREFLIGHT -ErrorAction SilentlyContinue
& $rscript (Join-Path $scripts "run_production_analysis.R")
if ($LASTEXITCODE -ne 0) { throw "The four-model workflow failed; inspect the preceding error." }
