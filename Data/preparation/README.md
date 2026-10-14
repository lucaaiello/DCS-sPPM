# Raw downloads to manuscript inputs

This directory contains the authors' supplied preparation sources, portable
entry points and a frozen, human-readable manuscript manifest. Raw downloads
and generated RDS files remain local and are not committed. No credentials,
posterior arrays or PM10 observations are included in the manifest.

## Status and scientific rules

The original `analysis_pm10_2024.R`, `analysis_ERA5_2024.R`, `PM10data.ipynb`
and `era5data.ipynb` preserve the supplied code as source records. Notebook
execution counts and saved console/error outputs were cleared for publication.
They are not
one-command production entry points: paths are relative, some exports are
commented out, and ERA5 quarter selection is manual. Use the scripts below.

`prepare_pm10_2024.R` preserves the supplied PM10 rules, rather than changing
the analysed observations: values must be positive and Validity=1; native
daily values additionally must be below 500; hourly values are averaged by UTC
day; where the products overlap, the lower value is retained. It does NOT
silently add Verification=1, an 18-hour threshold, or native-daily priority.
An original downloaded Verified dataset is not equivalent to an explicit
per-record Verification filter. The manuscript description of these rules
must be reconciled with the successful reconstruction check, not inferred
from the separate station-selection/provenance audit.

The audited final cohort is frozen in `manifests/selected_stations.csv`.
Its coordinates incorporate the approved location decisions. This includes
station and ERA5-cell identifiers, selection reasons and fingerprints, not
response values. The scaling constants are preserved in
`manifests/analysis_settings.R`; they are NOT recalculated on the final
2,989-location cohort. The land-grid coordinates are frozen in
`manifests/prediction_grid.csv`, avoiding dependence on a new country-polygon
release. `weather_fingerprints.csv` checks the meteorology actually used.

Reconstruction is accepted only when observed and imputed PM10 columns and
all matched weather time series agree with these manifests. Fingerprints
normalise numbers to 12 significant digits before hashing; this is a numerical
reproduction check, not a claim of byte-identical RDS serialization. Both the
observed missingness pattern and all 366 daily positions are checked. A newly
revised source download may legitimately differ and must not be advertised as
the original analysis without resolving the differences.

## 1. Retrieve PM10 and metadata

Official source: [EEA Air Quality download service](https://eeadmz1-downloads-webapp.azurewebsites.net/).
Request PM10 (pollutant 5), Verified data (dataset 2), all relevant countries,
daily and hourly products; the preparation scripts restrict observations to
2024. Save the download date and request details. The supplied
`PM10data.ipynb` documents the `airbase` Python client's retrieval commands.
Its final CSV export contains the native daily records. Retain the metadata
CSV with `Sampling.Point.Id`, `Longitude` and `Latitude` columns. Do not
commit raw files or notebook outputs. The notebook may download more years
than needed and its pandas conversion can require substantial memory; run
download/preparation deliberately, not as a smoke test.

Expected local layout (or supply alternative paths):

```
Data/preparation/raw/pm10_hourly/**/*.parquet
Data/preparation/raw/pm10_europe_daily_2024.csv
Data/preparation/raw/metadata.csv
```

## 2. Retrieve ERA5

Official source: [ERA5 daily statistics](https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics).
Register with the Climate Data Store, accept the dataset terms, and configure
`cdsapi` credentials outside the repository. The supplied `era5data.ipynb`
records the request for daily-mean t2m, tp, u10 and v10, UTC+00:00, sampled
hourly. Retrieve all four quarters of 2024 and unzip the NetCDF files into:

```
Data/preparation/raw/ERA5_2024/01_02_03_2024_euro/
Data/preparation/raw/ERA5_2024/04_05_06_2024_euro/
Data/preparation/raw/ERA5_2024/07_08_09_2024_euro/
Data/preparation/raw/ERA5_2024/10_11_12_2024_euro/
```

There must be exactly one array per variable per quarter. The original request
uses CDS area [72,-11,35,34] in NORTH,WEST,SOUTH,EAST order. This is not the
analysis bounding box [-11,46] x [34,72]. Do not enlarge the download and
silently rematch stations: the fixed station-to-weather-cell assignments and
weather fingerprints define the input used in the fitted analysis. Check
that a download covers every cell in the manifest. Regional source coverage
and nearest-cell matching need to be described accurately in the manuscript.

`prepare_era5_2024.R` decodes NetCDF time units and dimension order, converts
Kelvin to Celsius, retains the supplied tp statistic in metres, and computes
wind as sqrt(mean(u10)^2 + mean(v10)^2). This is NOT the daily mean of hourly
wind-speed magnitudes; tp is NOT converted to a daily precipitation total.
Preserve these definitions when reproducing the fitted covariates.

## 3. Prepare and validate (no MCMC)

From the repository root, using Rscript on PATH:

```powershell
Rscript Data/preparation/scripts/prepare_pm10_2024.R `
  --hourly Data/preparation/raw/pm10_hourly `
  --daily Data/preparation/raw/pm10_europe_daily_2024.csv `
  --metadata Data/preparation/raw/metadata.csv

Rscript Data/preparation/scripts/prepare_era5_2024.R `
  --input Data/preparation/raw/ERA5_2024

Rscript Data/preparation/scripts/build_manuscript_inputs.R `
  --pm10 Data/preparation/results/pm10_combined_2024.rds `
  --weather Data/preparation/results/era5_combined_2024.rds
```

Preparation requires R packages arrow, dplyr, data.table and ncdf4. Existing
downloads require no CDS/EEA client at this stage. Every script refuses to
overwrite its RDS outputs. On disagreement, the builder writes per-station
and per-cell reports and stops before creating analysis inputs. Add
`--check-only true` to validate without writing the reconstructed RDS files.

The builder applies the final station manifest and linear interpolation with
nearest-value endpoint extension. It recreates all three required production
inputs under `Data/preparation/results/reconstructed/`: the unique-location
PM10 object, compact preparation metadata and weather object. It does not
read historical fits or access the parent ST_clustering repository.

Only after successful validation, and on a fresh clone with no existing input
files, install them in the locations expected by the analysis scripts:

```powershell
foreach ($name in @('pm10_unique_locations_2024.rds',
                    'production_30000_online_grid_preparation.rds',
                    'era5_combined_2024.rds')) {
    $destination = Join-Path 'Data' $name
    if (Test-Path $destination) { throw "Input already exists: $destination" }
}
Copy-Item Data/preparation/results/reconstructed/*.rds Data/
Rscript Tools/scripts/validate_preparation_inputs.R
```

Refitting is expensive and must be launched separately; see
[reproduction instructions](../../docs/REPRODUCIBILITY.md).

## Author-only manifest maintenance

The manuscript-manifest check passed using the existing PM10 and ERA5 prepared
files: all observed/imputed station fingerprints and matched-weather
fingerprints agree. Syntax and path checks also pass. A new download-to-RDS
run has not been executed; arrow and ncdf4 were unavailable in the validation
runtime. This is a validated prepared-input bridge, not yet independent
certification of regeneration from a new raw download.

`export_manuscript_manifest.R` was used once to capture the existing final
input selection, grid, scaling and numerical fingerprints. Reviewers do not
run it: the text/CSV manifests are versioned. It reads prepared inputs only,
not MCMC results, refuses to replace existing manifests, and does not change
fitted results. A different scientific analysis requires a separate manifest
and documented settings, not replacing the manuscript reference silently.

The older `Real_data_analysis/preparation/scripts` selection/audit routines
are retained as provenance; they require historical EEA caches and metadata.
They are NOT required by the frozen-cohort reconstruction route above. The
historical mixing-benchmark warm start is also not a production dependency.
