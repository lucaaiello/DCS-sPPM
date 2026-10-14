# Frozen inputs and access

Start with [download-to-analysis preparation](preparation/README.md). Preparation
code and text/CSV manifests are versioned; raw observations, weather arrays and
generated RDS files are not. The frozen-cohort builder recreates the three
production inputs without historical fits or parent-repository caches, checking
numeric agreement before writing outputs. Source-download regeneration still
needs to pass those checks; having scripts alone does not establish validation.

The following local objects are intentionally excluded from ordinary Git:

- `pm10_combined_2024.rds`: daily observations before unique-location selection.
- `era5_combined_2024.rds`: cached weather data; too large for an ordinary GitHub file.
- `pm10_unique_locations_2024.rds`: audited final selection and imputed responses.
- `production_30000_online_grid_preparation.rds`: original station inventory,
  weather-cell matching, dates, scaling, geographic bounds and the fixed prediction
  grid. No posterior chains or predicted surfaces are included.
- `production_10000_single_chain_preparation.rds`: corresponding metadata and a
  single frozen set of initial values for historical mixing diagnostics only.

Refitting the final four models needs the unique-location object, weather cache
and 30,000-run preparation metadata. The earlier PM10 object, selection overrides,
EEA metadata and country-level observation caches are additionally required to
reproduce selection from the earlier inventory. Historical posterior fits are
no longer required by preparation or the mixing driver.

## GitHub distribution

Inputs may be attached to a versioned GitHub Release, subject to redistribution
permissions, without putting data or results into the code history. The complete
weather cache is approximately 163 MiB and exceeds GitHub's ordinary 100 MiB
per-file Git limit. Release assets or Git LFS are appropriate for that file.
The unique-location object is approximately 12 MiB; the two preparation objects
are each below 0.1 MiB. No release URL or permission to redistribute is assumed.
Publish source attribution, applicable data terms and checksums alongside inputs.

`Tools/scripts/export_preparation_inputs.R` is a one-time author migration utility:
it extracts metadata from existing reference fits, checks exact correspondence,
and refuses to overwrite frozen inputs. Reviewers receive its output and do not
need to run it. `Tools/scripts/validate_preparation_inputs.R` checks the final
data/design/grid construction without compiling or running a sampler.

Selection records and overrides remain in
`Real_data_analysis/preparation/results/unique_location_preparation`. The scripts in
`Real_data_analysis/preparation/scripts` use
existing EEA metadata and observation caches in the parent repository's
`Data/eea_pm10` tree. These inputs are outside this subproject and are not
supplied by Git. Do not rerun selection merely to regenerate figures.

A portable data release must supply frozen inputs and required source metadata
through an appropriate archive, subject to source redistribution terms, with
access dates, product names, selection overrides and checksums. The older
selection/audit route does not automatically recreate those inputs. The new
portable preparation route uses frozen station decisions and fingerprints to
rebuild the production inputs. No data-release URL is assumed here.
The final sample is the unique-location object, not an earlier inventory.
