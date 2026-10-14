# Preparation release check, 2026-10-14

The supplied preparation sources are included unchanged, alongside portable
entry points and a frozen manuscript input manifest. This is one code
repository; no separate historical-results repository is required.

## Checks completed / scope

- Versioned source paths and download/preparation instructions are present.
- Raw inputs, outputs, results, credentials and manuscript folders are ignored.
- Frozen final station decisions, scaling and grid are text/CSV configuration.
- Numeric fingerprints cover observed and imputed PM10 and matched meteorology.
- The builder refuses to create analysis objects if candidate inputs differ.
- Manifest creation reads existing inputs, not samplers or posterior fits.
- The existing PM10/weather inputs pass the checker for all 2,989 observed and
  imputed series and all matched weather cells. All 45 R sources parse, and
  the repository dependency/link audit passes.
- Raw parquet/NetCDF regeneration was not run; arrow and ncdf4 are absent in
  the validation runtime. This is explicit in the preparation README.

Passing a check on existing prepared inputs verifies manifest consistency,
not raw-download regeneration. Before the public release is described as
independently reproduced, run preparation on the downloaded raw files and
require the reconstruction reports to pass. Sources can revise historical
records, so retain acquisition dates, product request details and raw checksums.

## Manuscript provenance reconciliation

The supplied PM10 source filters Validity=1 and positive concentrations,
averages hourly records without an 18-hour eligibility threshold, filters native
daily observations below 500, then keeps the lower value when daily/hourly
products overlap. It does not explicitly filter Verification=1. The notebook
requests the Verified dataset, which should not be described as a separate
per-record filter without verifying the source contents. The final station
selection consumes these daily values; its provenance audits are not a new
aggregation or native-daily-priority operation.

The ERA5 source computes magnitude from daily-mean wind components, not a
daily average of hourly wind speeds. It retains the retrieved daily statistic
of tp, not an independently constructed daily precipitation total. The original
CDS box [72,-11,35,34] has EAST=34 and SOUTH=35; it differs from the PM10
study bounding box. The frozen weather-cell assignments reach Longitude=34
and Latitude=35. Preserve this matching when reproducing the analysis and
describe the covariate coverage/nearest-cell rule accurately.

Do not rewrite the model, observations or source rules solely to make them
agree with old manuscript wording. Resolve any differences using the raw
regeneration and fingerprint reports. No new MCMC is needed if the regenerated
input arrays agree. Coordinate selection and standardisation constants must
remain as in the final manifest.

## Public repository actions (authors)

Confirm data/metadata attribution and redistribution terms for the coordinate
manifest; choose a code licence; verify the raw-download reconstruction;
then commit and push the prepared code to the existing DCS-sPPM repository.
Provide reviewer access and replace manuscript repository placeholders.
No licence, commit, push or publication is performed by preparation scripts.
