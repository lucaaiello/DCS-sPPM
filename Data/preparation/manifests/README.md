# Frozen manuscript manifest

These files record the final analysed station decisions, fixed scaling,
prediction grid and input fingerprints. They contain no observed pollutant
values, meteorological time-series values, sampler chains or predictions.
Station identifiers/coordinates derive from EEA metadata and the documented
author selection. Grid coordinates derive from the original land mask.

The manifests allow a fresh clone to rebuild the final input objects without
historical posterior fits or the parent repository. Reconstruction uses
`../scripts/build_manuscript_inputs.R`. Numeric fingerprints compare values
rounded to 12 significant digits; source revisions trigger a mismatch report.
Preserve the source attribution and check metadata redistribution terms before
public release. The original PM10 and ERA5 downloads remain external inputs.
