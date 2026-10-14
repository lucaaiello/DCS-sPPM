# Frozen-input release workflow

Manuscript release supplied by the author:
[v1.0.0](https://github.com/lucaaiello/DCS-sPPM/releases/tag/v1.0.0).
Use this versioned page for the three RDS files, input README and checksums.
Confirm all five assets are visible to an authorised reader before submission;
the recorded link alone is not an independent upload/access verification.

## Scope

This repository reproduces reported analyses from the exact prepared inputs
used in the manuscript. It does not rebuild those inputs from fresh raw
downloads. The final data-selection, imputation and meteorological procedures
must still be documented accurately; distributing inputs does not remove
provenance obligations. In particular, the authors must reconcile source
filtering/aggregation descriptions before submission. No new MCMC is needed
merely to package these inputs, and no original RDS is rewritten.

Raw download notebooks, raw-preparation scripts, selection audits' source
code and newly developed raw reconstruction/manifests were archived locally
under `DCS_sPPM superseded experimental result/frozen_input_workflow_2026-10-14`.
Their historical audit results remain local. They are not production dependencies.

## Author packaging (no MCMC)

From the reproduction project root:

```powershell
python Tools/scripts/package_frozen_inputs.py --release-tag v1.0.0
python Tools/scripts/verify_frozen_inputs.py `
  --input-dir Tools/results/frozen_inputs_release `
  --checksums Tools/results/frozen_inputs_release/SHA256SUMS
Rscript Tools/scripts/validate_preparation_inputs.R
```

The output directory contains the three required RDS files, INPUTS_README.md
and SHA256SUMS. It is ignored by Git. Files are copied byte-for-byte and checked
against the originals. Existing output directories are not overwritten;
choose another --output when preparing a later release. --dry-run reports
paths and sizes without writing. The code checkout itself need not contain
data when uploaded; package from the author's local input folder using --source.

## Upload to the same GitHub repository

After confirming redistribution permissions/attribution and reviewing the
input README, commit the reproduction code. On the DCS-sPPM repository's
Releases page create the matching tag/release (e.g. v1.0.0) and attach all five
files. Do not git-add the RDS files or the ignored packaging directory. GitHub
release assets are separate from ordinary Git file storage; each asset must
be below 2 GiB. See [GitHub release documentation](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases#storage-and-bandwidth-quotas).

The packaging script does not grant redistribution permission or upload
anything. Private releases require repository read access. Arrange reviewer
access and public availability according to the journal's policy. Do not
invent a release URL before uploading; use the actual tag's URL in the
manuscript/Data Availability Statement, not an unversioned latest link.

## Reviewer installation

Download the matching assets and SHA256SUMS into a temporary directory and
verify them using Tools/scripts/verify_frozen_inputs.py before installing.
Copy the three RDS files to Data/ without overwriting existing inputs. Then:

```powershell
python Tools/scripts/audit_repository.py --require-inputs
Rscript Tools/scripts/validate_preparation_inputs.R
```

The checks validate access, matrix alignment, finite covariates, station/date
counts, grid and starting values. They do not fit a model or assess convergence.
Only after validation should the reader use the expensive run commands in
[reproduction instructions](REPRODUCIBILITY.md).
