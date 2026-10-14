# Software dependencies

The retained runs used R 4.4.1 on Windows. `SESSION_INFO.txt`, checkpoint
settings and `RUN_SETTINGS.txt` are authoritative environment records.
Cleanup does not create a dependency lockfile or prove platform equivalence.

Core packages include Rcpp, RcppArmadillo, Matrix, data.table, RANN, salso,
ggplot2, patchwork, maps and scales; raw preparation uses arrow, dplyr and ncdf4.
The original notebooks use Python airbase, pandas, pyarrow, nest_asyncio and
cdsapi. ERA5 retrieval requires a CDS account, accepted dataset terms and
credentials outside the repository. Never commit API keys or .cdsapirc. Consult each
script's package checks for optional requirements. Compile with a compatible
C++ toolchain and OpenMP (matching Rtools on Windows). Do not bypass the final
kernel's Matrix/CHOLMOD runtime check or reuse incompatible compiled libraries.

Figure composition uses Python, pypdf and Poppler's pdftoppm. The manuscripts
use latexmk, the supplied OUP class and preamble packages. `latexmkrc` handles
cross-document references. Scientific inputs are separate from package
dependencies; see `Data/README.md` and `docs/REPRODUCIBILITY.md`.

The manifest builder itself uses base R. It needs neither a historical fit nor
parent-repository caches once the candidate PM10/ERA5 RDS files are prepared.
The original notebooks are source records, not unattended production launchers;
use `Data/preparation/scripts` for portable preparation entry points.
