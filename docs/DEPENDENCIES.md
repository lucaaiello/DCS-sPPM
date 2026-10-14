# Software dependencies

The retained runs used R 4.4.1 on Windows. `SESSION_INFO.txt`, checkpoint
settings and `RUN_SETTINGS.txt` are authoritative environment records.
Cleanup does not create a dependency lockfile or prove platform equivalence.

Core packages include Rcpp, RcppArmadillo, Matrix, data.table, RANN, salso,
ggplot2, patchwork, maps and scales. Frozen-input reproduction does not require
raw-download clients, CDS credentials, arrow or ncdf4. Consult each
script's package checks for optional requirements. Compile with a compatible
C++ toolchain and OpenMP (matching Rtools on Windows). Do not bypass the final
kernel's Matrix/CHOLMOD runtime check or reuse incompatible compiled libraries.

Figure composition uses Python, pypdf and Poppler's pdftoppm. The manuscripts
use latexmk, the supplied OUP class and preamble packages. `latexmkrc` handles
cross-document references. Scientific inputs are separate from package
dependencies; see `Data/README.md` and `docs/REPRODUCIBILITY.md`.

Input-release packaging and SHA-256 verification use the Python standard
library. The input checker uses the retained R helpers to validate alignment,
scaling and the prediction design without compiling or launching a sampler.
Archived raw-preparation tools are not dependencies of production reproduction.
