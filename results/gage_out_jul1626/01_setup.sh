#!/usr/bin/env bash
# Set up R + gage and the Python side. Idempotent: skips what's already present.
set -euo pipefail

echo "== 1. R + gage =="
if ! command -v Rscript >/dev/null 2>&1; then
  # Debian/Ubuntu; on other distros install R your usual way.
  sudo apt-get update -qq && sudo apt-get install -y -qq r-base-core
fi
Rscript --version

# gage's CORE path (gagePrep -> gs.tTest -> gageSum) needs only base R + stats.
# Two ways to get gage; the script auto-detects either at run time:
#   (a) the real Bioconductor package (best if you have internet to Bioconductor)
#   (b) a clone of the source repo (what we use when Bioconductor is unreachable)
if ! Rscript -e 'quit(status = as.integer(!requireNamespace("gage", quietly=TRUE)))' 2>/dev/null; then
  echo "gage package not installed; trying BiocManager..."
  Rscript -e 'if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install("gage", update=FALSE, ask=FALSE)' \
    || { echo "Bioconductor unreachable -> cloning source instead"; \
         [ -d gage_R ] || git clone --depth 1 https://github.com/datapplab/gage.git gage_R; }
fi

echo "== 2. Python (pygage + deps) =="
# From the pygage repo root (adjust path if needed):
#   pip install -e .            # installs pygage + deps
# Or just the deps used by the comparison:
python3 -m pip install --quiet --break-system-packages polars numpy scipy matplotlib pandas pyarrow || true
python3 -c "import polars, numpy, scipy, matplotlib; print('python deps OK')"
echo "setup done."
