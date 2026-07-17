# Optional Rust/PyO3 kernel for the per-set statistic loop

**Status: design note, not yet implemented — deliberately gated on the regression test.**

## When it is worth it
The Python engine is fast for KEGG/GO-scale collections (hundreds of sets). The
inner per-gene-set loop in `GAGEAnalysis.run_gage` becomes the bottleneck only
for very large collections at many samples — e.g. **MSigDB C2 (~6,500 sets) ×
dozens of samples**, or permutation nulls with thousands of iterations. That is
the one clean target for a native kernel, matching the RAW Lab pure-Rust pattern
(MetaCerberus, MerCat2, NFixDB, DeGenPrime, SABER, Dagda, Merlin).

## What to port
Only the numeric hot path — everything else stays in Python:

- background moments per column (`mu`, `s`) — already O(genes × samples)
- for each set: gather rows, compute per-column `stat`, `p_up`, `p_down`
  (t/z/KS), then the Stouffer/Fisher meta-combine

Signature of the kernel (via PyO3, numpy arrays in/out):

```rust
// fn gage_kernel(
//     bg: ArrayView2<f64>,        // genes x samples (fold changes)
//     mu: ArrayView1<f64>, s: ArrayView1<f64>,
//     set_index: &[usize], n_present: usize,
//     test: Test, meta: Meta,
// ) -> (f64 /*stat_mean*/, f64 /*p_geo*/, f64 /*p_up*/, f64 /*p_down*/)
```

Use `ndarray` + `statrs` (Student-t, Normal, Gamma CDFs) to reproduce
`scipy.stats` exactly; parallelise over sets with `rayon`. Expose through
`maturin` as `pygage._kernel`, and select it at runtime:

```python
try:
    from pygage import _kernel          # Rust
    _HAVE_KERNEL = True
except ImportError:
    _HAVE_KERNEL = False                 # pure-Python fallback (this repo)
```

## Correctness gate (the important part)
**Do not** merge the kernel until it passes `tests/test_regression_gage.py`
byte-for-byte against the shipped gage-R fixtures (tolerance 1e-8). The pure
-Python path is the reference oracle; the Rust path must match it *and* gage R.
Add a CI matrix leg that builds the kernel and runs the same regression suite,
so any numeric drift in `statrs` vs `scipy` is caught immediately.

## Expected speedup
Set-parallel `rayon` over the loop typically gives a 5–20× wall-clock
improvement on C2-scale collections; permutation nulls scale ~linearly with
cores. Until a profile on real data shows the loop dominating, the Python
engine is the right default.
