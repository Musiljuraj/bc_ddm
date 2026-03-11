# bc_ddm — Octave implementation of FEM + FETI-DP + BDDC (+ spectral analysis)

This repository contains a **GNU Octave / MATLAB** implementation developed for a bachelor thesis on **domain decomposition methods**.

Core components live under `src/` (FEM, DDM, FETI-DP, BDDC), with tests grouped by thesis chapters in `tests/`.
The current main entry point focuses on **Chapter 6 spectral analysis + PCG diagnostics** comparing **FETI-DP vs BDDC**.

---

## Repository layout

- `src/` — library code
  - `fem/` — P1 FEM mesh + assembly routines
  - `ddm/` — subdomain decomposition + interface bookkeeping
  - `feti_dp/` — FETI-DP setup / operator / preconditioner / PCG solve
  - `bddc/` — BDDC setup / operator / preconditioner / PCG solve
  - `common/` — shared helpers/utilities
- `tests/` — tests grouped by thesis chapters (incl. Chapter 6 spectra)
- `main/` — runnable entry points
  - Chapter 6 (current): `run_ch6_spectral_experiments.m`, `ch6_run_case.m`, `ch6_define_cases.m`
  - `setup_paths.m` — adds `src/` and `tests/` to the Octave path
- `main/archive/` — older/demo runners preserved for reference (earlier chapters)

---

## Requirements

- **GNU Octave** (recommended) or MATLAB.
- No external packages are expected beyond standard sparse linear algebra.

---

## Quick start (Chapter 6: spectral experiments)

From the **repo root**:

### Option A — interactive Octave

```matlab
addpath('main');
run_ch6_spectral_experiments();
```

### Option B — command line

```bash
octave --eval "addpath('main'); run_ch6_spectral_experiments();"
```

### What it generates

The batch runner:
- saves **one `.mat` per case** to `output/mats/ch6/`
- saves **PDF figures** to `output/figures/ch6/`
  - sorted eigenvalue overlay
  - histogram overlay
  - residual history overlay
- exports a **LaTeX summary table** to `output/tables/ch6/table_ch6_summary.tex`

---

## Running a single case

```matlab
addpath('main');
setup_paths();

cfg = struct('n', 32, 'nSubX', 2, 'nSubY', 2, 'seed', 1);
out = ch6_run_case(cfg);
```

### Notes on parameters

- `n` must be divisible by `nSubX` and `nSubY` (structured partition assumption).
- The default Chapter 6 case set (Set A + Set B) is defined in `main/ch6_define_cases.m`.

---

## Spectral analysis controls (important)

Yes — this explanatory section belongs in `README.md`. It tells the reader when spectra are computed and which knobs to change.

`ch6_run_case.m` computes **full spectra** only when the method’s PCG-space dimension is small enough:

- FETI-DP: `nLambda` (multiplier-space dimension)
- BDDC: `nHat = size(R,2)` (assembled hat-space dimension)

Full spectra are computed only if `n_method <= cfg.nmax`.

Useful knobs:
- `cfg.do_spectra` (default `true`)
- `cfg.nmax` (default in `ch6_run_case`: 800; the batch runner sets it higher)

If spectra are skipped, the output contains a `spec_skipped` flag and a `spec_skip_reason`.

---

## Older runners (Chapters 2–4 demos / integration)

`main/archive/` contains earlier demo scripts and block integration runners for FEM/DDM/FETI-DP/BDDC, kept for reference.

Typical usage:

```matlab
addpath('main');
setup_paths();
addpath('main/archive');

% examples (see archive for the available scripts)
% run_fetidp_demo();
% run_bddc_demo();
% run_compare_solvers();
```
