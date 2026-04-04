# Reproducibility Guide (JOSS)

This document describes how to reproduce the results and artifacts used in the
JOSS paper for `FastTanhSinhQuadrature.jl`.

## Scope

This guide reproduces:

- package tests (including JET type-stability checks on supported Julia versions),
- benchmark result tables and CSV files,
- paper benchmark figure (`benchmark_summary.svg`),
- paper convergence figure (`convergence.svg`),
- docs build containing the same benchmark/convergence assets.

No external datasets are required. All inputs are defined in repository scripts.

## 1. Reproducible Snapshot (Tag/Commit)

Always reproduce from a fixed snapshot:

```bash
git clone https://github.com/svretina/FastTanhSinhQuadrature.jl.git
cd FastTanhSinhQuadrature.jl
git checkout <tag-or-commit>
```

For publication, use the archived release tag associated with the Zenodo DOI.

## 2. Software and System Requirements

- OS: Linux/macOS (CI is tested on `ubuntu-latest` x86_64).
- Julia:
  - `1.11` and `1.12` for full test reproducibility (includes JET checks),
  - `1.13` is used in this repo for docs/figure generation workflows.
- Optional (paper PDF build): Docker (or GitHub Actions artifact download).

Notes:

- Benchmarks are CPU/OS dependent in absolute timing; compare trends/orderings.
- On Julia `1.11`/`1.12`, tests install `JET` dynamically during `Pkg.test()`.

## 3. Dependency Specification

Dependency environments are pinned via `Project.toml` + `Manifest.toml`:

- package/runtime: `./Project.toml`, `./Manifest.toml`
- benchmark workflow: `./benchmark/Project.toml`, `./benchmark/Manifest.toml`
- docs/figure workflow: `./docs/Project.toml`, `./docs/Manifest.toml`
- tests: `./test/Project.toml` (with root project)

Instantiate all environments:

```bash
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=benchmark -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=docs -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

## 4. Exact Reproduction Commands

### 4.1 Run package tests (including JET on 1.11/1.12)

```bash
julia +1.11 --project -e 'using Pkg; Pkg.test()'
julia +1.12 --project -e 'using Pkg; Pkg.test()'
```

### 4.2 Regenerate benchmark tables/results used in paper/docs

```bash
julia --project=benchmark benchmark/benchmarks.jl
```

Outputs:

- `benchmark/results/timings.csv`
- `benchmark/results/timings_full.md`
- `benchmark/results/timings_summary.md`

### 4.3 Regenerate paper/docs benchmark summary figure

```bash
julia +1.13 --project=docs JOSS_paper/benchmark_summary_plot.jl
```

Outputs:

- `JOSS_paper/benchmark_summary.svg`
- `JOSS_paper/benchmark_summary.png`
- `docs/src/assets/benchmark_summary.svg`
- `docs/src/assets/benchmark_summary.png`
- `benchmark/results/benchmark_summary.svg`
- `benchmark/results/benchmark_summary.png`

### 4.4 Regenerate paper/docs convergence figure

```bash
julia +1.13 --project=docs docs/convergence_plots.jl
```

Outputs:

- `JOSS_paper/convergence.svg`
- `docs/src/assets/convergence.svg`
- `docs/src/assets/convergence.png`

### 4.5 Build documentation website

```bash
julia +1.13 --project=docs -e 'include("docs/make.jl")'
```

Primary output:

- `docs/build/`

### 4.6 Build JOSS draft PDF

Preferred reproducible route (matches CI in this repository):

1. Trigger workflow `.github/workflows/joss_paper.yml` on your snapshot.
2. Download the `paper` artifact produced by the workflow.

If using GitHub CLI:

```bash
gh workflow run joss_paper.yml --ref <tag-or-commit>
gh run list --workflow joss_paper.yml --limit 1
gh run download <run-id> --name paper --dir JOSS_paper/artifact
```

Expected output file:

- `JOSS_paper/artifact/paper.pdf`

Optional local Docker route (Open Journals `inara` image):

```bash
docker run --rm --platform linux/amd64 \
  -v "$PWD/JOSS_paper:/data" -w /data \
  openjournals/inara:latest paper.md -p -o pdf
```

This should write `JOSS_paper/paper.pdf`.

## 5. Mapping: Paper Artifact -> Script -> Output

- JOSS Figure 1 (benchmark speedup plot):
  - scripts: `benchmark/benchmarks.jl` then `JOSS_paper/benchmark_summary_plot.jl`
  - outputs: `benchmark/results/timings.csv` and `JOSS_paper/benchmark_summary.svg`
- JOSS convergence figure:
  - script: `docs/convergence_plots.jl`
  - output: `JOSS_paper/convergence.svg`
- Benchmark timing tables referenced in docs/discussion:
  - script: `benchmark/benchmarks.jl`
  - outputs: `benchmark/results/timings_summary.md`, `benchmark/results/timings_full.md`

## 6. Runtime and Platform Notes

- `Pkg.test()` on Julia `1.11`/`1.12` includes Aqua + JET and may take minutes.
- Benchmark execution time depends strongly on CPU and load.
- Absolute benchmark timings are not expected to match exactly across machines.
- Some benchmark methods may miss tolerance on specific functions; this is recorded
  in `status` columns and `*` markers in summary tables.

## 7. Known Reproducibility Limitations / Potential Blockers

- Paper PDF generation is tied to Open Journals tooling; local Docker command may
  vary by Docker/platform setup.
- Benchmark numbers are hardware-dependent and include stochastic behavior in some
  external algorithms; use relative trends and tolerance status for comparisons.
- The repository currently tracks generated artifacts (`benchmark/results/*`,
  `JOSS_paper/*.svg/png`), so users should regenerate them from scripts rather than
  treating checked-in binaries as authoritative.

## 8. Archival Release Checklist (Manual)

1. Ensure working tree is clean and CI passes on target commit.
2. Update version if needed in `Project.toml`.
3. Create a Git tag for the exact snapshot used by the paper:
   - `git tag -a vX.Y.Z -m "JOSS reproducibility snapshot"`
   - `git push origin vX.Y.Z`
4. Create a GitHub Release from that tag.
5. Ensure Zenodo-GitHub integration is enabled for this repository.
6. Confirm Zenodo archives the tagged release and mints DOI.
7. Copy DOI into:
   - JOSS review thread (required),
   - JOSS paper references,
   - README (recommended),
   - release notes (recommended).
8. Verify the commands in this file run from the release tag.
