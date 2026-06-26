# AGENTS.md

Agent orientation for the `lava` R package (Structural Equation Models).

## Commands

```sh
make test          # devtools::test() — standard unit tests
make check         # rcmdcheck with --ignore-vignettes; _R_CHECK_FORCE_SUGGESTS_=0
make check-cran    # builds tarball first, then R CMD check --as-cran --run-donttest
make lint          # lintr::lint_package() — no .lintr config, uses lintr defaults
make roxygen       # regenerate docs (alias: make doc)
make coverage      # covr::report() — opens HTML in browser (macOS only)
make pkgdown       # build pkgdown site
make clean         # removes vignette artifacts, tests/lib, build/
```

Run a focused test:
```r
devtools::test(filter = "estimate_default")  # substring match on test file name
```

## Non-obvious gotchas

- **`Rgraphviz`** (Bioconductor) is a Suggests dep for path diagrams. It is commented out in CI and treated as optional (`_R_CHECK_FORCE_SUGGESTS_=0`).
- `future.apply` and `progressr` are hard **Imports** (not Suggests) — parallelization in `sim.default()` is not optional infrastructure.
- use `devtools::load_all(".")` to test the package interactively

## Documentation

- Roxygen comments use **Markdown** (`Roxygen: list(markdown = TRUE)` in DESCRIPTION).
- `README.md` is generated from `README.Rmd` — edit the `.Rmd`, not the `.md` (`make readme`).
- `NEWS.md` is the hand-maintained user-facing changelog for CRAN.

## Code style

- 2-space indent, LF line endings, UTF-8 (enforced by `.editorconfig`).
- R and C/C++ files require a final newline; general files do not.
- Single hashes for comments and roxygen documentation

## CI

- CI (`R-CMD-check.yaml`) runs on `ubuntu-latest` with R `devel` only; triggers on push/PR to `main`, `master`, `dev`.
- pkgdown deploys to `gh-pages` on successful push to `main`.
- Coverage (`test-coverage.yaml`, Codecov) runs on push/PR to `main`/`master`; covers tests + `@examples`.
- Cross-platform R-hub checks are manual (`workflow_dispatch` only, requires `RHUB_TOKEN`).
