PKG ?= lava
R = R --silent --no-save --no-echo
BUILD_DIR = build
GETVER = $(shell cat DESCRIPTION | grep Version | cut -d":" -f2 | tr -d '[:blank:]')
make_build_dir = rm -Rf $(BUILD_DIR) && mkdir -p $(BUILD_DIR)
CLIFF_CFG = .cliff.toml
CHANGELOG = CHANGELOG.md

default: check

# generate changelog for unreleased commits
cliff-unreleased:
	@git cliff --unreleased -c $(CLIFF_CFG)

cliff-unreleased-prepend:
	@git cliff --unreleased -c $(CLIFF_CFG) -p $(CHANGELOG)

roxygen:
	@echo 'devtools::document(".")' | $(R)

doc: roxygen

clean:
	@find vignettes inst '(' -name "*.html" -o -name "*_cache" -o -name "*_files" ')' -exec rm -Rf {} +
	@# remove output files of code coverage
	@rm -rf tests/lib tests/coverage-report.html
	@rm -rf build

.PHONY: build
build:
	@$(make_build_dir)
	@echo 'pkgbuild::build(path=".", dest_path="$(BUILD_DIR)", args="--compact-vignettes=qpdf --resave-data=best")' | $(R)

install:
	@echo 'remotes::install_local(".", upgrade = "never", force = TRUE)' | $(R)

dependencies-install:
	@echo 'remotes::install_deps(".", dependencies = TRUE)' | $(R)

dependencies-upgrade:
	@echo 'remotes::install_local(".", upgrade = "always")' | $(R)

check-cran: build
	@$(R) CMD check $(BUILD_DIR)/$(PKG)_$(GETVER).tar.gz --timings --as-cran --no-multiarch --run-donttest

check:
	@_R_CHECK_FORCE_SUGGESTS_=0 echo 'res <- rcmdcheck::rcmdcheck(".", build_args=c("--no-build-vignettes"), args=c("--ignore-vignettes"))' | $(R)

lint:
	@echo 'lintr::lint_package(show_progress = TRUE)' | $(R)

.PHONY: pkgdown
pkgdown:
	@$(R) -q -e "pkgdown::build_site(install=FALSE)"

vignette:
	@$(R) -q -e "devtools::build_vignettes(clean=FALSE, install=FALSE, quiet=FALSE)"

test: test-installed
test-installed: # tests locally installed version package
	@echo 'testthat::test_package("$(PKG)")' | $(R)

test-loadall:
	@echo 'devtools::load_all("."); tinytest::testthat(".")' | $(R)

slowtest: test-slow
test-slow:
	@$(R) -f inst/slowtest.R

test-all: test test-slow

coverage:
	@echo 'covr::report(file="tests/coverage-report.html")' | $(R)
	@open tests/coverage-report.html

.PHONY: man
man:
	@$(make_build_dir)
	@echo 'devtools::build_manual(".", path = "$(BUILD_DIR)")' | $(R)
	@open build/$(PKG)_$(GETVER).pdf
