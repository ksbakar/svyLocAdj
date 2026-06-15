
# svyLocAdj

`svyLocAdj` is an R package for adjusting and analyzing Demographic and Health Survey (DHS) data using Bayesian modeling approaches implemented through Stan.

## Prerequisites

Before installing `svyLocAdj`, users must ensure that:

1. A working C++ compiler/toolchain is installed.
2. The `rstan` package is successfully installed and configured.

Since `svyLocAdj` uses Stan for Bayesian computation, proper compiler and `rstan` installation are required.

---

## Step 1: Verify Compiler Installation

### Windows

Install **Rtools** corresponding to your version of R:

https://cran.r-project.org/bin/windows/Rtools/

After installation, verify that the compiler is available:

```r
install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```

A return value of `TRUE` indicates that the build tools are correctly configured.

### macOS

Install Apple's Command Line Tools:

```bash
xcode-select --install
```

### Linux

For Ubuntu/Debian systems:

```bash
sudo apt update
sudo apt install build-essential
```

For other Linux distributions, install the equivalent GNU compiler tools.

---

## Step 2: Install RStan

Install `rstan` and its dependencies:

```r
install.packages(c("StanHeaders", "rstan"), dependencies = TRUE)
```

Load `rstan` and verify that it is working correctly:

```r
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

For detailed installation instructions and platform-specific guidance, see:

https://mc-stan.org/users/interfaces/rstan

---

## Step 3: Install devtools

If `devtools` is not already installed:

```r
install.packages("devtools")
```

---

## Step 4: Install svyLocAdj from GitHub

Install the latest development version from GitHub:

```r
devtools::install_github("USERNAME/svyLocAdj")
```

Replace `USERNAME` with the GitHub username or organization hosting the repository.

---

## Step 5: Load the Package

```r
library(svyLocAdj)
```

Confirm the installation:

```r
packageVersion("svyLocAdj")
```

---

## Troubleshooting

### Stan Model Compilation Errors

If Stan models fail to compile:

1. Verify that the compiler toolchain is installed correctly.
2. Restart R after installing Rtools (Windows) or Command Line Tools (macOS).
3. Reinstall or update `rstan` and `StanHeaders`:

```r
install.packages(c("StanHeaders", "rstan"))
```

4. Check compiler availability:

```r
pkgbuild::has_build_tools(debug = TRUE)
```

### Windows-Specific Issues

Ensure that:

* The correct version of Rtools is installed.
* Rtools is available on the system PATH.
* Your R version is compatible with the installed Rtools version.

---

## Example

```r
library(svyLocAdj)

# Example workflow
# ... will add package-specific examples here soon ...
```

## Reporting Issues

If you encounter bugs or installation problems, please open an issue in the GitHub repository.

## License

See the `LICENSE` file for licensing information.


