# Introduction to Econometrics Using R

**Level:** Master I — Introductory Econometrics  
**Software:** R (≥ 4.0) and RStudio  
**License:** GNU General Public License v3

---

## About This Course

This repository contains the teaching material for an introductory econometrics course designed for first-year Master students in Economics. The course combines theoretical foundations with immediate practical application in R, following the progression of a standard introductory econometrics syllabus (Wooldridge, 2020).

By the end of the course, students will be able to:
- Import, clean, and explore economic datasets in R using the `tidyverse` ecosystem
- Produce publication-quality descriptive statistics and visualisations with `ggplot2`
- Estimate and interpret OLS regression models and perform standard diagnostic tests
- Correct for endogeneity using instrumental variables (IV / 2SLS)
- Handle panel data with fixed and random effects models
- Estimate binary choice models (Logit and Probit) and compute average partial effects

**Assumed background:** introductory statistics (hypothesis testing, confidence intervals) and basic linear algebra. No prior R programming experience is required.

---

## Course Contents

| Chapter | Topic |
|---------|-------|
| 1 | Introduction to R: syntax, objects, packages, built-in datasets |
| 2 | Data import (CSV, Excel) and management (joins, pooling) |
| 3 | Exploratory analysis of the crimes dataset (descriptive statistics, correlation, chi-square) |
| 4 | Data visualisation with base R and `ggplot2` |
| 5 | Ordinary Least Squares regression, model diagnostics (VIF, Durbin-Watson, heteroskedasticity tests) |
| 6 | Log-linear models, elasticity interpretation, model comparison |
| 7 | Instrumental Variables estimation (2SLS, Hausman endogeneity test, `ivreg`) |
| 8 | Panel data models: Fixed Effects, Random Effects, Hausman specification test |
| 9 | Binary choice models: Linear Probability Model, Logit, Probit, marginal effects |

---

## Repository Structure

```
IntroductionToEconometricsUsingR.md   ← full tutorial (read this on GitHub)
rCode/
  IntroductionToEconometricsUsingR.R  ← R script with all code, chapter by chapter
  RequiredPackages.r                  ← run once to install all dependencies
data/
  crimes.xls / crimes_C.csv / crimes_SC.csv / crimes_T.csv
  f2010.xls / f2012.xls / f2014.xls
  var1.xls / var2.xls
```

---

## Getting Started

1. **Install R and RStudio** (if not already installed).
2. **Clone or download** this repository.
3. **Open RStudio**, set the working directory to the project root:  
   `Session > Set Working Directory > To Project Directory`
4. **Install required packages** (run once):
   ```r
   source("rCode/RequiredPackages.r")
   ```
5. **Open the tutorial** `IntroductionToEconometricsUsingR.md` on GitHub (renders directly in the browser) or in RStudio's preview panel.
6. **Follow the code** in `rCode/IntroductionToEconometricsUsingR.R`, chapter by chapter.

---

## Datasets

| File | Description |
|------|-------------|
| `crimes.xls` | Primary dataset — 95 French metropolitan departments, 2011. Variables: total crimes, GDP, poverty index, unemployment rate, population, number of pupils, kilometres of road network. Source: DCPJ and INSEE, processed by ONDRP. |
| `f2010/f2012/f2014.xls` | Panel data across three fiscal years — used for data-pooling exercises (Chapter 2) and panel regression (Chapter 8). |
| `var1.xls`, `var2.xls` | Auxiliary tables — used to illustrate all six types of `dplyr` joins (Chapter 2). |

The `mroz` dataset used in Chapter 9 (binary choice models) is loaded directly from the `wooldridge` R package and does not require a local file.

---

## License

Distributed under the [GNU General Public License v3](LICENSE).
