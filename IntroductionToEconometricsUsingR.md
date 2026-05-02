# Introduction to Econometrics Using R

**Level:** Master I — Introductory Econometrics  
**Software:** R (≥ 4.0) and RStudio  
**License:** GNU General Public License v3

---

This course provides a hands-on introduction to econometric methods using R, designed for first-year Master students in Economics. It follows the standard progression of an introductory econometrics syllabus — from data handling and descriptive analysis through OLS, instrumental variables, panel data, and binary choice models — with every concept illustrated by immediately runnable R code.

The primary empirical application uses original cross-sectional data on **95 French metropolitan departments** (2011), covering crime rates, GDP, poverty, unemployment, and population. The panel and binary choice chapters use datasets from Wooldridge (2020), the standard reference textbook at this level.

**Assumed background:** introductory statistics (hypothesis testing, confidence intervals) and basic linear algebra. No prior R experience is required — Chapter 1 covers R from scratch.

**How to use this document:** read the narrative, then run the corresponding section of `rCode/IntroductionToEconometricsUsingR.R` in RStudio. The script is organised in the same chapter order as this tutorial.

---

## Table of Contents

- [Prerequisites and Setup](#prerequisites-and-setup)
- [Chapter 1 — Introduction to R](#chapter-1--introduction-to-r)
- [Chapter 2 — Importing and Managing Data](#chapter-2--importing-and-managing-data)
- [Chapter 3 — Exploring the Crimes Dataset](#chapter-3--exploring-the-crimes-dataset)
- [Chapter 4 — Data Visualization](#chapter-4--data-visualization)
- [Chapter 5 — Ordinary Least Squares Regression](#chapter-5--ordinary-least-squares-regression)
- [Chapter 6 — Log-Linear Models](#chapter-6--log-linear-models)
- [Chapter 7 — Instrumental Variables Estimation](#chapter-7--instrumental-variables-estimation)
- [Chapter 8 — Panel Data Models](#chapter-8--panel-data-models)
- [Chapter 9 — Binary Choice Models: Logit and Probit](#chapter-9--binary-choice-models-logit-and-probit)
- [References](#references)

---

## Prerequisites and Setup

### How to Use This Tutorial

The file `rCode/IntroductionToEconometricsUsingR.R` contains all the code shown in this document, organized in the same chapter order. It is recommended to read this tutorial and apply the corresponding code section by section in RStudio.

### Required Packages

Run the following script **once** to install all packages needed for the course:

```r
source("rCode/RequiredPackages.r")
```

The packages installed are:

| Package | Purpose |
|---------|---------|
| `tidyverse` | Data manipulation (`dplyr`), visualization (`ggplot2`), Excel import (`readxl`), model utilities (`broom`) |
| `questionr` | Frequency tables (`freq()`) |
| `gmodels` | Cross-tabulations (`CrossTable()`) |
| `GGally` | Scatter-plot matrix with correlations (`ggpairs()`) |
| `car` | Regression diagnostics (`vif()`, `ncvTest()`, `outlierTest()`) |
| `wooldridge` | Wooldridge textbook datasets |
| `AER` | `ivreg()`, `coeftest()`, `vcovHC()` |
| `plm` | Panel data models |
| `stargazer` | Formatted regression output tables |
| `margins` | Average partial effects for logit/probit |
| `pscl` | McFadden pseudo-R² (`pR2()`) |

### Setting the Working Directory

All data files are in the `data/` subfolder. Before running the code, set the working directory to the **project root** (the folder containing `rCode/`, `data/`, etc.):

```r
# In RStudio: Session > Set Working Directory > To Project Directory
# Or manually:
setwd("/path/to/IntroductionToEconometricsUsingR")

getwd()  # verify
```

---

## Chapter 1 — Introduction to R

### 1.1 Package Management

R packages extend the base language. They are hosted on CRAN (Comprehensive R Archive Network).

```r
options(repos = "https://cran.rstudio.com")

nrow(installed.packages())   # how many packages are installed locally
nrow(available.packages())   # how many are on CRAN
```

### 1.2 Loading Packages

A package must be **installed once** and **loaded every session**:

```r
# install.packages("tidyverse")   # run once; remove # to install
library(tidyverse)               # load into the current session

library(help = "tidyverse")      # package documentation
tidyverse_packages()             # list packages included in tidyverse
```

### 1.3 Getting Help

```r
help(mean)   # open the documentation for mean()
?mean        # shorthand for help()
```

### 1.4 R as a Calculator

```r
2 * 2
100 / 5
2^10
```

### 1.5 Objects and Assignment

All results in R are stored in **objects** using the assignment operator `<-`:

```r
a <- 2^10
a          # print the object
log(a)     # apply a function to it
```

### 1.6 Vectors and Data Frames

```r
# Create vectors
countries   <- c("France", "Germany", "Spain", "Italy")
CovidCases  <- c(1502763, 580415, 1331756, 759829)
CovidDeaths <- c(38289, 10904, 36495, 39412)

# Combine into a data frame
covidata <- data.frame(countries, CovidCases, CovidDeaths)

# Inspection functions
class(covidata)   # "data.frame"
str(covidata)     # column names, types, and first values
covidata          # print to the console
View(covidata)    # open in RStudio spreadsheet viewer
names(covidata)   # column names
```

```r
# Saving and removing a file (illustrative)
save(covidata, file = "covidata.RData")
unlink("covidata.RData")
```

### 1.7 Built-in Datasets

R ships with the `datasets` package containing many example datasets:

```r
library(datasets)
data()           # list all available datasets

data(mtcars)     # load the Motor Trend Cars dataset
?mtcars          # variable descriptions

str(mtcars)
head(mtcars, 10) # first 10 rows (default is 6)
tail(mtcars)     # last 6 rows
ncol(mtcars)
nrow(mtcars)
```

---

## Chapter 2 — Importing and Managing Data

### 2.1 Reading CSV Files

The `crimes` dataset is provided in three CSV formats and one Excel format, all in the `data/` folder. The different formats illustrate how `read.table()` handles various separators.

```r
# Semicolon-separated
crimes_sc <- read.table(file = "data/crimes_SC.csv", sep = ";",
                        header = TRUE, dec = ".")

# Comma-separated
crimes_c  <- read.table(file = "data/crimes_C.csv", sep = ",",
                        header = TRUE, dec = ".")

# Tab-separated
crimes_t  <- read.table(file = "data/crimes_T.csv", sep = "\t",
                        header = TRUE, dec = ".", quote = "\"")
```

### 2.2 Reading Excel Files

```r
library(readxl)

crimes_xls <- read_xls("data/crimes.xls", sheet = "data")
```

### 2.3 Creating Data by Hand (Tibbles)

A **tibble** is the tidyverse's version of a data frame — it prints more cleanly and never converts strings to factors automatically.

```r
# Dataset 1 — 10 animals, their young, and whether they are wild or domestic
data1 <- tibble(
  Animals = c("Cat", "Chicken", "Dog", "Elephant", "Giraffe",
              "Gnu", "Lion", "Panda", "Penguin", "Rabbit"),
  Baby    = c("Kitten", "Chick", "Puppy", "Calf", "Calf",
              "Calf", "Cub", "Cub", "Chick", "Kitten"),
  Nature  = c("Domestic", "Domestic", "Domestic", "Wild", "Wild",
              "Domestic", "Wild", "Wild", "Wild", "Domestic")
)

# Dataset 2 — same 10 animals + Turkey and Duck, with birth weight in kg
data2 <- tibble(
  Animals = c("Cat", "Chicken", "Dog", "Elephant", "Giraffe", "Gnu",
              "Lion", "Panda", "Penguin", "Rabbit", "Turkey", "Duck"),
  Baby    = c("Kitten", "Chick", "Puppy", "Calf", "Calf", "Calf",
              "Cub", "Cub", "Chick", "Kitten", "Poult", "Duckling"),
  Weight  = c(0.08, 0.09, 0.45, 90, 99, 36, 8, 0.1, 0.32, 0.03, 0.05, 0.05)
)
```

### 2.4 Joining Datasets

`dplyr` provides six join functions to merge two datasets on a common key (here, `Animals`). `data1` has 10 animals; `data2` has 12 (10 in common + Turkey and Duck).

| Function | Rows kept |
|----------|-----------|
| `left_join(data1, data2)` | All rows from `data1`; `Weight` added where available |
| `right_join(data1, data2)` | All rows from `data2`; `Nature` added where available |
| `inner_join(data1, data2)` | Only rows present in **both** (10 animals) |
| `full_join(data1, data2)` | All rows from **both** (12 animals) |
| `semi_join(data1, data2)` | Rows from `data1` that have a match in `data2`; no new columns |
| `anti_join(data1, data2)` | Rows from `data1` with **no** match in `data2` |

```r
data_left  <- left_join(data1, data2)
data_right <- right_join(data1, data2)
data_inner <- inner_join(data1, data2)
data_full  <- full_join(data1, data2)
data_semi  <- semi_join(data1, data2)
data_anti  <- anti_join(data1, data2)
```

### 2.5 Binding (Pooling) Datasets

When several datasets share the same variables but cover different time periods, `bind_rows()` stacks them vertically into a single panel dataset.

```r
library(readxl)
library(dplyr)

f2010 <- read_xls("data/f2010.xls")
f2012 <- read_xls("data/f2012.xls")
f2014 <- read_xls("data/f2014.xls")

dim(f2010)          # rows × columns
slice(f2010, 1:10)  # first 10 rows
str(f2010)          # structure

fpanel <- bind_rows(f2010, f2012, f2014)

dim(fpanel)   # total rows = sum of rows from all three files
head(fpanel, 10)
str(fpanel)
```

---

## Chapter 3 — Exploring the Crimes Dataset

The primary dataset `crimes.xls` contains cross-sectional data on **95 French metropolitan departments** for the year 2011. Data were collected from national sources (DCPJ, INSEE) and processed by the ONDRP.

| Variable | Description |
|----------|-------------|
| `dep_code` | Department code |
| `dep_name` | Department name |
| `region_name` | Region name (22 regions, pre-2016 reform) |
| `big_region` | Aggregated into 6 groups: North-East, North-West, Île-de-France, Centre, South-East, South-West |
| `crimes` | Total number of crimes, incidents and offences recorded (2011) |
| `gdp_2011` | GDP (millions of euros, 2011) |
| `poverty_index` | Poverty intensity index |
| `population` | Total population |
| `unemp_rate` | Unemployment rate (%) |
| `pupils` | Number of registered students |
| `roads` | Total length of the departmental road network (kilometres) |

### 3.1 Loading and Inspecting the Data

```r
library(readxl)

crimes <- read_xls("data/crimes.xls", sheet = "data", col_names = TRUE)

dim(crimes)          # 95 rows × number of columns
slice(crimes, 1:10)  # first 10 observations
names(crimes)        # column names
str(crimes)          # data types
```

### 3.2 Descriptive Statistics for a Single Variable

```r
min(crimes$poverty_index,    na.rm = TRUE)
max(crimes$poverty_index,    na.rm = TRUE)
mean(crimes$poverty_index,   na.rm = TRUE)
sd(crimes$poverty_index,     na.rm = TRUE)
median(crimes$poverty_index, na.rm = TRUE)
quantile(crimes$poverty_index, na.rm = TRUE)   # min, Q1, Q2, Q3, max
```

### 3.3 Summary Statistics for Multiple Variables

```r
crimes_quant <- subset(crimes,
                       select = c(crimes, gdp_2011, poverty_index,
                                  population, unemp_rate))
summary(crimes_quant)
```

### 3.4 Frequency Table for a Categorical Variable

```r
# Count of departments per big region (base R)
table(crimes$big_region)

# Frequency table with relative frequencies and cumulative column
library(questionr)
freq(crimes$big_region, cum = TRUE, sort = TRUE, valid = FALSE,
     total = TRUE, na.last = TRUE)
```

### 3.5 Identifying Extreme Observations

```r
subset(crimes, crimes == min(crimes$crimes), select = c(dep_name, crimes))
subset(crimes, crimes == max(crimes$crimes), select = c(dep_name, crimes))
```

### 3.6 Correlation Analysis

The **Pearson** coefficient measures linear association; the **Spearman** coefficient is non-parametric and suitable when the relationship may be non-linear or when outliers are present.

```r
cor(crimes$crimes, crimes$poverty_index, method = "pearson")

# Test H0: rho = 0 (no linear correlation)
cor.test(crimes$crimes, crimes$poverty_index,
         alternative = "two.sided", method = "pearson", conf.level = 0.95)

# Non-parametric Spearman correlation
cor.test(crimes$crimes, crimes$poverty_index,
         alternative = "two.sided", method = "spearman", conf.level = 0.95)
```

### 3.7 Creating Binary (Dummy) Variables

```r
# 1 if crimes exceed the sample mean, 0 otherwise
crimes$crim_degree <- ifelse(crimes$crimes > mean(crimes$crimes), 1, 0)

# 1 if poverty_index exceeds the sample mean, 0 otherwise
crimes$poverty_cat <- ifelse(crimes$poverty_index > mean(crimes$poverty_index), 1, 0)
```

### 3.8 Cross-Tabulation and Chi-Square Test

The chi-square test of independence examines whether the level of crime and the level of poverty are associated across departments.

- $H_0$: No association between crime level and poverty level
- $H_1$: The two variables are associated

```r
library(gmodels)

# Descriptive cross-tabulation
CrossTable(crimes$crim_degree, crimes$poverty_cat,
           digits = 2, prop.chisq = FALSE, chisq = FALSE,
           fisher = FALSE, mcnemar = FALSE,
           missing.include = FALSE, format = "SPSS")

# With chi-square test
CrossTable(crimes$crim_degree, crimes$poverty_cat,
           digits = 2, prop.chisq = FALSE, chisq = TRUE,
           fisher = FALSE, mcnemar = FALSE,
           missing.include = FALSE, format = "SPSS")
```

---

## Chapter 4 — Data Visualization

### 4.1 Base R Scatter Plot

```r
plot(crimes$poverty_index,
     main = "Poverty Index by Department",
     xlab = "Department Index",
     ylab = "Poverty Index Value",
     type = "p", pch = 20, cex = 1, col = "blue")
```

Department labels can be added using `text()`, called on a **separate line** — unlike `ggplot2`, base R graphics cannot be extended with `+`:

```r
plot(crimes$poverty_index,
     main = "Poverty Index by Department",
     xlab = "Department Index",
     ylab = "Poverty Index Value",
     type = "p")
text(crimes$poverty_index,
     labels = crimes$dep_name,
     cex = 0.7, pos = 1)
```

### 4.2 Histogram and Kernel Density Estimation

A **histogram** divides the range into bins and counts observations; a **kernel density** estimates the underlying continuous distribution.

```r
hist(crimes$poverty_index,
     main = "Distribution of the Poverty Index",
     xlab = "Poverty Index")

dens <- density(crimes$poverty_index, bw = "nrd0",
                kernel = "epanechnikov", na.rm = TRUE)
plot(dens, frame = TRUE, col = "steelblue",
     main = "Kernel Density — Poverty Index")
```

### 4.3 Scatter Plots with ggplot2

`ggplot2` builds plots in layers. `geom_smooth()` adds a regression fit.

```r
library(ggplot2)

# Basic scatter plot
ggplot(crimes, aes(x = poverty_index, y = crimes)) +
  geom_point() +
  labs(title    = "Scatterplot",
       subtitle = "Relationship between Crimes and Poverty",
       x        = "Poverty Index",
       y        = "Number of Crimes",
       caption  = "Source: DCPJ, processed by ONDRP.")

# OLS regression line with 95% confidence band
ggplot(crimes, aes(x = poverty_index, y = crimes)) +
  geom_point() +
  geom_smooth(method = lm, se = TRUE) +
  labs(title = "Scatterplot with OLS Fit",
       x = "Poverty Index", y = "Number of Crimes",
       caption = "Source: DCPJ, processed by ONDRP.")

# LOESS smoothing — captures non-linear patterns
ggplot(crimes, aes(x = poverty_index, y = crimes)) +
  geom_point() +
  geom_smooth(method = loess, se = TRUE) +
  labs(title = "Scatterplot with LOESS Fit",
       x = "Poverty Index", y = "Number of Crimes",
       caption = "Source: DCPJ, processed by ONDRP.")
```

### 4.4 Correlation Matrix and Scatter-Plot Matrix

```r
# Working dataset — include big_region for grouping
crim_tab <- subset(crimes,
                   select = c(crimes, gdp_2011, poverty_index,
                              unemp_rate, population, big_region))

summary(crim_tab, digits = 2)

# Correlation matrix: use numeric columns only (big_region is character)
crim_num <- select(crim_tab, where(is.numeric))
round(cor(crim_num, method = "pearson"), 2)

# Pairs scatter plot
plot(crim_num)

# Combined scatter-plot matrix, histograms, and correlation coefficients
library(GGally)
ggpairs(crim_tab, columns = 1:5, aes(colour = big_region, alpha = 0.5))
```

> **Note on `attach()`:** Some older R tutorials use `attach(dataset)` to make column names available as standalone variables. This is **not recommended** in modern R because it can cause unpredictable name conflicts. Use `dataset$column` or `with(dataset, expression)` instead.

---

## Chapter 5 — Ordinary Least Squares Regression

### 5.1 Model Specification

We estimate the following linear model relating the total number of crimes to four departmental characteristics:

$$\text{crimes}_i = \beta_0 + \beta_1 \cdot \text{gdp\_2011}_i + \beta_2 \cdot \text{poverty\_index}_i + \beta_3 \cdot \text{unemp\_rate}_i + \beta_4 \cdot \text{population}_i + \varepsilon_i$$

```r
MyModel <- lm(crimes ~ gdp_2011 + poverty_index + unemp_rate + population,
              data = crimes)
names(MyModel)   # components of a lm object (coefficients, residuals, fitted.values, …)
```

### 5.2 Results and Interpretation

```r
options(scipen = 999)   # disable scientific notation in output
summary(MyModel)
```

The output reports:
- **Coefficients** — estimate, standard error, t-statistic, and p-value for each regressor
- **R²** — share of the variance of `crimes` explained by the model
- **Adjusted R²** — penalises for additional regressors
- **F-statistic** — tests joint significance of all regressors ($H_0$: all $\beta = 0$)

### 5.3 Confidence Intervals and Fitted Values

```r
confint(MyModel, level = 0.95)   # 95% CIs for each coefficient

crimes_hat <- fitted(MyModel)    # in-sample predicted values
head(MyModel$residuals)          # first residuals ê = y − ŷ
```

### 5.4 Multicollinearity — Variance Inflation Factor (VIF)

High correlation among regressors inflates standard errors, making individual coefficients imprecisely estimated.

```r
library(car)

vif(MyModel)

# Rule of thumb: VIF > 4 warrants inspection; VIF > 10 indicates a serious problem
sqrt(vif(MyModel)) > 2
```

### 5.5 Autocorrelation — Durbin-Watson Test

The Durbin-Watson statistic tests for first-order autocorrelation in the residuals. Values near 2 indicate no autocorrelation; values near 0 (positive) or 4 (negative) indicate a problem.

```r
durbinWatsonTest(MyModel)
```

### 5.6 Heteroskedasticity — Non-Constant Variance Test

OLS standard errors assume constant error variance (homoskedasticity). If the variance of $\varepsilon$ changes with the fitted values, inference is unreliable.

- $H_0$: Constant error variance (homoskedasticity)  
- $H_1$: Error variance changes with fitted values (heteroskedasticity)

```r
ncvTest(MyModel)
```

A significant result (small p-value) rejects $H_0$. The remedy is to use **heteroskedasticity-consistent (HC) standard errors**, which correct the variance-covariance matrix of the estimator without changing the point estimates:

```r
library(AER)

# HC1-robust standard errors (equivalent to Stata's robust option)
coeftest(MyModel, vcov = vcovHC, type = "HC1")
```

The coefficients are identical to those of `summary(MyModel)`, but the standard errors, t-statistics, and p-values are now robust to heteroskedasticity.

### 5.7 Residual Diagnostic Plots

```r
# Spread-level plot: residual magnitude vs. fitted values
spreadLevelPlot(MyModel)

# QQ plot of studentised residuals — assess normality of errors
qqPlot(MyModel, main = "QQ Plot — Studentised Residuals")

# Added-variable plots — detect influential observations
leveragePlots(MyModel)

# Bonferroni test: identifies the single most extreme outlier
outlierTest(MyModel)
```

---

## Chapter 6 — Log-Linear Models

### 6.1 Motivation for Log Transformation

The total number of crimes and GDP both vary enormously across the 95 departments (a few large cities vs. many small rural departments). This scale heterogeneity tends to produce heteroskedasticity and non-linearity. Expressing variables in **per-capita logarithms**:

1. Compresses the right tail and reduces the influence of outliers
2. Controls for population size automatically
3. Produces **elasticity** coefficients: $\hat{\beta}$ is the expected percentage change in the dependent variable for a 1% change in the regressor

### 6.2 Creating Log-Transformed Variables

```r
crimes$lcrimes_cap <- log(crimes$crimes / crimes$population)
crimes$lgdp_cap    <- log(crimes$gdp_2011 * 1e6 / crimes$population)
crimes$lpop        <- log(crimes$population)
```

### 6.3 Estimation

$$\ln\!\left(\frac{\text{crimes}}{pop}\right)_i = \beta_0 + \beta_1 \ln\!\left(\frac{\text{gdp}}{pop}\right)_i + \beta_2 \cdot \text{poverty\_index}_i + \beta_3 \cdot \text{unemp\_rate}_i + \beta_4 \ln(pop)_i + \varepsilon_i$$

```r
MyModel_log <- lm(lcrimes_cap ~ lgdp_cap + poverty_index + unemp_rate + lpop,
                  data = crimes)
summary(MyModel_log)
```

**Interpretation of $\hat{\beta}_1$:** A 1% increase in GDP per capita is associated with a $\hat{\beta}_1$% change in the crime rate per capita, holding other variables constant.

### 6.4 Model Comparison

```r
library(broom)

glance(MyModel)       # R², AIC, BIC, sigma for the level model
glance(MyModel_log)   # same for the log model
```

> **Caution:** R² cannot be compared between a model with `crimes` and one with `lcrimes_cap` as the dependent variable — they measure explained variance on different scales. Use **AIC** or **BIC** for model selection between nested or non-nested models with the same dependent variable.

---

## Chapter 7 — Instrumental Variables Estimation

### 7.1 The Endogeneity Problem

In the log-log model, GDP per capita (`lgdp_cap`) may be **endogenous** due to an **omitted variable problem**. Specifically, the model fails to capture unobserved structural factors, such as the fundamental nature of the departmental economy (e.g., whether a department is predominantly rural or heavily industrial) and the specific sectoral composition of its activities. Because these omitted characteristics simultaneously influence both the region's economic performance and its crime rates, they violate the OLS exogeneity assumption. Consequently, the OLS estimator would be biased and inconsistent.

### 7.2 Instrumental Variables: Identification Strategy

To identify the causal effect of GDP per capita on crime in the presence of omitted variable bias, we propose **two instrumental variables** for `lgdp_cap`: the student population (`pupils`) and the kilometres of departmental road networks (`roads`).

The validity of these instruments rests on two standard identification conditions:

**1. Relevance** — $\text{Cov}(Z, X) \neq 0$: each instrument must be correlated with the endogenous regressor `lgdp_cap`. This is theoretically supported by dual channels:
- *Human capital accumulation:* student enrolment (`pupils`) is linked to economic output through the accumulation of skilled labour, which raises departmental productivity and GDP per capita.
- *Infrastructural capacity:* road networks (`roads`) facilitate regional trade, reduce transport costs, and attract economic activity, thereby positively affecting GDP per capita.

**2. Exclusion restriction** — $\text{Cov}(Z, u) = 0$: each instrument must be orthogonal to the structural error term, i.e., it must affect the crime rate **only through** GDP per capita, and not through the unobserved omitted factors.
- For `pupils`: student enrolment is largely independent of the unobserved regional characteristics — such as the specific rural or industrial sectoral composition — that jointly confound GDP and crime.
- For `roads`: the historical layout of departmental road infrastructure affects crime exclusively through its impact on current economic activity, rather than being determined by the contemporary unobserved structural factors that drive local crime rates.

> **Overidentification:** with two instruments and one endogenous variable, the model is **overidentified** (degree of overidentification = 1). This allows the **Sargan–Hansen test** to partially assess instrument validity: conditional on one instrument being valid, it tests whether the other is also exogenous. The test is automatically reported by `summary(ivreg_est, diagnostics = TRUE)`.

### 7.3 Two-Stage Least Squares (2SLS) — Manual Approach

**Stage 1 — Reduced form:** regress the endogenous variable on all exogenous variables **and both instruments**.

```r
library(AER)

reducedForm <- lm(lgdp_cap ~ pupils + roads + unemp_rate + poverty_index + lpop,
                  data = crimes)
summary(reducedForm)

# Test the joint significance of both instruments (H0: their coefficients are jointly = 0)
# A significant F-statistic confirms instrument relevance
coeftest(reducedForm, vcov = vcovHC, type = "HC1")

summary(reducedForm)$r.squared   # first-stage R²: should be reasonably high

lgdp_pred <- fitted(reducedForm)   # predicted values of lgdp_cap
```

**Stage 2 — Structural equation:** replace `lgdp_cap` with its Stage 1 prediction.

```r
structuralEq <- lm(lcrimes_cap ~ lgdp_pred + unemp_rate + poverty_index + lpop,
                   data = crimes)
summary(structuralEq)
coeftest(structuralEq, vcov = vcovHC, type = "HC1")
```

> **Important:** The standard errors from manual 2SLS are biased (they ignore the uncertainty from Stage 1). Always use `ivreg()` for correct inference.

### 7.4 Hausman Endogeneity Test

The Hausman test compares OLS and IV estimates: if they differ significantly, the OLS assumption of exogeneity is rejected.

The **regression-based version** adds the Stage 1 **residuals** (not fitted values) to the structural equation. A significant coefficient on the residuals confirms endogeneity.

```r
HausmanTest   <- lm(lgdp_cap ~ pupils + roads + unemp_rate + poverty_index + lpop,
                    data = crimes)
gdp_residuals <- residuals(HausmanTest)   # first-stage residuals

endoTest <- lm(lcrimes_cap ~ lgdp_cap + unemp_rate + poverty_index + lpop +
                 gdp_residuals, data = crimes)
coeftest(endoTest, vcov = vcovHC, type = "HC1")
```

### 7.5 IV Estimation with `ivreg()`

`ivreg()` implements 2SLS in one step with correct standard errors. The formula syntax `y ~ endogenous + exogenous | instruments + exogenous` separates structural regressors (left of `|`) from instruments (right of `|`). Exogenous regressors appear on **both sides**. Both `pupils` and `roads` appear on the right side only.

```r
ivreg_est <- ivreg(lcrimes_cap ~ lgdp_cap + unemp_rate + poverty_index + lpop |
                                 pupils + roads + unemp_rate + poverty_index + lpop,
                   data = crimes)

# diagnostics = TRUE reports: weak-instrument F-test, Wu-Hausman endogeneity test,
# and the Sargan overidentification test (valid because we have 2 instruments for 1
# endogenous variable)
summary(ivreg_est, diagnostics = TRUE)

coeftest(ivreg_est, vcov = vcovHC, type = "HC1")
```

**Reading the diagnostic tests:**

| Diagnostic | $H_0$ | Reject $H_0$ means |
|------------|--------|---------------------|
| Weak instruments (F-test) | Instruments are weak ($F < 10$) | Instruments are relevant |
| Wu-Hausman | `lgdp_cap` is exogenous | Endogeneity confirmed → prefer IV |
| Sargan | Both instruments are exogenous | At least one instrument is invalid |

---

## Chapter 8 — Panel Data Models

### 8.1 Introduction

Panel data follow the same individuals over multiple time periods. They allow controlling for **unobserved time-invariant heterogeneity** (individual fixed effects), which would otherwise bias cross-sectional OLS estimates.

We use the **`wagepan`** dataset from Wooldridge (Chapter 14): 545 U.S. men observed annually from 1980 to 1987 (4,360 observations). The dependent variable is log hourly wage (`lwage`).

### 8.2 Fixed Effects (FE) Model

The fixed effects model includes an individual-specific intercept $\alpha_i$ that captures all time-invariant characteristics of individual $i$ (ability, family background, etc.). Estimation uses **within-group demeaning** (subtracting individual means), which eliminates $\alpha_i$.

$$\text{lwage}_{it} = \alpha_i + \beta_1 \text{exper}_{it} + \beta_2 \text{expersq}_{it} + \beta_3 \text{union}_{it} + \cdots + \varepsilon_{it}$$

```r
library(wooldridge)
library(plm)

data(wagepan)

panel.fe <- plm(lwage ~ educ + exper + expersq + union + south + married + black,
                data  = wagepan,
                index = c("nr", "year"),
                model = "within")
summary(panel.fe)
```

> **Note:** Time-invariant variables (`educ`, `black`) cannot be estimated by FE — they are absorbed into $\alpha_i$.

### 8.3 Random Effects (RE) Model

The random effects model treats $\alpha_i$ as a random variable drawn from a distribution, **assumed uncorrelated with the regressors**. It is more efficient than FE and can estimate time-invariant coefficients.

```r
panel.re <- plm(lwage ~ educ + exper + expersq + union + south + married + black,
                data  = wagepan,
                index = c("nr", "year"),
                model = "random")
summary(panel.re)
```

### 8.4 Hausman Specification Test

The Hausman test chooses between FE and RE by testing whether the individual effects are correlated with the regressors.

- $H_0$: No correlation → RE is consistent and more efficient
- $H_1$: Correlation exists → Only FE is consistent

```r
phtest(panel.fe, panel.re)
```

Reject $H_0$ (small p-value) → prefer **Fixed Effects**.

### 8.5 Regression Output with `stargazer`

```r
library(stargazer)

stargazer(MyModel, MyModel_log, structuralEq,
          title = "Regression Results",
          align = TRUE,
          type  = "text")   # type = "latex" for a LaTeX document
```

---

## Chapter 9 — Binary Choice Models: Logit and Probit

### 9.1 Introduction

When the dependent variable $y_i$ takes only two values — 0 or 1 — the **Linear Probability Model (LPM)** estimated by OLS has well-known limitations:

- Predicted probabilities can fall outside the $[0, 1]$ interval
- The linear specification is implausible: the marginal effect of a regressor on the probability should diminish as the probability approaches 0 or 1
- Residuals are heteroskedastic by construction

**Logit** and **Probit** models resolve these issues by modelling the probability through a monotone S-shaped (sigmoid) function $F(\cdot)$ bounded in $(0, 1)$:

$$P(y_i = 1 \mid \mathbf{x}_i) = F(\mathbf{x}_i' \boldsymbol{\beta})$$

| Model | Link function $F(z)$ |
|-------|----------------------|
| Logit | Logistic CDF: $\displaystyle\frac{e^z}{1 + e^z}$ |
| Probit | Standard normal CDF: $\Phi(z)$ |

Both models are estimated by **Maximum Likelihood (MLE)**.

### 9.2 The `mroz` Dataset

We use the **`mroz`** dataset from Wooldridge (1987), available in the `wooldridge` R package. It contains data on **753 married women** from the 1975 wave of the Panel Study of Income Dynamics (PSID).

The binary outcome of interest is **`inlf`** (= 1 if the woman was in the labor force in 1975, 0 otherwise).

| Variable | Description |
|----------|-------------|
| `inlf` | **Dependent variable** — 1 if in labor force |
| `nwifeinc` | Non-wife household income (thousands of USD) |
| `educ` | Years of schooling |
| `exper` | Actual labor market experience (years) |
| `expersq` | `exper²` (captures diminishing returns) |
| `age` | Age in years |
| `kidslt6` | Number of children under 6 years old |
| `kidsge6` | Number of children aged 6–18 |

**Source:** Mroz, T. A. (1987). "The sensitivity of an empirical model of married women's hours of work to economic and statistical assumptions." *Econometrica*, 55(4), 765–799. Reproduced in Wooldridge (2002), *Introductory Econometrics*, Example 17.1.

```r
library(wooldridge)
data(mroz)

str(mroz)
head(mroz)

mean(mroz$inlf)   # ~0.57 — about 57% of women in the sample are in the labor force

summary(mroz[, c("inlf", "nwifeinc", "educ", "exper", "age", "kidslt6", "kidsge6")])
```

### 9.3 Linear Probability Model (Baseline)

Despite its limitations, the LPM is a useful benchmark and is easy to interpret.

```r
lpm <- lm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
          data = mroz)
summary(lpm)
```

**Interpretation:** $\hat{\beta}_j$ is the estimated **change in the probability of labor force participation** for a one-unit increase in $x_j$, holding all other variables constant.

```r
pred_lpm <- fitted(lpm)
sum(pred_lpm < 0 | pred_lpm > 1)   # count predictions outside [0, 1]
```

### 9.4 Logit Model

```r
logit_model <- glm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
                   data   = mroz,
                   family = binomial(link = "logit"))
summary(logit_model)
```

The estimated coefficients are on the **log-odds** scale:

$$\ln\!\left(\frac{P(y_i = 1)}{1 - P(y_i = 1)}\right) = \mathbf{x}_i' \boldsymbol{\hat{\beta}}$$

A **positive** coefficient means the log-odds — and therefore the probability — of $y = 1$ increase with the regressor. The size of the probability effect depends non-linearly on all other variables (see Section 9.5).

### 9.5 Probit Model

```r
probit_model <- glm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
                    data   = mroz,
                    family = binomial(link = "probit"))
summary(probit_model)
```

The probit coefficients are on the **standard normal index** scale. The signs and significance of coefficients should agree with those of the logit. The logit coefficients are approximately $\pi / \sqrt{3} \approx 1.81$ times larger in magnitude than the probit ones, due to the different scales of the two distributions. To interpret the **magnitude** of effects, compute marginal effects (Section 9.6).

### 9.6 Marginal Effects

Because the relationship between $\mathbf{x}$ and $P(y = 1)$ is non-linear, the marginal effect of a regressor $x_j$ on the probability varies across individuals. The most commonly reported quantity is the **Average Partial Effect (APE)** — the derivative of $P(y = 1 | \mathbf{x})$ with respect to $x_j$, averaged across all observations:

$$\text{APE}_j = \frac{1}{n} \sum_{i=1}^{n} f(\mathbf{x}_i' \hat{\boldsymbol{\beta}}) \cdot \hat{\beta}_j$$

where $f(\cdot)$ is the density corresponding to $F(\cdot)$ (logistic density for logit; standard normal density for probit).

#### Using the `margins` package

```r
library(margins)

ape_logit <- margins(logit_model)
summary(ape_logit)

ape_probit <- margins(probit_model)
summary(ape_probit)
```

#### Manual computation for education (logit) — for understanding

For the logit model, the logistic density at observation $i$ is $f(\hat{z}_i) = \hat{p}_i(1 - \hat{p}_i)$.

```r
p_hat     <- fitted(logit_model)
f_logit   <- p_hat * (1 - p_hat)

beta_educ <- coef(logit_model)["educ"]
ape_educ  <- mean(f_logit * beta_educ)
ape_educ
```

**Interpretation:** `ape_educ` is the estimated average change in the probability of labor force participation for one additional year of schooling.

### 9.7 Model Evaluation

#### McFadden Pseudo-R²

OLS's R² is not defined for MLE models. The **McFadden pseudo-R²** is the most widely used alternative:

$$R^2_{\text{McFadden}} = 1 - \frac{\ln \hat{L}_{\text{full}}}{\ln \hat{L}_{\text{null}}}$$

where $\hat{L}_{\text{full}}$ is the log-likelihood of the estimated model and $\hat{L}_{\text{null}}$ is the log-likelihood of an intercept-only model. Values between 0.2 and 0.4 are generally considered to indicate good fit.

```r
library(pscl)

pR2(logit_model)    # McFadden and other pseudo-R² measures
pR2(probit_model)
```

#### AIC and BIC

```r
AIC(logit_model, probit_model)
BIC(logit_model, probit_model)
```

Lower AIC or BIC indicates a better model. These criteria can be compared across logit and probit (both estimated by MLE), but **not** between an OLS model (LPM) and a MLE model.

#### Confusion Matrix and Classification Accuracy

```r
# Classify using a 0.5 probability threshold
class_logit  <- ifelse(fitted(logit_model)  > 0.5, 1, 0)
class_probit <- ifelse(fitted(probit_model) > 0.5, 1, 0)

# Confusion matrices
table(Predicted = class_logit,  Actual = mroz$inlf)
table(Predicted = class_probit, Actual = mroz$inlf)

# Percentage correctly classified
mean(class_logit  == mroz$inlf)
mean(class_probit == mroz$inlf)
```

### 9.8 Odds Ratios (Logit Only)

Exponentiating a logit coefficient gives the **odds ratio** — the multiplicative change in the odds of $y = 1$ for a one-unit increase in $x_j$, holding other variables constant:

$$\text{OR}_j = e^{\hat{\beta}_j}$$

```r
exp(coef(logit_model))       # odds ratios
exp(confint(logit_model))    # 95% confidence intervals for odds ratios
```

**Interpretation example:** If $\widehat{OR}_{\text{kidslt6}} = 0.35$, having one additional child under the age of 6 **multiplies the odds** of labor force participation by 0.35 — a 65% reduction in the odds — holding other variables constant.

### 9.9 Comparing LPM, Logit, and Probit

```r
library(stargazer)

stargazer(lpm, logit_model, probit_model,
          title           = "Labor Force Participation — Binary Choice Models",
          dep.var.labels  = "In Labor Force (1 = Yes)",
          column.labels   = c("LPM (OLS)", "Logit", "Probit"),
          covariate.labels = c("Non-wife income", "Education", "Experience",
                               "Experience²", "Age", "Kids < 6", "Kids 6–18",
                               "Intercept"),
          type = "text")   # type = "latex" for a LaTeX document
```

**Summary of model properties:**

| Feature | LPM | Logit | Probit |
|---------|-----|-------|--------|
| Estimation method | OLS | MLE | MLE |
| Predictions bounded in [0,1] | No | Yes | Yes |
| Coefficient interpretation | Marginal effect (direct) | Log-odds | Normal index |
| Marginal effects | Equal to coefficients | Non-linear, use `margins()` | Non-linear, use `margins()` |
| Goodness-of-fit measure | R² | Pseudo R², AIC/BIC | Pseudo R², AIC/BIC |
| Handles rare events well | No | Reasonably | Yes |

**General guidance:**

- The LPM is adequate when predicted probabilities stay well within (0.1, 0.9) and the focus is on estimating average partial effects.
- Logit and Probit give nearly identical results in most empirical applications; the choice between them is largely conventional.
- Use the LPM as a robustness check alongside Logit or Probit.

---

## References

- Wooldridge, J. M. (2020). *Introductory Econometrics: A Modern Approach* (7th ed.). Cengage Learning.
- Kleiber, C., & Zeileis, A. (2008). *Applied Econometrics with R*. Springer-Verlag. (`AER` package)
- Fox, J., & Weisberg, S. (2019). *An R Companion to Applied Regression* (3rd ed.). Sage Publications. (`car` package)
- Wickham, H., et al. (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686.
- Croissant, Y., & Millo, G. (2008). Panel Data Econometrics in R: The plm Package. *Journal of Statistical Software*, 27(2).
- Mroz, T. A. (1987). The sensitivity of an empirical model of married women's hours of work to economic and statistical assumptions. *Econometrica*, 55(4), 765–799.

---

*Distributed under the GNU General Public License v3. See `LICENSE` for details.*
