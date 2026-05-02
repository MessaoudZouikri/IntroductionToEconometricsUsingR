
# =============================================================================
# Introduction to Econometrics Using R
# Master I — Introductory Econometrics Course
# =============================================================================
# Before running this script, set the working directory to the project root:
#   Session > Set Working Directory > To Project Directory   (RStudio menu)
# All data files are expected in the data/ subfolder.
# Run rCode/RequiredPackages.r once to install all dependencies.
# =============================================================================


# =============================================================================
# CHAPTER 1 — Introduction to R
# =============================================================================

# Set the CRAN mirror (avoids interactive prompts in batch mode)
options(repos = "https://cran.rstudio.com")

# Count installed packages
nrow(installed.packages())

# Count packages currently available on CRAN
nrow(available.packages())

# Install tidyverse (remove the # to run)
# install.packages("tidyverse")

# Load tidyverse into memory
library(tidyverse)

# Package documentation
library(help = "tidyverse")
tidyverse_packages()

# Getting help on a function
help(mean)
?mean

# R as a calculator
2 * 2
100 / 5
2^10

# Assignment operator <-
a <- 2^10
a
log(a)

# Creating vectors
countries   <- c("France", "Germany", "Spain", "Italy")
CovidCases  <- c(1502763, 580415, 1331756, 759829)
CovidDeaths <- c(38289, 10904, 36495, 39412)

# Combine into a data frame
covidata <- data.frame(countries, CovidCases, CovidDeaths)

# Save and remove (for illustration)
save(covidata, file = "covidata.RData")
unlink("covidata.RData")

# Inspect the data frame
class(covidata)
str(covidata)
covidata
View(covidata)
names(covidata)

# Built-in datasets
library(datasets)
library(help = "datasets")
data()

data(mtcars)
?mtcars
str(mtcars)
head(mtcars, 10)
tail(mtcars)
ncol(mtcars)
nrow(mtcars)


# =============================================================================
# CHAPTER 2 — Importing and Managing Data
# =============================================================================

# --- 2.1  Reading CSV files ---------------------------------------------------

# Semicolon-separated
crimes_sc <- read.table(file = "data/crimes_SC.csv", sep = ";",
                        header = TRUE, dec = ".")

# Comma-separated
crimes_c  <- read.table(file = "data/crimes_C.csv", sep = ",",
                        header = TRUE, dec = ".")

# Tab-separated
crimes_t  <- read.table(file = "data/crimes_T.csv", sep = "\t",
                        header = TRUE, dec = ".", quote = "\"")

# --- 2.2  Reading Excel files -------------------------------------------------

library(readxl)
help(package = "readxl")

crimes_xls <- read_xls("data/crimes.xls", sheet = "data")

# --- 2.3  Creating data by hand (tibble) -------------------------------------

# data1: 10 animals, their young, and whether they are wild or domestic
data1 <- tibble(
  Animals = c("Cat", "Chicken", "Dog", "Elephant", "Giraffe",
              "Gnu", "Lion", "Panda", "Penguin", "Rabbit"),
  Baby    = c("Kitten", "Chick", "Puppy", "Calf", "Calf",
              "Calf", "Cub", "Cub", "Chick", "Kitten"),
  Nature  = c("Domestic", "Domestic", "Domestic", "Wild", "Wild",
              "Domestic", "Wild", "Wild", "Wild", "Domestic")
)

# data2: same animals + Turkey and Duck, with birth weight (kg)
data2 <- tibble(
  Animals = c("Cat", "Chicken", "Dog", "Elephant", "Giraffe", "Gnu",
              "Lion", "Panda", "Penguin", "Rabbit", "Turkey", "Duck"),
  Baby    = c("Kitten", "Chick", "Puppy", "Calf", "Calf", "Calf",
              "Cub", "Cub", "Chick", "Kitten", "Poult", "Duckling"),
  Weight  = c(0.08, 0.09, 0.45, 90, 99, 36, 8, 0.1, 0.32, 0.03, 0.05, 0.05)
)

# --- 2.4  Joining (merging) datasets -----------------------------------------

data_left  <- left_join(data1, data2)   # all data1 rows; Weight added where matched
data_right <- right_join(data1, data2)  # all data2 rows; Nature added where matched
data_inner <- inner_join(data1, data2)  # only rows present in both
data_full  <- full_join(data1, data2)   # all rows from both
data_semi  <- semi_join(data1, data2)   # data1 rows that have a match (no new columns)
data_anti  <- anti_join(data1, data2)   # data1 rows with no match in data2

# --- 2.5  Binding (pooling) datasets -----------------------------------------

f2010 <- read_xls("data/f2010.xls")
f2012 <- read_xls("data/f2012.xls")
f2014 <- read_xls("data/f2014.xls")

dim(f2010)
slice(f2010, 1:10)
str(f2010)

library(dplyr)
fpanel <- bind_rows(f2010, f2012, f2014)

dim(fpanel)
head(fpanel, 10)
str(fpanel)


# =============================================================================
# CHAPTER 3 — Exploring the Crimes Dataset
# =============================================================================

crimes <- read_xls("data/crimes.xls", sheet = "data", col_names = TRUE)

dim(crimes)
slice(crimes, 1:10)
names(crimes)
str(crimes)

# --- Descriptive statistics for a single variable ----------------------------

min(crimes$poverty_index,    na.rm = TRUE)
max(crimes$poverty_index,    na.rm = TRUE)
mean(crimes$poverty_index,   na.rm = TRUE)
sd(crimes$poverty_index,     na.rm = TRUE)
median(crimes$poverty_index, na.rm = TRUE)
quantile(crimes$poverty_index, na.rm = TRUE)

# --- Summary statistics for quantitative variables ---------------------------

crimes_quant <- subset(crimes,
                       select = c(crimes, gdp_2011, poverty_index,
                                  population, unemp_rate))
summary(crimes_quant)

# --- Frequency table for the categorical variable big_region -----------------

table(crimes$big_region)

library(questionr)
freq(crimes$big_region, cum = TRUE, sort = TRUE, valid = FALSE,
     total = TRUE, na.last = TRUE)

# --- Identify extreme observations -------------------------------------------

subset(crimes, crimes == min(crimes$crimes), select = c(dep_name, crimes))
subset(crimes, crimes == max(crimes$crimes), select = c(dep_name, crimes))

# --- Correlation analysis -----------------------------------------------------

cor(crimes$crimes, crimes$poverty_index, method = "pearson")

cor.test(crimes$crimes, crimes$poverty_index,
         alternative = "two.sided", method = "pearson", conf.level = 0.95)

cor.test(crimes$crimes, crimes$poverty_index,
         alternative = "two.sided", method = "spearman", conf.level = 0.95)

# --- Create binary (dummy) variables -----------------------------------------

crimes$crim_degree <- ifelse(crimes$crimes > mean(crimes$crimes), 1, 0)
crimes$poverty_cat <- ifelse(crimes$poverty_index > mean(crimes$poverty_index), 1, 0)

# --- Cross-tabulation and chi-square test ------------------------------------

library(gmodels)

CrossTable(crimes$crim_degree, crimes$poverty_cat,
           digits = 2, prop.chisq = FALSE, chisq = FALSE,
           fisher = FALSE, mcnemar = FALSE,
           missing.include = FALSE, format = "SPSS")

CrossTable(crimes$crim_degree, crimes$poverty_cat,
           digits = 2, prop.chisq = FALSE, chisq = TRUE,
           fisher = FALSE, mcnemar = FALSE,
           missing.include = FALSE, format = "SPSS")


# =============================================================================
# CHAPTER 4 — Data Visualization
# =============================================================================

# --- 4.1  Base R plots -------------------------------------------------------

plot(crimes$poverty_index,
     main = "Poverty Index by Department",
     xlab = "Department Index",
     ylab = "Poverty Index Value",
     type = "p", pch = 20, cex = 1, col = "blue")

# Add department labels — note: text() is called separately, not with +
plot(crimes$poverty_index,
     main = "Poverty Index by Department",
     xlab = "Department Index",
     ylab = "Poverty Index Value",
     type = "p")
text(crimes$poverty_index,
     labels = crimes$dep_name,
     cex = 0.7, pos = 1)

# --- 4.2  Histogram and kernel density ---------------------------------------

hist(crimes$poverty_index,
     main = "Distribution of the Poverty Index",
     xlab = "Poverty Index")

dens <- density(crimes$poverty_index, bw = "nrd0",
                kernel = "epanechnikov", na.rm = TRUE)
plot(dens, frame = TRUE, col = "steelblue",
     main = "Kernel Density — Poverty Index")

# --- 4.3  Scatter plots with ggplot2 -----------------------------------------

library(ggplot2)

ggplot(crimes, aes(x = poverty_index, y = crimes)) +
  geom_point() +
  labs(title    = "Scatterplot",
       subtitle = "Relationship between Crimes and Poverty",
       x        = "Poverty Index",
       y        = "Number of Crimes",
       caption  = "Source: DCPJ, processed by ONDRP.")

ggplot(crimes, aes(x = poverty_index, y = crimes)) +
  geom_point() +
  geom_smooth(method = lm, se = TRUE) +
  labs(title    = "Scatterplot with OLS Fit",
       subtitle = "Relationship between Crimes and Poverty",
       x        = "Poverty Index",
       y        = "Number of Crimes",
       caption  = "Source: DCPJ, processed by ONDRP.")

ggplot(crimes, aes(x = poverty_index, y = crimes)) +
  geom_point() +
  geom_smooth(method = loess, se = TRUE) +
  labs(title    = "Scatterplot with LOESS Fit",
       subtitle = "Relationship between Crimes and Poverty",
       x        = "Poverty Index",
       y        = "Number of Crimes",
       caption  = "Source: DCPJ, processed by ONDRP.")

# --- 4.4  Correlation matrix and scatter-plot matrix -------------------------

# Working dataset — quantitative variables only
crim_tab <- subset(crimes,
                   select = c(crimes, gdp_2011, poverty_index,
                              unemp_rate, population, big_region))

summary(crim_tab, digits = 2)

# Numeric subset for correlation and base plot (big_region is non-numeric)
crim_num <- select(crim_tab, where(is.numeric))

plot(crim_num)

round(cor(crim_num, method = "pearson"), 2)

library(GGally)
ggpairs(crim_tab, columns = 1:5, aes(colour = big_region, alpha = 0.5))


# =============================================================================
# CHAPTER 5 — Ordinary Least Squares Regression
# =============================================================================

MyModel <- lm(crimes ~ gdp_2011 + poverty_index + unemp_rate + population,
              data = crimes)

names(MyModel)

options(scipen = 999)
summary(MyModel)

confint(MyModel, level = 0.95)

crimes_hat <- fitted(MyModel)
head(MyModel$residuals)

# --- Multicollinearity (VIF) --------------------------------------------------

library(car)

vif(MyModel)
sqrt(vif(MyModel)) > 2   # rule of thumb: VIF > 4 warrants investigation

# --- Autocorrelation (Durbin-Watson) -----------------------------------------

durbinWatsonTest(MyModel)

# --- Heteroskedasticity (non-constant variance test) -------------------------
# H0: constant error variance   H1: variance changes with fitted values

ncvTest(MyModel)

# --- Residual diagnostic plots -----------------------------------------------

spreadLevelPlot(MyModel)
qqPlot(MyModel, main = "QQ Plot — Studentised Residuals")
leveragePlots(MyModel)
outlierTest(MyModel)


# =============================================================================
# CHAPTER 6 — Log-Linear Models
# =============================================================================

# Create per-capita log-transformed variables
crimes$lcrimes_cap <- log(crimes$crimes / crimes$population)
crimes$lgdp_cap    <- log(crimes$gdp_2011 * 1e6 / crimes$population)
crimes$lpop        <- log(crimes$population)

# Log-log model
MyModel_log <- lm(lcrimes_cap ~ lgdp_cap + poverty_index + unemp_rate + lpop,
                  data = crimes)
summary(MyModel_log)

# Model comparison (use AIC/BIC — R² cannot be compared across different dep. vars.)
library(broom)
glance(MyModel)
glance(MyModel_log)


# =============================================================================
# CHAPTER 7 — Instrumental Variables Estimation
# =============================================================================

library(AER)

# --- Stage 1: reduced-form equation ------------------------------------------

reducedForm <- lm(lgdp_cap ~ pupils + unemp_rate + poverty_index + lpop,
                  data = crimes)
summary(reducedForm)

coeftest(reducedForm, vcov = vcovHC, type = "HC1")
summary(reducedForm)$r.squared

lgdp_pred <- fitted(reducedForm)

# --- Stage 2: structural equation --------------------------------------------

structuralEq <- lm(lcrimes_cap ~ lgdp_pred + unemp_rate + poverty_index + lpop,
                   data = crimes)
summary(structuralEq)
coeftest(structuralEq, vcov = vcovHC, type = "HC1")

# --- Hausman endogeneity test (regression-based) -----------------------------
# Add Stage-1 residuals to the structural equation.
# A significant coefficient on gdp_residuals confirms endogeneity.

HausmanTest  <- lm(lgdp_cap ~ pupils + unemp_rate + poverty_index + lpop,
                   data = crimes)
gdp_residuals <- residuals(HausmanTest)   # residuals, not fitted values

endoTest <- lm(lcrimes_cap ~ lgdp_cap + unemp_rate + poverty_index + lpop +
                 gdp_residuals, data = crimes)
coeftest(endoTest, vcov = vcovHC, type = "HC1")

# --- IV estimation with ivreg() ----------------------------------------------
# Formula: structural regressors | exogenous variables + external instrument

ivreg_est <- ivreg(lcrimes_cap ~ lgdp_cap + unemp_rate + poverty_index + lpop |
                                 pupils   + unemp_rate + poverty_index + lpop,
                   data = crimes)

summary(ivreg_est, diagnostics = TRUE)
coeftest(ivreg_est, vcov = vcovHC, type = "HC1")


# =============================================================================
# CHAPTER 8 — Panel Data Models
# =============================================================================

library(wooldridge)
library(plm)

data(wagepan)

# Fixed effects (within estimator — removes time-invariant individual effects)
panel.fe <- plm(lwage ~ educ + exper + expersq + union + south + married + black,
                data  = wagepan,
                index = c("nr", "year"),
                model = "within")
summary(panel.fe)

# Random effects
panel.re <- plm(lwage ~ educ + exper + expersq + union + south + married + black,
                data  = wagepan,
                index = c("nr", "year"),
                model = "random")
summary(panel.re)

# Hausman test: H0 = RE consistent and efficient; H1 = FE preferred
phtest(panel.fe, panel.re)

# Formatted results table
library(stargazer)
stargazer(MyModel, MyModel_log, structuralEq,
          title = "Regression Results",
          align = TRUE,
          type  = "text")   # change to type="latex" for LaTeX output


# =============================================================================
# CHAPTER 9 — Binary Choice Models: Logit and Probit
# =============================================================================

library(wooldridge)

data(mroz)
str(mroz)
head(mroz)

# Share of women in the labor force
mean(mroz$inlf)

summary(mroz[, c("inlf", "nwifeinc", "educ", "exper", "age", "kidslt6", "kidsge6")])

# --- 9.1  Linear Probability Model (baseline OLS) ----------------------------

lpm <- lm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
          data = mroz)
summary(lpm)

# Check for out-of-range predictions
pred_lpm <- fitted(lpm)
sum(pred_lpm < 0 | pred_lpm > 1)

# --- 9.2  Logit model --------------------------------------------------------

logit_model <- glm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
                   data   = mroz,
                   family = binomial(link = "logit"))
summary(logit_model)

# --- 9.3  Probit model -------------------------------------------------------

probit_model <- glm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
                    data   = mroz,
                    family = binomial(link = "probit"))
summary(probit_model)

# --- 9.4  Average Partial Effects (marginal effects) -------------------------

library(margins)

ape_logit  <- margins(logit_model)
summary(ape_logit)

ape_probit <- margins(probit_model)
summary(ape_probit)

# Manual APE for education (logit) — illustrates the formula
p_hat      <- fitted(logit_model)
f_logit    <- p_hat * (1 - p_hat)           # logistic density at each observation
beta_educ  <- coef(logit_model)["educ"]
ape_educ   <- mean(f_logit * beta_educ)
ape_educ

# --- 9.5  Model evaluation ---------------------------------------------------

# McFadden pseudo-R² and other pseudo-R² measures
library(pscl)
pR2(logit_model)
pR2(probit_model)

# AIC and BIC
AIC(logit_model, probit_model)
BIC(logit_model, probit_model)

# Confusion matrix (threshold = 0.5)
class_logit  <- ifelse(fitted(logit_model)  > 0.5, 1, 0)
class_probit <- ifelse(fitted(probit_model) > 0.5, 1, 0)

table(Predicted = class_logit,  Actual = mroz$inlf)
table(Predicted = class_probit, Actual = mroz$inlf)

# Percentage correctly classified
mean(class_logit  == mroz$inlf)
mean(class_probit == mroz$inlf)

# --- 9.6  Odds ratios (logit only) -------------------------------------------

exp(coef(logit_model))
exp(confint(logit_model))

# --- 9.7  Comparison table: LPM, Logit, Probit -------------------------------

library(stargazer)
stargazer(lpm, logit_model, probit_model,
          title          = "Labor Force Participation — Binary Choice Models",
          dep.var.labels = "In Labor Force (1 = yes)",
          type           = "text")
