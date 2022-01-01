
# Check the number of installed packages.

options(repos="https://cran.rstudio.com")
nrow(installed.packages())

# The number of current packages available at the CRAN website.

nrow(available.packages())

# Installing and loading "tidyverse" package.

# install.packages("tidyverse")  # Remove the "#" symbol to enable the installation.

# Loading "tidyverse" package into the active memory.

library("tidyverse")

# Getting information on the "tidyverse" package.

library(help="tidyverse")


# Check the content of "tidyverse" package.

tidyverse_packages()

# Getting help on the "mean" function.

help(mean)

# Another way to ask for help. 
? mean

# We can use R as a simple calculator.

2*2
100/5
2^10

# Store the result in the object "a" using the R's assignement symbol <-

a <- 2^10

# Apply the log function to the object "a".

log(a)

# Creating variables.

countries <- c("France", "Germany", "Spain", "Italy")
CovidCases <- c(1502763, 580415, 1331756, 759829)
CovidDeaths <- c(38289,10904,36495,39412)

# Creating 'covidata' set as a data frame.

covidata <- data.frame(countries,CovidCases,CovidDeaths)

# Saving a copyt of 'covidata' in the local working directory.

save(covidata,file="covidata.RData")

# Deleting the "covidata" file from the local directory.

unlink("covidata.RData")

# Checking the class category of 'covidata'.

class(covidata)

# Looking at the structure of 'covidata'.

str(covidata)


# Showing the 'covidata' content.

covidata

# Visualising data in a sheet. 

View (covidata)

# Showing the names of all the columns (variables) of "covidata" which are three in this case.

names(covidata)

# Loading the "datasets" package (installed with R).

library("datasets")

# Getting information on the "datasets" package.

library(help="datasets")

# Showing the list of available internal datasets.

data()

# Loading the "mtcars" data.

data(mtcars)

# First, to have more information on the mtcars data columns (variables), we can ask for help on this dataset.

?mtcars

# Sructure of mtcars

str(mtcars)


# Show the first 10 rows of the data. The default is 6 rows.

head(mtcars, 10)


# Show the last rows of data.

tail(mtcars)

# To show how many columns the data contains.

ncol(mtcars)

# To show the number of rows in the dataset.

nrow(mtcars)

# Read the "crimes_SC.csv" dataset where every single data is separated by a Semicolon.

crimes1 <- read.table(file="crimes_SC.csv", sep=";", header=TRUE, dec=".")


# Read the "crimes_C.csv" dataset where data is separated by a Comma.

crimes2 <- read.table(file="crimes_C.csv", sep=",", header=TRUE, dec=".")


# Read the "crimes_SC.csv" dataset where data is separated by a Tabulation.

crimes3 <- read.table(file="crimes_T.csv", sep="\t", header=TRUE, dec=".", quote = "\"")

#Load the "readxl" package.

library("readxl")

#Getting help
help(package="readxl")

# Read into R the "crimes.xls" data file.

crimes_xls <- read_xls("crimes.xls", sheet=("data"))


# Creating data by insertion: creating the first dataset "data1" as a tibble format.

# This first data file (data1) contains three variables: Animals, their babies, and their nature as wild or domestic.

data1 <-tibble(Animals=c("Cat","Chicken","Dog","Elephant","Giraffe","Gnu","Lion",
                         "Panda","Penguin","Rabbit"),
Baby=c("Kitten","Chick","Puppy","Calf","Calf","Calf","Cub","Cub","Chick","Kitten"),
Nature=c("Domestic","Domestic","Domestic","Wild","Wild","Domestic","Wild","Wild",
         "Wild","Domestic"))


# Creating the second data file "data2" as a tibble format.

# The second data file (data2) contains also three variables: Animals, their babies, and a new variable, on the 'weight of babies' at birth measured in a kilogram. We notice that the variable 'Animals' includes now two new animals: 'Turkey' and 'Duck'. Likewise, for their babies and weight at birth. 

data2<- tibble(Animals=c("Cat","Chicken","Dog","Elephant","Giraffe","Gnu",
                         "Lion","Panda","Penguin","Rabbit","Turkey","Duck"),
Baby=c("Kitten","Chick","Puppy","Calf","Calf","Calf","Cub","Cub","Chick",
       "Kitten","Poult","Duckling"), 
Weight=c(0.08,0.9,0.45,90,99,36,8,0.1,0.32,0.03,0.05,0.05))


# Joining (merging) data files.

data_left<- left_join(data1, data2)

data_right<- right_join(data1, data2)

data_inner <- inner_join(data1, data2)

data_full <- full_join(data1, data2)

data_semi <- semi_join(data1,data2)

data_anti <- anti_join(data1,data2)

# The following code read three Excel files (f2010.xls, f2012.xls and f2014.xls) and assign them to three data frames: f2010, f2012 and f2014.

library("readxl")
f2010 <- read_xls("f2010.xls")
f2012 <- read_xls("f2012.xls")
f2014 <- read_xls("f2014.xls")

# Exploring The f2010 data dimension (how many rows and columns contains this data).

dim(f2010)

# Showing the first ten rows and all columns. 

slice(f2010,1:10)


# Examining the content of f2010 data (see the columns name).

str(f2010)



# pooling together the three files containing the same variables observed for three different years, f2010, f2010, f2014 using 'blind_rows()' function.

library("dplyr")
fpanel <- bind_rows(f2010,f2012,f2014)

# Exploring the dimension of the entire data after pooling the three files.

dim(fpanel)

# Showing the first ten rows using the "head()" function instead.

head(fpanel,10)

# Exploring the content of the data which equals now: 10 same columns (or variables) and 321 rows (or observations). 

str(fpanel)


# Loading in memory the "readxl" library.

library("readxl")


# Reading "crimes.xls" file using "read_xls" function.

crimes <- read_xls("crimes.xls", sheet="data", col_names=TRUE)


# The dimension of the crimes data (rows x columns).

dim (crimes)



# Showing the first ten rows of the data using "slice" function from "dpyr" package. The "head(crimes, 10)" function from the R Base package could be also used.  

slice(crimes,1:10)

# Getting the columns names of "crimes" data using the "names()" function. 

names(crimes) 


# Data structure of "crimes" object.

str(crimes)

# Computing the Minimum of poverty_index.

min(crimes$poverty_index, na.rm=TRUE)

# Computing the Maximum of poverty_index.

max(crimes$poverty_index, na.rm=TRUE)

# Computing the Mean of poverty_index.

mean(crimes$poverty_index, na.rm=TRUE)

# Computing the Sdandard Deviation of poverty_index.

sd(crimes$poverty_index, na.rm=TRUE)

# Computing the Median of poverty_index.

median(crimes$poverty_index, na.rm=TRUE)

# Computing the quantile of poverty_index (Min, Q1, Q2, Q3, Max).

quantile(crimes$poverty_index)

# Creating a data for quantitative variables from "crimes" using "subset()" and "select()" functions.

crimes_quant <- subset(crimes, select= c(crimes, gdp_2011, poverty_index, population, unemp_rate))

# Summary statistics 

summary(crimes_quant)

# Installing and loading the "tufte" package.

install.packages("tufte")
library(tufte)

# A simple plot of poverty index

plot(crimes$poverty_index, main="Poverty Index by department", xlab="Department Index", 
     ylab="Poverty index value", type="p", pch=20, cex=1, col="blue")

# Plot of poverty index using "department name" as labels on the observations points.

plot(crimes$poverty_index, main="Poverty index by department", xlab="Department index", 
     ylab="Poverty index value", type="p") +
text(crimes$poverty_index, labels=crimes$dep_name, cex=0.7, pos=1)

# Plotting of poverty index distribution using the Histogram plot   

hist(crimes$poverty_index, main="Poverty distribution using Histogram", 
     xlab="Poverty index")

# Computing the density (Kernel Density Estimation) 

dens <- density(crimes$poverty_index, bw="nrd0", kernel="epanechnikov", na.rm = TRUE)

# Plotting of poverty distribution using the calculated density

plot(dens, frame = TRUE, col = "steelblue", cex=1.5,
     main = "Poverty distribution using Density estimation")

# Counting the number of departments in each region. The variable "big_region" is computed by aggregating 22 French regions (following the division of the regions before 2016) in six big regions (1. North-East, 2. North-West, 3. Ile-de-France, 4.Center, 5.South-East, 6. South-West)

# We use the "table()" function for counting.
table(crimes$big_region)

# install.packages("questionr") # remove the # symbol at the beginning to install the # package

library("questionr") 

freq(crimes$big_region, cum=TRUE, sort=TRUE, valid=FALSE, total=TRUE, na.last=TRUE)

# Scatter plot between "crimes" and "poverty_index"

plot(x=crimes$poverty_index, y=crimes$crimes, 
     main="Relashioship between Crimes and Poverty", xlab="Poverty Index", 
     ylab="Crimes", type="p", pch=20, cex=1, col="blue")

# Computing the correlation coefficient between "crimes" and "poverty_index' using the Pearson method. 

cor(crimes$crimes, crimes$poverty_index, method="pearson")

# Performing the statistical significance test of the correlation between "crimes" and "poverty_index" using "pearson" method.

cor.test(crimes$crimes, crimes$poverty_index, alternative="two.sided",
         method="pearson", conf.level=0.95)

# Performing the test for statistical significance of the correlation between "crimes" and "poverty_index" using "spearman" method, which is a non-parametric one. This method is suited when the supposed link between variables is not linear.

cor.test(crimes$crimes, crimes$poverty_index, alternative="two.sided", 
         method="spearman", conf.level=0.95)

# Creating a dummy for "crimes", named 'crim_degree' based on the mean threshold of "crimes".

crimes$crim_degree <- ifelse(crimes$crimes>mean(crimes$crimes), 1,0 )

# Creating a dummy for "poverty_index", named 'poverty_cat' based on the mean threshold of "poverty_index".

crimes$poverty_cat <- ifelse(crimes$poverty_index>mean(crimes$poverty_index),1,0)

# Performing cross-tabulation between "crim_degree" and "poverty_cat" using "CrossTable()" function.

# Insalling the package "gmodels" to get "CrossTable()" function.
#install.packages("gmodels")  # Remove # symbol to install.

# Loading the package "gmodels" in memory.
library("gmodels")
 
CrossTable(crimes$crim_degree, crimes$poverty_cat, digits=2, 
          prop.chisq=FALSE, chisq=FALSE,fisher=FALSE,mcnemar=FALSE, 
          missing.include = FALSE, format="SPSS")

# Testing the independence between the crimes levels ("cri_degree") and poverty intensity (poverty_cat) expressed as dummy variables using Chi-square test.

CrossTable(crimes$crim_degree,crimes$poverty_cat, digits=2, prop.chisq=FALSE, chisq=TRUE,fisher=FALSE,mcnemar=FALSE, missing.include = FALSE,format="SPSS")

# Scatter plot of "crimes" and "poverty_index" using 'ggplot2'.
# First, we load 'ggplot2' in memory using 'library()'.

library("ggplot2")

ggplot(crimes, aes(x=poverty_index, y=crimes)) +
  geom_point() +   # Show dots 
  labs(subtitle="The link between Crimes and Poverty", 
       x="Poverty index", 
       y="Crimes", 
       title="Scatterplot", 
     caption = "Source: DCPJ and data processing by ONDRP.") 


# Adding a linear adjustment (regression line) with confidence interval around the line.

ggplot(crimes, aes(x=poverty_index, y=crimes)) +
  geom_point() +   # Show dots 
  labs(subtitle="The link between Crimes and Poverty", 
       x="Poverty index", 
       y="Crimes", 
       title="Scatterplot", 
     caption = "Source: DCPJ and data processing by ONDRP.") +
  geom_smooth(method=lm)


# Adding a linear adjustment (regression line) without the confidence interval draw.

ggplot(crimes, aes(x=poverty_index, y=crimes)) +
  geom_point() +   # Show dots 
  labs(subtitle="The link between Crimes and Poverty", 
       x="Poverty index", 
       y="Crimes", 
       title="Scatterplot", 
     caption = "Source: DCPJ and data processing by ONDRP.") +
  geom_smooth(method=lm, se=FALSE)

# Adding a non-linear adjustment using  the LOcally Estimated Scatter plot Smoothing (loess) and showing confidence interval draw.

ggplot(crimes, aes(x=poverty_index, y=crimes)) +
  geom_point() +   # Show dots 
  labs(subtitle="The link between Crimes and Poverty", 
       x="Poverty index", 
       y="Crimes", 
       title="Scatterplot", 
     caption = "Source: DCPJ and data processing by ONDRP.") +
  geom_smooth(method=loess, se=TRUE)

# Creating a temporary data "crim_tab" containing the main variables, using the functions "subset()" and "select()".

crim_tab <- subset(crimes, select=c(crimes, gdp_2011,poverty_index,unemp_rate,
                  population, big_region) )

# Checking the columns names in "crim_tab".
names(crim_tab)

# Showing the variables type using "str()" command.
str(crim_tab)

# Taking the database "crim_tab" as the main dataset reference. 
attach(crim_tab)

# Performing descriptive statistics for all variables in one step.

summary(crim_tab, digits=2)

# Showing which department has the minimum crime number in the sample.
# It should be noted that the database used here is "crimes" which contains all variables.
 
subset(crimes, crimes==min(crimes), select=c(dep_name, crimes))

# Showing which department has the maximum crime number in the sample.

subset(crimes, crimes==max(crimes), select=c(dep_name, crimes))

# Scatter plots of the variables of interest. 

corr_plot <- plot(crim_tab)


# Correlation coefficients matrix.

cor(crim_tab, method="pearson")

# Correlation coefficients matrix, while limiting the decimals' number to 2 digits using the _round()_ function.

round(cor(crim_tab, method="pearson"), 2)


# Installing "GGally{}" package.

# install.packages("GGally")  # Remove the # symbol to start installing.

# Load in memory the "GGally{}" package

library("GGally")

# Performing scatter plots and correlation coefficients matrix in one step using "ggpairs()" function.

ggpairs(crim_tab, columns =1:5) 

# Estimating regression coefficients using "lm()" function.

MyModel <- lm(crimes ~ gdp_2011 + poverty_index + unemp_rate + population, data=crimes)


# Checking the cotent of "MyModel" object.

names(MyModel)

# Summarising the regression results.
# The option "(scipen=999)" is added to avoid the output of numbers with scientific notation.

options(scipen=999)
summary(MyModel)

# Getting the 95% confidence intervals of the estimated coefficients.

confint(MyModel, level=0.95)

# Computing the predicted values of y (crimes).

crimes_hat <- fitted(MyModel)

# Viewing the residuals (estimation of model errors).

View(MyModel$residuals)

# Installing the "car" package made by Fox and S. Weisberg as a companion for applied regression and diagnostics.
# https://cran.r-project.org/web/packages/car/index.html

# install.packages("car")    # Remove the # symbol to start installing.

# Loading "car" package.

library("car")

# Compute the Variance Inflation Factor (VIF)

vif(MyModel)

# Logic test of the presence of Variance Inflation problem. We may have multi-collinearity problem if the VIF statistic is superior to 2.  

sqrt(vif(MyModel)) > 2

# Test for Autocorrelated errors 

durbinWatsonTest(MyModel)

# Computing non-constant error variance test.

# Null hypothesis H0:  Constant error variance.
# Alternative hypothesis H1: The error variance changes with the level of the response (fitted values), or with a linear combination of predictors.

ncvTest(MyModel)


# Plotting residuals versus fitted values
spreadLevelPlot(MyModel)


# Assessing outliers. Compute the Bonferonni p-value for most extreme observations.
outlierTest(MyModel)

# Assessing outliers. Draw the qq plot for studentised residuals.

qqPlot(MyModel, main="QQ Plot")

# Assessing the outliers and influential observations.

leveragePlots(MyModel)

# Compute new log scale variables

crimes$lcrimes_cap <- log(crimes$crimes/crimes$population)
crimes$lgdp_cap <- log(crimes$gdp_2011*1000000/crimes$population)
crimes$lpop <- log(crimes$population)

# Estimating regression coefficients using "lm()" function.

MyModel_log <- lm(lcrimes_cap ~ lgdp_cap + poverty_index + unemp_rate + lpop, data=crimes)

# Showing results of the new estimation

summary(MyModel_log)

# Install "broom" package
# "broom" package is a component of the "tidyverse". If "tidyverse" is 
# already installed, we don't need to install again "broom". We only load it.
# Otherwise we must install it.

# Load "broom" package
library(broom)

# Using 'glance' function to compare OLS estimates

# Statistics for the model with the original non-transformed variables
glance(MyModel)

# Statistics for the model with logarithmic scale variables
glance(MyModel_log)

# Install the "AER" package: Kleiber, Ch., Zeileis, A. (2008). Applied Econometrics with R, Springer-Verlag

# install.packages("AER")

# Load "AER" , "lm" packages

library(AER)

# 1. Perform the first stage regression: reduced form equation
reducedForm <- lm(lgdp_cap ~ pupils + unemp_rate + poverty_index + lpop, data = crimes)

# Inspecting the regression results

summary(reducedForm)

# Statistical significance test of pupils variable
# H0 assumption: the beta(pupils) coefficient is equals to zero

coeftest(reducedForm, vcov = vcovHC, type = "HC1")

# inspect the R^2 of the first stage regression

summary(reducedForm)$r.squared

# store the predicted values to be used in the second stage regression as 
# independent variable, istead of the original variable 'pupils' 


lgdp_pred <- reducedForm$fitted.values

# A second way to predict the fitted values 
lgdp_pred2 <- predict(reducedForm)

# Run the second stage regression: structural equation
structuralEq <- lm(lcrimes_cap ~ lgdp_pred +  unemp_rate + poverty_index + lpop, data = crimes)

# Test of statistical significance of pupils variable
coeftest(structuralEq, vcov = vcovHC)


# Show the results of the structural equation
summary(structuralEq)

# Hausman procedure to test endogeneity assumption Regression 

HausmanTest <- lm(lgdp_cap ~ pupils +  unemp_rate + poverty_index + lpop, data = crimes)

# Predict residuals
gdp_residuals <- predict(HausmanTest)

# Test the significance of residuals term via a regression within the structural equation

endoTest <- lm(lcrimes_cap ~ lgdp_cap + unemp_rate + poverty_index + lpop +
                 gdp_residuals, data = crimes)

# Residuals term significance test

coeftest(endoTest, vcov = vcovHC)


# IV estimator using 'ivreg' function

ivreg_est <- ivreg(lcrimes_cap ~  lgdp_cap + unemp_rate + poverty_index + lpop | 
                                             unemp_rate + poverty_index + lpop + pupils, data = crimes)

# Schow the estimates results of using the 'ivreg'  function.

summary(ivreg_est)

# The linear coefficients test hypothesis

coeftest(ivreg_est, vcov = vcovHC, type = "HC1")

# Load the "Wooldridge" package
library(wooldridge)

# Load the 'plm' package on 'Panel data econometrics in R', by Yves Croissant and Giovanni Millo.

library(plm)

# Create the wage dataset

wage <- data(wagepan)

# Fixed effects panel estimate 

panel.fe <- plm(lwage ~ educ + exper + expersq + union + south + married +
                  black, data=wagepan, index=c("nr", "year"), model="within" )

# Results summary of fe
summary(panel.fe)

# Random effects panel estimate 
panel.re <- plm(lwage ~ educ + exper + expersq + union + south + married +
                  black, data=wagepan, index=c("nr", "year"), model="random" )

# Results summary of re 
summary(panel.re)

# Hausman specification test, to choose between models (Fe vs. Re)
# H0: There is no difference between the two models

phtest(panel.fe, panel.re)


# Insgtall and load the "stargazer" package install.packages("stargazer")

library(stargazer)

# Format regression results as a LateX table

stargazer(MyModel, MyModel_log, structuralEq, title = "Regression Results", align = TRUE)












