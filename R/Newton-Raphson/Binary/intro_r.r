
################################################
# Tips

# Use google to find names of functions
# Use ? for help documentation once you find the name of the function 
# that you need, for example:
?mean

# Use scripts instead of the command line
# Use control + enter to send code from a script to the command line
# Comment your code
# Don't use attach()

################################################
# The topics that we covered during the first three workshops:

# data types and variables
# data frames
# univariate summaries
# comparing 2 variables
# plots
# writing functions


################################################
# data types and variables

# Numbers, characters, and strings are entered as follows:
34
"a"
"abc"

# There are two assignment operators that are used to assign a value to a variable.
# They are basically equivalent:
x = "abc"
y <- "def"

# R has functions for all of the usual arithmetic operations, for example:
4*5
4 - 5
4 + 5
4^5
sqrt()
exp()
log()

# Here are four different ways to create vectors in R:
c(1, 2, 3)
1:20
rep(4, 10)
a = seq(2, 15, .5)

# Specific elements in a vector can be called by using the bracket notation:
a[2]
a[c(2, 3, 4, 10)]

# There are many ways to make matrices in R.  One way is to start with vectors,
# and then combine them with the row-bind or column-bind functions:
a = c(1, 2, 3)
b = c(3, 4, 5)
c = c(6, 7, 8)
rbind(a, b, c)
x = cbind(a, b, c)

# Entries, rows, and columns of a matrix can be called using bracket notation:
x[2, 3]
x[2, ]
x[ , 3]
# Note that R does not distinguish between row and column vectors.

################################################
# data frames

# To load data in R, we first need to tell R which directory to use.
# We'll use quotation marks, since R will see the path as a string.
# On windows computers, remember to change all the backslashes to forward-slashes.
# Alternatively, as Rolondo described, you can use 
# a double backslash instead of a forward-slash.
setwd("/home/josh/Desktop/Dropbox/Summer 2015/")

# If you want to verify that you're in the correct directory, you can list
# all the files that are available:
list.files()
# Or you can ask R which directory it's using:
getwd()

# To load a csv file as a data frame:
df = read.table("practice_data_1.csv", header = T)

# Here we're loading a second data frame, and then creating a third
# data frame by merging.
df2 = read.table("practice_data_2.csv", header = T)
data = merge(df, df2, by = "ID")

# We can save the merged data frame to a csv file:
write.table(data, file = "merged_data.csv", sep = " ", row.names = F)
# If you want to save your data in a different location, just do a setwd()
# before you do write.table().

# Finally, clean up the workspace by removing data frames that we will no longer be using:
rm(df, df2)

# There are three ways to call a variable from a data frame.  You'll probably
# find this one to be the most useful:
data$weight

# A data frame is like a matrix, so you can also use matrix notation to ask for
# the corresponding column vector:
data[ , 4]

# The double-bracket notation is very useful when you write functions or 
# write scripts to automate data cleaning:
data[["weight"]]

# Here is the formula for body mass index:
# BMI = 703 * weight/ height^2
# Here we're creating a new variable in our data frame for BMI.
# It's not a bad idea to initialize your variable with an empty string of NA's,
# (though this is often not necessary):
data$BMI = NA
data$BMI = 703* data$weight / (data$height^2)
# Remember that RStudio does not automatically refresh the view of your data frame.
# To verify that a new variable has been created, you can click on the name
# of the data frame in the global environment window.

# As Jennifer pointd out, the definitions for BMI categories that I gave in class were
# not so good.  Here's a better definition:
# underweight: BMI < 18.5
# normal:      18.5 <= BMI < 25.0
# overweight:  25.0 <= BMI < 30.0
# obese:       BMI >= 30.0

# First let's create an indicator variable for underweight status.
# The ifelse() function is very useful:
data$BMI_under = NA
data$BMI_under = ifelse(data$BMI < 18.5, 1, 0)

# Next let's create a categorical variable with levels for each of the
# four BMI categories: Under, Normal, Over, Obese.
data$BMI_cat = NA
data$BMI_cat = ifelse(data$BMI < 18.5, "Under", data$BMI_cat)
data$BMI_cat = ifelse(data$BMI >= 18.5 & data$BMI < 25, "Normal", data$BMI_cat)
data$BMI_cat = ifelse(data$BMI >= 25   & data$BMI < 30, "Over", data$BMI_cat)
data$BMI_cat = ifelse(data$BMI >= 30, "Obese", data$BMI_cat)


################################################
# univariate summaries

# R has a summary function that provides some useful information:
summary(data)

# numerical variables
# Individual summary statistics for numerical variables are easy to obtain:
max(data$weight)
min(data$weight)
mean(data$weight)
sd(data$weight)

# You can make a histogram for a graphical summary:
hist(data$weight, breaks = 20)

# Note that the hist() function contains a bug that occurs when
# the variable has integer values.  There are only 10 bins
# in this histogram (there should be 11), and the left-most bin is too high:
hist(data$age)
# To test for this bug, re-run the histogram with something like 0.1 added to
# each value, and see if you get a different plot:
hist(data$age + .1)
# To fix this bug, specify exactly where you want your bin breaks:
hist(data$age, breaks = seq(34.5, 45.5, 1))


# categorical variables
# The summary() function will provide some information about categorical variables.
# You can determine the levels of a categorical function:
unique(data$sex)
# And count how many observations you have for each level:
length(which(data$sex == "M"))
length(which(data$sex == "F"))
# And also find out the corresponding proportions:
length(which(data$sex == "M")) / nrow(data)
length(which(data$sex == "F")) / nrow(data)
# There are 79 males and 71 females.  You can make a barplot as follows:
s = c(79, 71)
names(s) = c("M", "F")
barplot(s)

# summarizing with missing data
# Missing values are indicated with "NA".  You can find out how many observations
# are missing for a specific variable:
length(which(is.na(data$Assay)))
# You can obtain numerical summaries in the presence of missing data 
# by adding the option "na.rm = T".  
mean(data$Assay, na.rm = T)
sd(data$Assay, na.rm = T)
# Remember that estimates obtained in this way may be biased.


################################################
# Comparing 2 variables

# numerical vs. numerical
# You can get a correlation coefficient for two variables:
cor(data$height, data$weight)
# A scatterplot is a useful graphical summary:
plot(data$height, data$weight)

# numerical vs. categorical
# First we can obtain numerical summaries for each level of the categorical
# variable.  For example, mean heights for males and females:
mean(data$height[which(data$sex == "M")])
mean(data$height[which(data$sex == "F")])
# You can use boxplots for graphical summaries:
boxplot(height ~ sex, data = data)
boxplot(weight ~ Outcome_2, data = data)
# It's often necessary to reorder factor levels.  The following command sets
# the sequence of factor levels without changing the data:
data$Outcome_2 = factor(data$Outcome_2, levels = c("Lo", "Mid", "Hi"))
# Redraw the boxplot to see the difference:
boxplot(weight ~ Outcome_2, data = data)

# categorical vs. categorical
# Use table objects to compare two categorical variables:
t1 = table(data$Treatment, data$Outcome_1)
t2 = table(data$Treatment, data$Outcome_2)
# We can inspect the tables:
t1
t2
# And make mosaic plots:
mosaicplot(t1, color = T)
mosaicplot(t2, color = T)


################################################
# plots
# The following website has some useful recipes for basic plots:
# http://www.harding.edu/fmccown/r/
# Remember that you can scroll through past plots with the arrow buttons in the 
# plot window.  You can also export a png of your plot, which you can then insert
# into a report or slide presentation.

# Example 1
# For this example, we'll redraw the scatterplot comparing height and weight.
# We're changing the dot style (pch), the plot title (main), 
# the axis labels (xlab and ylab), and the axis dimensions (xlim and ylim).
plot(data$height, data$weight,
     pch = 20,
     main = "Weight vs. Height",
     xlab = "height", ylab = "weight",
     xlim = c(50, 90), ylim = c(95, 280))
# Next we'll overlay a line of best fit obtained from 
# fitting a linear regression model:
fit = lm(weight ~ height, data = data)
abline(coef(fit)[1], coef(fit)[2], col = "red")
# And we'll also add a loess curve.  Thanks to Natalie for finding
# my typo in this code:
local_fit = loess.smooth(data$height, data$weight)
lines(local_fit[[1]], local_fit[[2]], 
      col = "blue", 
      lty = 2, 
      lwd = 2)
# Lastly, we're adding a legend.  As Andrea pointed out, the legend might
# hide some of your data points.  You can scale down the legend 
# with the "cex" option.
legend("bottomright", 
       c("Regression Line", "Loess Curve"), 
       lty = c(1, 2), 
       lwd = c(1, 2),
       col = c("red", "blue"),
       cex = 0.8)

# Example 2
# We also looked at a histogram example.  First we plot a histogram showing 
# densities rather than frequencies:
hist(data$weight, breaks = 20, freq = F)
# Next we can add a vertical line for the sample mean:
abline(v = mean(data$weight), col = "red")
# Lastly, we can overlay a normal density curve:
x = 80:300
y = dnorm(x, mean = mean(data$weight), 
          sd = sd(data$weight))
lines(x, y, col = "blue", lty = 2)


################################################
# functions

# We first looked at a function for counting missing values in a data frame:

count_NA = function(df){
  # find the number of columns (and therefore the number of variables)
  cols = ncol(df)
  # create an empty vector to store counts
  counts = rep(NA, cols)
  # iterate through the columns, counting NA's and storing the counts
  for (i in 1:cols){
    counts[i] = length(which(is.na(df[ , i])))
  }
  # add names to the vector entries for convenience
  names(counts) = names(data)
  # output the vector of counts
  return(counts)
}

# Now we can run the function:
count_NA(data)

# In class we talked about the difference between a complete case (a subject for whom
# we have data for all variables) and an incomplete case (a subject for whom we are
# missing one or more data points).  The class assignment was to write a function that
# identifies incomplete cases.  Rolondo's solution was as follows:

incomplete_cases = function(df){
  num = nrow(df)
  fi = rep(NA, num)
  for (i in 1:num) {
    fi[i] = length(which(is.na(df[i, ])))
  }
  end = data.frame("nrows" = 1:num, "NA" = fi)
  return(end)
}

# This is a very good solution.  It returns a data frame with a count of missing data
# for each subject:
incomplete_cases(data)

# I also showed this solution, which is similar to the version of the function that
# I use:  
incomplete_cases_2 = function(df){
  num = nrow(df)
  fi = rep(NA, num)
  for (i in 1:num) {
    if (length(which(is.na(df[i, ]))) > 0) {fi[i] = 1}    
  }
  return(which(fi == 1))
}

# This version returns a vector of row indices for incomplete cases:
incomplete_cases_2(data)
# We can use this to count incomplete cases:
length(incomplete_cases_2(data))
# And also to create a "complete case" data frame, which might be used 
# for a crude preliminary analysis:
ic = incomplete_cases_2(data)
cc_data = data[-ic, ]
# Note that in general one must be very careful about dropping incomplete cases,
# and there are other (much better) ways of dealing with missing data than
# simply conducting a complete case analysis.

