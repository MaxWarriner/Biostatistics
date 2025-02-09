################## Hypothesis Test Code for Biostatistics #################

########################################################################
# One-Sample T-Test
#######################################################################

xbar=mean(DATA)     # Mean of the hg data
s=sd(DATA)          # Standard deviation of the hg data
n=length(DATA)      # Sample size for the hg data
DoF = n-1           # Degrees of freedom for the hg data
alpha=0.05
sem = qt(1 - alpha / 2, DoF)*s/sqrt(n)   # SEM for the hg data


# Verify Assumptions

#Histogram
hist(DATA,freq = F, xlab = "Mercury Concentration (ppb)", main = "")
curve(dnorm(x, mean(DATA), sd(DATA)), min(DATA)-1,max(DATA)+1,col = "red",lwd = 2, add = TRUE)

#qqplot
qqnorm(DATA,main = "", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(DATA)

#Shapiro Test for Normality
shapiro.test(DATA)

#Outliers
boxplot(DATA)
outliers = boxplot.stats(DATA)$out
outliers_locations = which(DATA %in% outliers)
DATA[outliers_locations]
# 
# Perform T-test
t.test(DATA, mu = "FILL IN", alternative = "two.sided")
# 
# Power Analysis
xbar = mean(DATA)
s = sd(DATA)
mu = "FILL IN"
n = length(DATA)
cohen_d = (xbar - "FILL IN") / s
#
#
# Post-Hoc
power_analysis = pwr.t.test(n = n,d = cohen_d,sig.level = 0.05,type = "one.sample",alternative = "two.sided")
#
# Apriori
sample_size = pwr.t.test(power = 0.8, d = cohen_d,sig.level = 0.05, type = "one.sample", alternative = "two.sided")
#
# 
#
#############################################################################################
# Two Sample T-Test
############################################################################################


alpha = 0.05

# Assumptions

#Histograms
hist(DATA_A, main="", freq=FALSE, xlab=xlabel, ylab="Relative Frequency A")
curve(dnorm(x, mean=mean(DATA_A), sd=sd(DATA_A)), from=min(DATA_A)-1, to=max(DATA_A)+1, add=TRUE, col="red", lwd=2)

hist(DATA_B, main="", freq=FALSE, xlab=xlabel, ylab="Relative Frequency B")
curve(dnorm(x, mean=mean(DATA_B), sd=sd(DATA_B)), from=min(DATA_B)-1, to=max(DATA_B)+1, add=TRUE, col="red", lwd=2)

#qqplots
qqnorm(DATA_A, main="", xlab="Theoretical Quantiles A", ylab="Sample Quantiles B")
qqline(DATA_A, col="blue")

qqnorm(DATA_B, main="", xlab="Theoretical Quantiles A", ylab="Sample Quantiles B")
qqline(DATA_B, col="blue")

#Shapiro test for normality
shapiro.test(DATA_A)

shapiro.test(DATA_B)

#Boxplots for Outliers
boxplot(DATA_A)                                   # A boxplot to visually inspect for outliers, which are represented by circles.
outliers = boxplot.stats(DATA_A)$out              # Calculating the statistical outliers.
outliers_locations = which(DATA_A %in% outliers)  # Identifying the locations of the outliers.
DATA_A[outliers_locations]                        # Display the outliers in the data.

boxplot(DATA_B)                                   
outliers = boxplot.stats(DATA_B)$out              
outliers_locations = which(DATA_B %in% outliers)  
DATA_B[outliers_locations]                        

#Variance Test
var.test(DATA_A, DATA_B)


# Choose the right test
# Change alternative = "two.sided"  to "greater" or "less" as needed
t.test(DATA_A, DATA_B, var.equal = TRUE, alternative = "two.sided")  # Student's Two-sample T-test
t.test(DATA_A, DATA_B, var.equal = FALSE, alternative = "two.sided")  # Welch's Two-sample T-test


#Power for student's two-sample t-test

sP <- sqrt(((length(Data_A) - 1) * sd(Data_A)^2 + (length(Data_B) - 1) * sd(Data_B)^2) / (length(Data_A) + length(Data_B) - 2))
cohens_d = (mean(Data_A) - mean(Data_B)) / sP    # Calculating Cohen's d, a measure of effect size
library(pwr) # Load necessary library

n = length(Data_A) + length(Data_B)                        # Calculating total sample size
# If the total sample size is less than 20, we can correct Cohen's d using Hedges's g formula
# if (n < 20) {cohens_d = cohens_d * ((n - 3) / (n - 2.25)) * sqrt((n - 2) / n)}

# Post-hoc power
pwr.t2n.test(n1 = length(Data_A), n2 = Length(Data_B), d = cohens_d, sig.level = 0.05, power = NULL, alternative = "two.sided")


# #Apriori Power

desired_power = 0.8 # Replace with your desired power

# Define the sample size for the first group. You can play with this value to attain equal sample sizes.
nsA = Length(Data_A)

# Calculating required sample size n2 for given power and n1
pwr.t2n.test(n1 = nsA, n2 = NULL, d = cohens_d, sig.level = 0.05, power = desired_power, alternative = "two.sided")


#Post-hoc for Welch's
d = (xbarA - xbarB)/sqrt(sA^2/2 + sB^2/2) # Effect size: cohens_d
source("pwr.welch.test.R")
pwr.welch.test(nA = Length(Data_A), nB = length(Data_B), xbarA = mean(Data_A), xbarB = mean(Data_B), sA = sd(Data_A), sB = sd(Data_B), n.simulations = 1000, alternative = "two.sided")

# Apriori
# # Using the custom pwr.welch.test function to find required sample size n2 for given power and n1
pwr.welch.test(nA = Length(Data_A), nB = NULL, xbarA = mean(Data_A), xbarB = mean(Data_B), sA = sd(Data_A), sB = sd(Data_B), power = 0.8, sig.level = 0.05, n.simulations = 1000, alternative = "two.sided")
# 
# 
###########################################################################################
################### Sign Test(Non Parametric One-Sample T-Test) ###########################
###########################################################################################
#
library(BSDA)
SIGN.test(DATA,md="FILL IN",alternative = "two.sided")
# 
# 
##########################################################################################
################## Wilcoxon Signed Rank Test(Non Parametric One-Sample T-Test)
#########################################################################################


# Assumptions

#Histogram
hist(DATA, main="", freq=F,)

#Skewness
library(e1071) # Loading the necessary library
skewness(DATA) # If the skewness is between -1 and 1, we can assume the data is symmetric

# Perform the test
wilcox.test(DATA, mu = "FILL IN", alternative = "two.sided")

# Power analysis: See Mann Whitney U Test
##########################################################################################
################## Mann-Whitney U Test (Non Parametric Two-Sample T-Test)#############
###########################################################################################

# Assumptions

#Histograms
hist(DATA_A, main="", freq=F, ylab="Relative Frequency", col=rgb(1/2,1/2,1/2,1/2))
curve(dnorm(x, mean(DATA_A), sd(DATA_A)), min(DATA_A)-1, max(DATA_A)+1, add=TRUE, col="red", lwd=2)

hist(DATA_B, main="", freq=F, ylab="Relative Frequency", col=rgb(1/2,1/2,1/2,1/2))
curve(dnorm(x, mean(DATA_B), sd(DATA_B)), min(DATA_B)-1, max(DATA_B)+1, add=TRUE, col="blue", lwd=2)

#Similarity of Distribution
ks.test(DATA_A, DATA_B) #If similar distributions, comparing medians. If not, comparing mean ranks

# Perform Test
wilcox.test(DATA_A, DATA_B, alternative = "two.sided")

# Power Analysis

# Fit the distributions

library(fitdistrplus) # Load the required package
x <- DATA_A
distributions <- c("norm", "weibull", "gamma", "logis", "exp", "pois", "geom", "nbinom")
fit_list <- lapply(distributions, function(dist) {
  tryCatch(fitdist(x, dist), error = function(e) NULL)
})
names(fit_list) <- distributions
fit_list <- Filter(Negate(is.null), fit_list)
# Extract fit statistics and parameters into a data frame
fit_stats <- do.call(rbind, lapply(names(fit_list), function(dist) {
  fit <- fit_list[[dist]]
  params <- paste(names(coef(fit)), "=", round(coef(fit), 4), collapse = ", ")
  data.frame(Distribution = dist, AIC = fit$aic, BIC = fit$bic, LogLikelihood = fit$loglik, Parameters = params)
}))
print(fit_stats) # Display fit statistics

#Run the Power Analysis

source("pwr.wilcox.test.R")
nA=length(Data_A)
nB=length(Data_B)
distA = "FILL IN"
distB = "FILL IN"
paramsA = list("FILL IN")
paramsB = list("FILL IN")
pwr.wilcox.test(nA, nB, dist1 = distA, dist2 = distB, params1 = paramsA, params2 = paramsB)

############################################################################################
########################## One Way ANOVA - Fixed Effect #####################################
############################################################################################

library(car)
library(pwr)
library(emmeans)
library(rstatix)

#Fit the model
my.anova = lm("DEPENDENT VARIABLE" ~ "INDEPENDENT VARIABLE", data = DATA)

# Check Assumptions

#Histogram of the residuals
hist(my.anova$residuals, freq=FALSE, xlab="Residuals", main="Residuals Histogram")
curve(dnorm(x, mean(my.anova$residuals), sd(my.anova$residuals)), add=TRUE, col='red')

#qqplot of the residuals
qqnorm(my.anova$residuals); qqline(my.anova$residuals)

#Shapiro test for normality of residuals
shapiro.test(my.anova$residuals)

# Boxplot for outliers
boxplot("DEPENDENT VARIABLE" ~ "INDEPENDENT VARIABLE", data=DATA, ylab="FILL IN", xlab="FILL IN")
get_outliers = function(data) {which(data %in% identify_outliers(data.frame(data))[[1]])}
EX: GROUP1.OUTLIERS = get_outliers(Group1)

# Variance Test
bartlett.test("DEPENDENT VARIABLE" ~ "INDEPENDENT VARIABLE", data = DATA) # Bartlett's Test

# Follow-Up Tukey Test
emmeans_results = emmeans(my.anova, ~ "INDEPENDENT VARIABLE")
pairs(emmeans_results, adjust = "tukey")

#Power Analysis
# 
# Post-Hoc
# # Extract SSB and SSW for effect size calculation
SSB = anova_df["INDEPENDENT VARIABLE", "Sum Sq"]
SSW = anova_df["Residuals", "Sum Sq"]
eta_squared = SSB / (SSB + SSW)
fvalue = sqrt(eta_squared / (1 - eta_squared))
# Calculate harmonic mean of group sizes for the power analysis
group_sizes = table(DATA$INDEPENDENT.VARIABLE)
harmonic_mean_n = length(group_sizes) / sum(1 / group_sizes)

# Perform post hoc power analysis
pwr.anova.test(k = length(group_sizes), n = harmonic_mean_n, f = fvalue, sig.level = 0.05)

# Apriori
desired_power = 0.8
estimated_effect_size = fvalue
# Calculate required sample size for each group
pwr.anova.test(k = length(group_sizes), f = estimated_effect_size, sig.level = 0.05, power = 0.8)
# 
##############################################################################################
############################### One Way ANOVA - Random Effect ################################
###########################################################################################3
# 
# Same assumptions as Fixed Effect ANOVA

#Histogram of the residuals
hist(my.anova$residuals, freq=FALSE, xlab="Residuals", main="Residuals Histogram")
curve(dnorm(x, mean(my.anova$residuals), sd(my.anova$residuals)), add=TRUE, col='red')

#qqplot of the residuals
qqnorm(my.anova$residuals); qqline(my.anova$residuals)

#Shapiro test for normality of residuals
shapiro.test(my.anova$residuals)

# Boxplot for outliers
boxplot("DEPENDENT VARIABLE" ~ "INDEPENDENT VARIABLE", data=DATA, ylab="FILL IN", xlab="FILL IN")
get_outliers = function(data) {which(data %in% identify_outliers(data.frame(data))[[1]])}
EX: GROUP1.OUTLIERS = get_outliers(Group1)

# Variance Test
bartlett.test("DEPENDENT VARIABLE" ~ "INDEPENDENT VARIABLE", data = DATA) # Bartlett's Test


#Fit the model and run the test
my.anova = lmer("DEPENDENT VARIABLE" ~ (1|"INDEPENDENT VARIABLE"), data = DATA)
ranova(my.anova)


#Follow Up Variance Test
random_effect = "INDEPENDENT VARIABLE"
var_components = VarCorr(my.anova)
random_variance = as.numeric(attr(var_components[[random_effect]], "stddev"))^2
residual_variance = attr(var_components, "sc")^2
total_variance = random_variance + residual_variance
percent_random_variance = (random_variance / total_variance) * 100
percent_residual_variance = (residual_variance / total_variance) * 100
cat("Variance explained by", random_effect, "(random factor):", round(percent_random_variance, 2), "%\n")
cat("Residual variance (unexplained):", round(percent_residual_variance, 2), "%\n")

# Post-Hoc Power
source("pwr.model2_oneway_anova.R")
calculated_power = pwr.model2_oneway_anova(data = DATA, response_var = "DEPENDENT VARIABLE", grouping_var = "DEPENDENT VARIABLE",  n_sim = 100,max_singular_rate = 0.05)
cat("Post hoc power: ", calculated_power, "\n")

# Aprioro
summary = DATA %>%group_by("INDEPENDENT VARIABLE") %>% summarise(Mean = mean("DEPENDENDENT VARIABLE", na.rm = TRUE),SD = sd("DEPENDENT VARIABLE", na.rm = TRUE))
required_nA = pwr.model2_oneway_anova(n_groups = length(unique(dung$Species)), group_means = summary$Mean, within_sd = summary$SD, desired_power = 0.8, n_sim = 100, start_n=3)
cat("Required sample size per group to achieve power = 0.8: ", required_nA, "\n")

##########################################################################################
##################### Kruskwal Wallis (Non-Parametric One Way ANOVA)#####################
##########################################################################################

# Check Assumptions

#Similarity of distributions
boxplot("DEPENDENT VARIABLE" ~ "INDEPENDENT VARIABLE", data=DATA, col=c("lightgray", "white", "darkgray"), cex.lab=1.5, cex.axis=1.5, cex=1.5)

#Pairwise ks tests
p1 = ks.test(DATA_A, DATA_B)$p.value
p2 = ks.test(DATA_A, DATA_C)$p.value
p3 = ks.test(DATA_B, DATA_C)$p.value
p.adjust(c(p1, p2, p3), method="BH")

# Perform Test
kruskal.test("DEPENDENT VARIABLE" ~ "INDEPENDENT VARIABLE", data=DATA)


# Follow-up: Dunn Test
dunnTest(DATA$"DEPENDENT.VARIABLE", DATA$"INDEPENDENT VARIABLE", method="bh")
# 
# Power Analysis
#
# Post-Hoc:
# 
#Fit distributions
library(fitdistrplus) # Load the library
x = DATA_A  # Replace with each data set
distributions = c("norm", "weibull", "gamma", "logis", "exp", "pois", "geom", "nbinom") # Distributions to fit
fit_list = list() # Fit distributions
for (dist in distributions) {
  fit = tryCatch({
    fitdist(x, dist)
  }, error = function(e) {
    cat("\nError in fitting", dist, ":", e$message, "\n")
    return(NULL)
  })
  if (!is.null(fit)) fit_list[[dist]] = fit
}

fit_stats = data.frame(Distribution = character(), AIC = numeric(), BIC = numeric(), LogLikelihood = numeric(), Parameters = character(), stringsAsFactors = FALSE)
for (dist in names(fit_list)) {
  fit = fit_list[[dist]]
  params_str = paste(names(coef(fit)), "=", round(coef(fit), 4), collapse = ", ")
  fit_stats = rbind(fit_stats, data.frame(Distribution = dist, AIC = fit$aic, BIC = fit$bic, LogLikelihood = fit$loglik, Parameters = params_str))
}
print(fit_stats) # Display fit stats

# Replace with all the fit distributions

#Post Hoc
source("pwr_kruskal_test.R")
pwr_kruskal_test(group_sizes = c(length(DATA_A), length(DATA_B), length(DATA_C)),
                 dist_list = list("DISTRIBUTION A", "DISTRIBUTION B", "DISTRIBUTION C"),
                 params_list = list(list("PARAMETERS FOR A"),
                                    list("PARAMETERS FOR B"),
                                    list("PARAMETERS FOR C")))

#Apriori
pwr_kruskal_test(desired_power = 0.8,
                 dist_list = list("DISTRIBUTION A", "DISTRIBUTION B", "DISTRIBUTION C"),
                 params_list = list(list("PARAMETERS FOR A"),
                                    list("PARAMETERS FOR B"),
                                    list("PARAMETERS FOR C")),
                 start_size = "FILL IN")


##########################################################################################
########################### Two Way ANOVA - Fixed Effects ################################
##########################################################################################

my.anova = lm(DependentVariable ~ factor1 * factor2, data = DATA)

#Histogram of the residuals
hist(my.anova$residuals, freq=FALSE, xlab="Residuals", main="Residuals Histogram")
curve(dnorm(x, mean(my.anova$residuals), sd(my.anova$residuals)), add=TRUE, col='red')

#qqplot of the residuals
qqnorm(my.anova$residuals); qqline(my.anova$residuals)

#Shapiro test for normality of residuals
shapiro.test(my.anova$residuals)

#Variance Test
bartlett.test("DEPENDENT VARIABLE" ~ interaction("INDEPENDENT VARIABLE 1", "INDEPENDENT VARIABLE 2"), data = DATA)

#Boxplot for outliers
library(rstatix)
boxplot("Dependent Variable" ~ interaction(Factor1, Factor2), data = DATA, ylab = "ylab", xlab = "xlab", cex.lab = 0.75, cex.axis = 0.75, cex = 1)
get_outliers = function(data) {
outliers = identify_outliers(data.frame(data))[[1]]
return(outliers)
}
interaction_outliers = split(DATA$dependent.variable, interaction(DATA$factor1, DATA$factor2))
outliers_list = lapply(interaction_outliers, get_outliers)
outliers_list # Display the outliers for each interaction pair


# Perform Test
library(car)
options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly"))
my.anova = lm(DependentVariable ~ factor1 * factor2, data = DATA)
Anova(my.anova, type = "III")

# Run Interaction Plot if Necessary
interaction.plot(x.factor = DATA$factor1, trace.factor = DATA$factor2, response = DATA$dependent.variable, type = "b", ylab = "ylab", xlab = "xlab", lty = c(1, 2, 3), lwd = 2, pch = c(0, 19, 2), legend = FALSE)
legend(x = "topright", legend = c("factorgroup", "factorgroup"), bty = "n", lty = c(1, 2, 3), lwd = 2, pch = c(0, 19, 2), title = "", inset = 0.02)

# Follow Up Tukey-Test
emm_results = emmeans(my.anova, ~ factor1 * factor2)
summary(emm_results)
pairs(emm_results, adjust = "tukey")

# Power Analysis
#library(pwr)
# Fit the two-way ANOVA model using lm
lm_model = lm(DependentVariable ~ factor1 * factor2, data = DATA)
anova_results = Anova(lm_model, type = "III")
# 
# Parameters for power and sample size analysis
alpha = 0.05            # Significance level
power_target = 0.8      # Desired power

# Extracting sums of squares
SS_A = anova_results["factor1", "Sum Sq"]
SS_B = anova_results["factor2", "Sum Sq"]
SS_AB = anova_results["factor1:factor2", "Sum Sq"]
SS_Error = anova_results["Residuals", "Sum Sq"]

# Calculate total sum of squares
SS_Total = SS_A + SS_B + SS_AB + SS_Error

# Calculate overall model effect size
f2_overall = (SS_A + SS_B + SS_AB) / SS_Error

# Calculate effect sizes for main effects and interaction
f2_A = SS_A / SS_Error
f2_B = SS_B / SS_Error
f2_AB = SS_AB / SS_Error

# Output effect sizes
cat("Effect size (Cohen's f²) for Overall Model:", f2_overall, "\n")
cat("Effect size (Cohen's f²) for Factor 1:", f2_A, "\n")
cat("Effect size (Cohen's f²) for Factor 2:", f2_B, "\n")
cat("Effect size (Cohen's f²) for Interaction:", f2_AB, "\n")

################################################################################
# Mixed Model Two Way ANOVA
################################################################################

library(lme4) # Load required libraries
library(car) 
library(dplyr)
library(rcompanion)
library(FSA)
library(emmeans)

# load as DATA and set independent variables as factors

# Set contrast options for Type III SS (necessary for correct Type III ANOVA results)
options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly"))
my.anova = lmer("DEPENDENT VARIABLE" ~ "FIXED FACTOR" + (1 | "RANDOM FACTOR") + (1 | "FIXED FACTOR:RANDOM FACTOR"), data = DATA) # Fit the model

# Assumptions are the same as Fixed Effects Two Way ANOVA

#Histogram of the residuals
hist(my.anova$residuals, freq=FALSE, xlab="Residuals", main="Residuals Histogram")
curve(dnorm(x, mean(my.anova$residuals), sd(my.anova$residuals)), add=TRUE, col='red')

#qqplot of the residuals
qqnorm(my.anova$residuals); qqline(my.anova$residuals)

#Shapiro test for normality of residuals
shapiro.test(my.anova$residuals)

#Variance Test
bartlett.test("DEPENDENT VARIABLE" ~ interaction("INDEPENDENT VARIABLE 1", "INDEPENDENT VARIABLE 2"), data = DATA)

#Boxplot for outliers
library(rstatix)
boxplot("Dependent Variable" ~ interaction(Factor1, Factor2), data = DATA, ylab = "ylab", xlab = "xlab", cex.lab = 0.75, cex.axis = 0.75, cex = 1)
get_outliers = function(data) {
  outliers = identify_outliers(data.frame(data))[[1]]
  return(outliers)
}
interaction_outliers = split(DATA$dependent.variable, interaction(DATA$factor1, DATA$factor2))
outliers_list = lapply(interaction_outliers, get_outliers)
outliers_list # Display the outliers for each interaction pair

#Run the Test
AnovaTable = Anova(my.anova, type = "III") # Perform Type III ANOVA
AnovaTable

# Fit the reduced model without the random effect
my.anova.no.random = lmer("DEPENDENT VARIABLE" ~ "FIXED FACTOR" + (1 | "FIXED FACTOR:RANDOM FACTOR"), data = DATA, REML = FALSE)
# Perform likelihood ratio test for the random effect
anova(my.anova.no.random, my.anova)

#Test significance of the interaction term
# Fit the reduced model without the interaction term 
my.anova.no.interaction = lmer("DEPENDENT VARIABLE" ~ "FIXED FACTOR" + (1 | "RANDOM FACTOR"), data = DATA, REML = FALSE)
# Perform likelihood ratio test for the interaction term
anova(my.anova.no.interaction, my.anova)

#Follow up test for fixed factor
emmeans(my.anova, pairwise ~ "Fixed FACTOR", adjust = "tukey")


#Follow up for random factor
anova_results = Anova(my.anova, type = "III") # Perform Type III ANOVA
variance_components = VarCorr(my.anova) # Extract variance components
residual_variance = attr(variance_components, "sc")^2 # Extract residual variance

# Fixed effect variance for Fixed effect (derived from model residuals)
fixed_effect_variance = sum((predict(my.anova) - mean(predict(my.anova)))^2) / length(predict(my.anova))
# Extract variance for random effects
random_variance = attr(variance_components$"RANDOM FACTOR", "stddev")^2
interaction_variance = attr(variance_components$"FIXED FACTOR:RANDOM FACTOR", "stddev")^2
# Calculate total variance (fixed effect, random effects, and residual variance)
total_variance = fixed_effect_variance + random_variance + interaction_variance + residual_variance

# Calculate percent variance contribution only if total variance is not zero
percent_variance = c("FIXED FACTOR" = (fixed_effect_variance / total_variance) * 100,
                     "RANDOM FACTOR" = (random_variance / total_variance) * 100,
                     "FIXED FACTOR:RANDOM FACTOR" = (interaction_variance / total_variance) * 100,
                     "Residuals" = (residual_variance / total_variance) * 100)
percent_variance

# Power analysis

#Post hoc
source("pwr_mixedmodel_twoway_anova.R")
nA =   # Number of levels for factor1
nB =  # Number of levels for factor2

statistical_summary = DATA %>%
  group_by("RANDOM FACTOR", "FIXED FACTOR") %>%
  summarize(means = mean("DEPENDENT VARIABLE", na.rm = TRUE), sds = sd("DEPENDENT VARIABLE", na.rm = TRUE),ncells=n(), .groups = "drop")
means = matrix(statistical_summary$means, nrow=nA)
sds = matrix(statistical_summary$sds, nrow=nA)
n_cells = matrix(statistical_summary$ncells, nrow=nA)

pwr_mixedmodel_twoway_anova(means = means, sds = sds, nA = nA, nB = nB, desired_power = NULL, n_cells = n_cells, min_n = 5, max_n = 100, n_sim = 100, alpha = 0.05, verbose = FALSE) 

#Apriori
desired_power = 0.80
pwr_mixedmodel_twoway_anova(means = means, sds = sds, nA = nA, nB = nB, desired_power = desired_power, n_cells = NULL, min_n = 20, max_n = 21,n_sim = 100, alpha = 0.05, verbose = FALSE) 


################################################################################
# Random-Random Two Way ANOVA
################################################################################
library(lme4) # Load required libraries
library(car) 
library(dplyr)
library(rcompanion)
library(FSA)
library(emmeans)

# Fit the mixed model treating both as random effects
my.anova = lmer("DEPENDENT VARIABLE" ~ (1 | "RANDOM FACTOR 1") + (1 | "RANDOM FACTOR 2") + (1 | "RANDOM FACTOR 1:RANDOM FACTOR 2"), data = DATA)

# Same assumptions as other Two Way ANOVA

#Histogram of the residuals
hist(my.anova$residuals, freq=FALSE, xlab="Residuals", main="Residuals Histogram")
curve(dnorm(x, mean(my.anova$residuals), sd(my.anova$residuals)), add=TRUE, col='red')

#qqplot of the residuals
qqnorm(my.anova$residuals); qqline(my.anova$residuals)

#Shapiro test for normality of residuals
shapiro.test(my.anova$residuals)

#Variance Test
bartlett.test("DEPENDENT VARIABLE" ~ interaction("INDEPENDENT VARIABLE 1", "INDEPENDENT VARIABLE 2"), data = DATA)

#Boxplot for outliers
library(rstatix)
boxplot("Dependent Variable" ~ interaction(Factor1, Factor2), data = DATA, ylab = "ylab", xlab = "xlab", cex.lab = 0.75, cex.axis = 0.75, cex = 1)
get_outliers = function(data) {
  outliers = identify_outliers(data.frame(data))[[1]]
  return(outliers)
}
interaction_outliers = split(DATA$dependent.variable, interaction(DATA$factor1, DATA$factor2))
outliers_list = lapply(interaction_outliers, get_outliers)
outliers_list # Display the outliers for each interaction pair

#Test factor 1
# Fit the reduced model without the random effect 1
model_no_factor_1 = lmer("DEPENDENT VARIABLE" ~ (1 | "FACTOR 2") + (1 | "FACTOR 1:FACTOR 2"), data = DATA, REML = FALSE)
# Perform likelihood ratio test for the random effect 1
anova(model_no_factor_1, my.anova)

#Test factor 2
# Fit the reduced model without the random effect 2
model_no_factor_2 = lmer("DEPENDENT VARIABLE" ~ (1 | "FACTOR 1") + (1 | "FACTOR 1:FACTOR 2"), data = DATA, REML = FALSE)
# Perform likelihood ratio test for the random effect 2
anova(model_no_factor_2, my.anova)

#Test interaction
# Fit the reduced model without the interaction term
model_no_interaction = lmer("DEPENDENT VARIABLE" ~ (1 | "FACTOR 1") + (1 | "FACTOR 2"), data = DATA, REML = FALSE)
# Perform likelihood ratio test for the interaction effect
anova(model_no_interaction, my.anova)


#Calculate the variance for each effects
anova_results = Anova(my.anova, type = "III") # Perform Type III ANOVA
variance_components = VarCorr(my.anova) # Extract variance components
residual_variance = attr(variance_components, "sc")^2 # Extract residual variance

# Calculate the variance for random effects
factor_1_variance = attr(variance_components$Factor1,"stddev")^2
factor_2_variance = attr(variance_components$Factor2, "stddev")^2
interaction_variance = attr(variance_components$Factor1:Factor2, "stddev")^2
# Calculate total variance (random effects and residual variance)
total_variance = factor_1_variance + factor_2_variance + interaction_variance + residual_variance

# Calculate percent variance contribution only if total variance is not zero
percent_variance = c("Factor 1" = (factor_1_variance / total_variance) * 100,
                     "Factor 2" = (factor_2_variance / total_variance) * 100,
                     "Interaction" = (interaction_variance / total_variance) * 100,
                     "Residuals" = (residual_variance / total_variance) * 100)
percent_variance

#Power analysis
source("pwr_model2_twoway_anova.R")
nA =  # Number of levels for factor1
nB =  # Number of levels for factor2

statistical_summary = DATA %>%
  group_by(factor1, factor2) %>%
  summarize(means = mean(dependentvariable, na.rm = TRUE), sds = sd(dependentvariable, na.rm = TRUE), ncells=n(), .groups = "drop")
means = matrix(statistical_summary$means, nrow=nA)
sds = matrix(statistical_summary$sds, nrow=nA)
n_cells = matrix(statistical_summary$ncells, nrow=nA)

pwr_model2_twoway_anova(means = means, sds = sds, nA = nA, nB = nB, desired_power = NULL, n_cells = n_cells, min_n = 5, max_n = 100, n_sim = 100, alpha = 0.05, verbose = FALSE) 

################################################################################
#Scheirer-Ray-Hare Test (Non-Parametric Two Way ANOVA)
################################################################################

#load dataset and set factors as as.factor()

#Perform test
testresult=scheirerRayHare(DependentVariable ~ Factor1 * Factor2, data = DATA)
testresult

#Test for for medians vs. mean ranks
# Create a boxplot
boxplot(DependentVariable ~ interaction(Factor1, Factor2), data = DATA, 
        ylab = "ylab", xlab = "xlab", cex.lab = 1.2, cex.axis = 1.2)

# Filter the dataset
G1 = dplyr::filter(DATA, factor1 == "level 1" & factor2 == "level 1")$Dependent.Variable
G2 = dplyr::filter(DATA, factor1 == "level 1" & factor2 == "level 2")$Dependent.Variable
#repeat for each combination

# Perform Kolmogorov-Smirnov tests to compare distributions between groups
ks1 = ks.test(G1, G2)$p.value
#repeat for each comparison

pvalues = c(ks1, ks2, ks3, ks4, ks5, ks6) # Combine p-values from all tests into a vector
p.adjust(pvalues, method = 'BH')           # Adjust the p-values for multiple testing



#create interaction plot is interaction term is significant
interaction.plot(
  x.factor = DATA$factor1,
  trace.factor = DATA$factor2,
  response = DATA$DependentVariable,
  fun = median,
  type = "b",
  ylab = "ylab",
  xlab = "xlab",
  lty = c(1, 2, 3),
  lwd = 2,
  pch = c(0, 19, 2),
  legend = FALSE
)

# Add legend to the bottom left of the plot
legend(
  "bottomleft",
  legend = c("level1", "level2"),
  bty = "n",
  lty = c(1, 2, 3),
  lwd = 2,
  pch = c(0, 19, 2),
  title = "",
  inset = .02
)


#Follow up Dunn Test
dunnTest(DATA$DependentVariable, DATA$factor1, method="bh") # Gives error if you use it with 2 groups
dunnTest(DATA$DependentVariable, DATA$factor2, method="bh") # Gives error if you use it with 2 groups
dunnTest(DATA$DependentVariable, interaction(DATA$factor1, DATA$factor2), method="bh")


#Power Analysis
#Fit distributions
library(fitdistrplus) # Load the library
x = G1  # Replace with actual dataset
distributions = c("norm", "weibull", "gamma", "logis", "exp", "pois", "geom", "nbinom") # Distributions to fit

fit_list = list() # Fit distributions
for (dist in distributions) {
  fit = tryCatch({
    fitdist(x, dist)
  }, error = function(e) {
    cat("\nError in fitting", dist, ":", e$message, "\n")
    return(NULL)
  })
  if (!is.null(fit)) fit_list[[dist]] = fit
}

# Gather goodness-of-fit stats
fit_stats = data.frame(Distribution = character(), AIC = numeric(), BIC = numeric(), LogLikelihood = numeric(), Parameters = character(), stringsAsFactors = FALSE)
for (dist in names(fit_list)) {
  fit = fit_list[[dist]]
  params_str = paste(names(coef(fit)), "=", round(coef(fit), 4), collapse = ", ")
  fit_stats = rbind(fit_stats, data.frame(Distribution = dist, AIC = fit$aic, BIC = fit$bic, LogLikelihood = fit$loglik, Parameters = params_str))
}

print(fit_stats) # Display fit stats

#Post-Hoc - replace as needed

source("pwr_srh_test.R")
# Example usage for 2x2 factorial design:
gs=3 # Group size

pwr_srh_test(group_sizes = c(gs, gs, gs, gs), desired_power = NULL,
             dist_list = list("norm", "norm", "norm", "norm"),
             params_list = list(
               list(mean = 695.6667, sd = 12.4722),
               list(mean = 642.6667, sd = 35.3679),
               list(mean = 535.3333, sd = 47.3943),
               list(mean = 517.3333, sd = 15.3695)
             ), n.simulations = 100, alpha = 0.05, n_factors = 2)


#Apriori - replace as needed
pwr_srh_test(group_sizes = NULL, desired_power = 0.99,
             dist_list = list("norm", "norm", "norm", "norm"),
             params_list = list(
               list(mean = 695.6667, sd = 12.4722),
               list(mean = 642.6667, sd = 35.3679),
               list(mean = 535.3333, sd = 47.3943),
               list(mean = 517.3333, sd = 15.3695)
             ), n.simulations = 100, alpha = 0.05, start_size = 1, n_factors = 2)


################################################################################
#Randomized Block ANOVA
################################################################################

#Assumptions

# Fit the Randomized Block Design ANOVA model
my.rbanova = lm(DependentVariable ~ Block + Factor, data = DATA)

#Histogram for normality
hist(my.rbanova$residuals, freq=FALSE, main="", xlab="Residuals", cex.lab=1.2,cex.axis=1.2) #Histogram
curve(dnorm(x,mean(my.rbanova$residuals), sd(my.rbanova$residuals)), add=TRUE,col="red")

#qqplot
qqnorm(my.rbanova$residuals, main="",cex.lab=1.2,cex.axis=1.2) # QQ Plot
qqline(my.rbanova$residuals)

#Shapiro test
shapiro.test(my.rbanova$residuals)

#Boxplot for outliers
boxplot(DependentVariable ~ Factor, data = DATA, names = c("Group1","Group2"), ylab = "ylab", cex.lab=1.5, cex.axis=1.5, cex=1)
# Identify and display outliers for each treatment group
outliers_df = DATA %>% group_by(Factor) %>%
  summarise(Outliers = list(boxplot.stats(DependentVariable)$out)) %>% unnest(Outliers)
outliers_df # Display the data frame of outliers by treatment group

#Levene Test for Variance
leveneTest(DependentVariable ~ factor, data = DATA) # Note: The blocking factor is not included.

#Run the ANOVA test
AnovaTable = anova(my.rbanova)
AnovaTable

#Test for Relative Efficiency
# Extract Mean Squares for Blocks (MS_B) and Mean Squares for Error (MS_E)
MS_B = AnovaTable["Block", "Mean Sq"]  # Mean Square for Block
MS_E = AnovaTable["Residuals", "Mean Sq"]  # Mean Square for Error
b =   # Number of blocks
a =   # Number of treatments
relative_efficiency = ((b - 1) * MS_B + b * (a - 1) * MS_E) / ((a * b - 1) * MS_E) # Calculate Relative Efficiency
print(paste("Relative Efficiency:", relative_efficiency)) # Print the result

#Follow up Tukey Test
# Calculate the estimated marginal means for the variable using the ANOVA model (my.rbanova)
marginal_means = emmeans(my.rbanova, ~ Factor)
marginal_means
# Perform pairwise comparisons between the levels
pairs(marginal_means, adjust="tukey")

#Power Analysis
# Fit the two-way ANOVA model using lm
lm_model = lm(DependentVariable ~ Block + Factor, data = DATA)
anova_results = Anova(lm_model, type = "III")  
# Parameters for power and sample size analysis
alpha = 0.05            # Significance level
power_target = 0.8      # Desired power
# Extracting sums of squares
SS_A = anova_results["Block", "Sum Sq"]
SS_B = anova_results["Factor", "Sum Sq"]
SS_Error = anova_results["Residuals", "Sum Sq"]
# Calculate overall model effect size
f2_overall = (SS_A + SS_B) / SS_Error
# Calculate effect sizes for main effects and interaction
f2_B = SS_B / SS_Error
# Output effect sizes
cat("Effect size (Cohen's f²) for Overall Model:", f2_overall, "\n")
cat("Effect size (Cohen's f²) for Factor 2:", f2_B, "\n")

#Post-Hoc Power Analysis
# Calculate power based on initial sample size for each factor, interaction term, and overall model
# Calculate u for the entire model
my.u =  sum(anova_results$Df[2:3]) 
my.v = anova_results$Df[4]  
power_calculated_overall = pwr.f2.test(u = my.u, v = my.v, f2 = f2_overall, sig.level = alpha)
my.u2 =  sum(anova_results$Df[3]) 
power_calculated_main2 = pwr.f2.test(u = my.u2, v = my.v, f2 = f2_B, sig.level = alpha)
# Display calculated power
cat("Calculated power for the sample size for Overall Model is:", power_calculated_overall$power, "\n")
cat("Calculated power for the sample size for main effect of Factor 2 is:", power_calculated_main2$power, "\n")

#Apriori
# Sample size calculation for overall model
power_target = 0.8
my.u =  sum(anova_results$Df[2:3]) 
my.v = anova_results$Df[4] 
sample_size_result_overall = pwr.f2.test(u = my.u, v = NULL, f2 = f2_overall, sig.level = alpha, power = power_target)        
# Sample size calculation for main effect of Factor 2
my.u2 =  sum(anova_results$Df[3]) 
sample_size_result_main2 = pwr.f2.test(u = my.u2, v = NULL, f2 = f2_B, sig.level = alpha, power = power_target)    
# Calculate required total sample size (N) for each calculation
required_sample_size_overall = ceiling(sample_size_result_overall$v) + my.u +1
required_sample_size_main2 = ceiling(sample_size_result_main2$v) + my.u +1
# Output results for main factors, interaction term, and overall model
cat("Required sample size per group to achieve power of 0.8 for Overall Model:", required_sample_size_overall, "\n")
cat("Required sample size per group to achieve power of 0.8 for main effect of Factor 2:", required_sample_size_main2, "\n")


########################################################################################################
# One-Way Repeated Measures ANOVA
#######################################################################################################

#Assumptions

# Check for Normality
# Loop through each and perform Shapiro-Wilk test
# Apply Shapiro-Wilk normality test within each group and adjust p-values
normality_results = DATA %>%
  group_by(Factor) %>%
  summarise(
    W = shapiro.test(DependentVariable)$statistic,  # Extract W value
    p_value = shapiro.test(DependentVariable)$p.value,  # Extract p-value
    .groups = 'drop'  # Drop grouping after summarise to avoid warning
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))  # Adjust p-values
normality_results

# Plot QQ plots for each level to visually assess normality
DATA %>% ggplot(aes(sample = DependentVariable)) + facet_wrap(~ Factor) + stat_qq() + stat_qq_line()
DATA %>% ggplot(aes(x = DependentVariable)) + facet_wrap(~ Factor) + geom_histogram(aes(y = ..density..), binwidth , color = "black", fill = "lightblue") + stat_function(fun = dnorm, args = list(mean = mean(DATA$DependentVariable, na.rm = TRUE), sd = sd(DATA$DependentVariable, na.rm = TRUE)), color = "red") + labs(x = "DependentVariable", y = "Density") + theme_minimal()

#Boxplot for Outliers
# Plot a boxplot for across each group
boxplot(DependentVariable ~ Factor, data = DATA,  ylab = "DependentVariable", cex.lab=1.5, cex.axis=1.5, cex=1)
# Identify and display outliers for each treatment group
outliers_df = DATA %>% group_by(Factor) %>% summarise(Outliers = list(boxplot.stats(DependentVariable)$out)) %>% unnest(Outliers)
outliers_df # Display the data frame of outliers by treatment group

#Run the ANOVA test
# Perform one-way repeated measures ANOVA with afex, ignoring Age
rm_oneway_result <- aov_ez(
  id = "Person",
  dv = "DependentVariable",
  within = "Factor",            
  data = DATA,
  type = 3,                    # Use type 3 sums of squares
  anova_table = list(correction = "none", es = "ges")
)
# Display the ANOVA results
summary(rm_oneway_result)

# If the sphericity condition had been violated, a correction would be necessary. 
# The choice of correction depends on the Greenhouse-Geisser epsilon (GGe).If GGe 
# is greater than 0.75 or if the sample size is small (e.g., 10), the Huynh-Feldt 
# epsilon should be used.

#Pairwise Post-Hoc Tests
# Conduct post-hoc pairwise comparisons for Factor levels
pairwise_results = emmeans(rm_oneway_result, pairwise ~ Factor, adjust = "BH") # Bonferroni adjustment
pairwise_results$contrasts

#Power Analysis
# Set the seed for reproducibility
# set.seed(0)
source("pwr_rm_anova.R")
gs=length(FactorGroup1)
# Define the means and standard deviations for the three temperature groups
means = c(mean(FactorGroup1), mean(FactorGroup2), mean(FactorGroup3), mean(FactorGroup4))
sds = c(sd(FactorGroup1), sd(FactorGroup2), sd(FactorGroup3), sd(FactorGroup4))
# Estimate the power using the pwr.rm.anova function
pwr_rm_anova(means = means, sds = sds, n = gs, sig.level = 0.05, n.simulations = 100)
# Use the pwr.rm.anova function to find the required sample size for a power of 0.8
pwr_rm_anova(means = means, sds = sds, n = NULL, sig.level = 0.05, n.simulations = 100, desired_power = 0.8, start.n = 5)

#######################################################################################################
# Friedman Test (Non-Parametric Version of Randomized Block and Repeated Measures ANOVA)
#######################################################################################################

#Perform Test
friedman.test(DependentVariable ~ Factor | Block, data = DATA)

#Medians vs. Mean Ranks
A = dplyr::filter(DATA, Factor == "A")$DependentVariable
B = dplyr::filter(DATA, Factor == "B")$DependentVariable
C = dplyr::filter(DATA, Factor == "C")$DependentVariable
D = dplyr::filter(DATA, Factor == "D")$DependentVariable
# Plot boxplots for ratings of each restaurant
boxplot(A, B, C, D, names = c("A", "B", "C", "D"))
# Conduct Kolmogorov-Smirnov tests to compare ratings distributions between restaurants
ks1 = ks.test(A, B)$p
ks2 = ks.test(A, C)$p
ks3 = ks.test(A, D)$p
ks4 = ks.test(B, C)$p
ks5 = ks.test(B, D)$p
ks6 = ks.test(C, D)$p
pvalues = c(ks1, ks2, ks3, ks4, ks5, ks6) # Combine p-values from the tests
p.adjust(pvalues, method = 'BH')    # Adjust the p-values for multiple comparisons using the Benjamini-Hochberg method

#Follow Up Nemenyi Test
frdAllPairsNemenyiTest(DATA$DependentVariable, groups = DATA$Factor, blocks = DATA$Block)

#Power Analysis
library(fitdistrplus) # Load the library
x = A  # Replace with actual dataset
distributions = c("norm", "weibull", "gamma", "logis", "exp", "pois", "geom", "nbinom") # Distributions to fit
fit_list = list() # Fit distributions
for (dist in distributions) {
  fit = tryCatch({
    fitdist(x, dist)
  }, error = function(e) {
    cat("\nError in fitting", dist, ":", e$message, "\n")
    return(NULL)
  })
  if (!is.null(fit)) fit_list[[dist]] = fit
}
# Gather goodness-of-fit stats
fit_stats = data.frame(Distribution = character(), AIC = numeric(), BIC = numeric(), LogLikelihood = numeric(), Parameters = character(), stringsAsFactors = FALSE)
for (dist in names(fit_list)) {
  fit = fit_list[[dist]]
  params_str = paste(names(coef(fit)), "=", round(coef(fit), 4), collapse = ", ")
  fit_stats = rbind(fit_stats, data.frame(Distribution = dist, AIC = fit$aic, BIC = fit$bic, LogLikelihood = fit$loglik, Parameters = params_str))
}
print(fit_stats) # Display fit stats

#Post-Hoc
#Replace as needed
source("pwr_friedman_test.R")
gs= # Number of Blocks
posthoc_power = pwr_friedman_test(
  group_sizes = c(gs, gs, gs, gs),
  dist_list = list("norm", "norm", "norm", "norm"),
  params_list = list(list(mean = 77.5, sd = 4.2328),
                     list(mean = 66.6667, sd = 4.4222),
                     list(mean = 91, sd = 5.2599),
                     list(mean = 79.3333, sd = 4.4222)
  ),
  n.simulations = 1000, alpha = 0.05
)
print(posthoc_power)

#Apriori
#Replace as needed
sample_size = pwr_friedman_test(group_sizes = NULL, desired_power = 0.8,
                                dist_list = list("norm", "norm", "norm", "norm"),
                                params_list = list(list(mean = 77.5, sd = 4.2328),
                                                   list(mean = 66.6667, sd = 4.4222),
                                                   list(mean = 91, sd = 5.2599),
                                                   list(mean = 79.3333, sd = 4.4222)
                                ), n.simulations = 1000, alpha = 0.05, start_size = 2)
sample_size


##########################################################################################################
# Two-Way Repeated Measures ANOVA (1 between AND 1 Within)
##########################################################################################################

#Assumptions

#Normality

# Apply Shapiro-Wilk normality test within each group and adjust p-values
normality_results = DATA %>%
  group_by(Factor1, Factor2) %>%
  summarise(
    W = shapiro.test(DATA)$statistic,  # Extract W value
    p_value = shapiro.test(DATA)$p.value,  # Extract p-value
    .groups = 'drop'  # Drop grouping after summarise to avoid warning
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))  # Adjust p-values

normality_results

# Visual normality check with Histograms
DATA %>%
  ggplot(aes(x = DependentVariable)) +
  facet_wrap(~ Factor1 * Factor2) +
  geom_histogram(aes(y = ..density..), binwidth = 1, color = "black", fill = "lightblue") +
  stat_function(fun = dnorm, args = list(mean = mean(DATA$DependentVariable, na.rm = TRUE), 
                                         sd = sd(DATA$DependentVariable, na.rm = TRUE)), color = "red") +
  labs(x = "DependentVariable", y = "Density") +
  theme_minimal()

# Visual normality check with QQ plots
DATA %>%
  ggplot(aes(sample = DependentVariable)) +
  facet_grid(Factor1 ~ Factor2) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plots by Factor1 and Factor2")

#Outliers
# Plot a boxplot for 'DependentVariable' across different 'Factor1'
boxplot(DependentVariable ~ interaction(Factor1, Factor2), data = DATA,  ylab = "DependentVariable", cex.lab=0.5, cex.axis=0.5, cex=0.5)

# Identify and display outliers for each treatment group
outliers_df = DATA %>% 
  group_by(Factor1, Factor2) %>% 
  summarise(Outliers = list(boxplot.stats(DependentVariable)$out), .groups = 'drop') %>% 
  unnest(Outliers)
outliers_df

#Run ANOVA and test for Sphericity
two_way_result = aov_ez(
  id = "Subject",
  dv = "DependentVariable",
  within = "Factor1",         
  between = "Factor2",         
  data = DATA,
  anova_table = list(correction = "none", es = "ges") #If Mauchly's test is significant (p < 0.05) then use correction
)

# Display the ANOVA results
summary(two_way_result)

#Follow Up Test
# Follow-up tests for Factor1 (within-subjects factor)
pairwise_results_Factor1 = emmeans(two_way_result, pairwise ~ Factor1, adjust = "BH")
# Display pairwise comparisons for Factor1
pairwise_results_Factor1$contrasts


###################################################################################################
#ART (Aligned Rank Transformation) Model (Non-Parametric Version of Repeated Measures Two-Way ANOVA)
###################################################################################################

# Fit the ART model for the two-way mixed design
art_model = art(DependentVariable ~ Factor1 * Factor2 + Error(Subject/Factor1), data = DATA)
anova(art_model) # Summary of main effects and interaction
art.con(art_model, "Factor1") # Pairwise comparisons for within each group


####################################################################################################
#Pearson Correlation
####################################################################################################

#Assumptions

#Linear Relationship Between the Two Variables

#Plot x vs. y
plot(DATA$independentvariable, DATA$DependentVariable, xlab="IndependentVariable", ylab="DependentVariable", cex.axis=1.5, cex.lab=1.5, cex=1.5)
abline(lm(DATA$DependentVariable~DATA$IndependentVariable),col='red',lwd=2)

#Residuals vs. Fitted Line: Should see non-linear pattern
line = lm(DependentVariable~IndependentVariable, data=DATA) # Plot of residuals vs fitted values
plot(line, which=1, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) 


#Bivariate Normality
library(MVN) # Royston's test of bivariate normality
mvn(data=DATA, mvnTest = "royston", univariatePlot = "histogram", multivariatePlot = "contour")


#No Significant Outliers: Cook's Distance Plot
line = lm(DependentVariable~IndependentVariable, data=DATA)
plot(line, which=4, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Cut offs people use are 4/n, and 1. n: number of samples.
plot(line, which=5, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Points beyond dashed red curves are problematic


#Homoscedasticity
par(mar = c(4, 5, 2, 1)) 
line = lm(DependentVariable~IndependentVariable, data=DATA)
plot(line, which=3, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Spread Location Plot
#Should see randomly scattered points

library(lmtest)
bptest(line)    #perform Breusch-Pagan Test


#Perform the Correlation Test: Significance and Confidence Interval
cor.test(DATA$InpedendentVariable, DATA$DependentVariable) 
PCC=cor.test(DATA$IndependentVariable, DATA$DependentVariable)$estimate


#Power Analysis
library(pwr)
pwr.r.test(n = nrow(DATA), r = PCC, sig.level = 0.05, power = NULL, alternative = "two.sided") #Post-Hoc
pwr.r.test(n = NULL, r = PCC, sig.level = 0.05, power = 0.8, alternative = "two.sided") #Apriori


##################################################################################################
#Spearman Rank Correlation (Non-Parametric Version of Pearson's Correlation)
##################################################################################################

#Assumptions

#Monotonic Relationship Between Two Variables
plot(DATA$IndependentVariable, DATA$DependentVariable, xlab="IndependentVariable", ylab="DependentVariable",cex.axis=1.5, cex.lab=1.5, cex=1.5)

#Perform the Test
cor.test(DATA$IndependentVariable, DATA$DependentVariable, method="spearman", exact = TRUE) 
library(DescTools)
SpearmanRho(DATA$IndependentVariable, DATA$DependentVariable, conf=0.95)

#Power Analysis
source("pwr.spearman.test.R")
# Calculate the power
pwr.spearman.test(data_sample_size, true_rho, calculation_type = "power")
# Calculate the sample size for desired power
pwr.spearman.test(initial_sample_size, true_rho, desired_power, calculation_type = "sample_size")

###################################################################################################
#Simple Linear Regression
###################################################################################################

IndependentVariable = DATA$IndependentVariable
DependentVariable = DATA$DependentVariable

#Perform Correlation Test
plot(IndependentVariable, DependentVariable, xlab="IndependentVariable", ylab="DependentVariable",cex.axis=1.5, cex.lab=1.5, cex=1.5)
cor(IndependentVariable, DependentVariable)

#Assumptions

#Linear Relationship Between Variables: Plot of residuals, should be randomly scattered
plot(IndependentVariable, DependentVariable, xlab="IndependentVariable", ylab="DependentVariable", cex.axis=1.5, cex.lab=1.5, cex=1.5)
cor(IndependentVariable, DependentVariable)

#Calculate regression line for the next assumptions
line = lm(DependentVariable~IndependentVariable, data=DATA) # Run the lm() function, and save its results to line variable.

#Independent Residuals
set.seed(13)           # Set the seed for the random number generator for reproducibility.
library(car)           # Durbin-Watson test uses randomly generated data. 
durbinWatsonTest(line) 

#Normality of the Residuals
hist(line$residuals, freq=FALSE, xlab = "Residuals", main="",cex.axis=1.5,cex.lab=1.5)
curve(dnorm(x,mean(line$residuals), sd(line$residuals)), -1.0, 1.2, add=TRUE,col=c('red'))
qqnorm(line$residuals, main="")
qqline(line$residuals)
shapiro.test(line$residuals)

#Homoscedasticity
plot(line, which=1, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Plot of residuals vs fitted values
library(lmtest) #load lmtest library
bptest(line)    #perform Breusch-Pagan Test

#Outliers
plot(line, which=4, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Cut off we use is 4/n or 1.
plot(line, which=5, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Points beyond dashed red curves are problematic

#Perform Test
line = lm(DependentVariable~IndependentVariable, data=DATA) 
summary(line)                      # Provides the summary statistics for the regression line
confint(line,'IndependentVariable',level=0.95) # 95% Confidence Interval for the Slope

#How to perform prediction
newdata = data.frame(IndependentVariable=Prediction)
predict(line, newdata, interval="confidence")  # Our prediction, and its confidence interval

#Power Analysis
library(pwr)
r2  = summary(line)$r.squared # R-squared for our linear model
my.f2 = r2 / (1-r2)           # Effect Size
my.u = 1                      # u = 1 (for simple linear regression)
N=nrow(DATA)                  # Total sample size is N
my.v = N-my.u-1               # v = N - u - 1
pwr.f2.test(u = my.u, v = my.v, f2 = my.f2, sig.level = 0.05, power = NULL) # Power of the linear regression
pwr.f2.test(u = my.u, v = NULL, f2 = my.f2, sig.level = 0.05, power = 0.8)  # Total samples needed to get a power of 0.8.


##############################################################################################
#Siegel Regression (Non-Parametric Version of Simple Linear Regression)
##############################################################################################

#Run the Test
library(mblm)    # Load the 'mblm' package for nonparametric linear regression
model.s = mblm(DependentVariable ~ IndependentVariable, data = DATA) # Set 'repeated=FALSE' to use Theil-Sen method
summary(model.s) # Display a summary of the regression model

# Plot 'DependentVariable' against 'IndependentVariable'
plot(DATA$IndependentVariable, DATA$DependentVariable, xlab = "IndependentVariable", ylab = "DependentVariable", 
     pch = 16, cex.lab = 1.5, cex.axis = 1.5)
# Add the regression line from the Siegel model to the plot.
abline(model.s, col = "blue", lwd = 2, lty = 3)
# The legend indicates that the blue dotted line represents the Siegel regression line.
legend("topleft", legend = c("Siegel"), col = c("blue"), lty = 1:3, cex = 1.3)

#How to run prediction
# Create a new data frame for making predictions
newdata = data.frame(IndependentVariable = Value)
predict(model.s, newdata, interval = "confidence")

#Power Analysis
source("pwr.siegel.test.R") 
desired_power = 0.8               # Desired power level
n = length(IndependentVariable)   # Size of each simulated dataset
sd = sd(DependentVariable)        # Standard deviation of the noise
m = SLOPE                         # Slope from Siegel Regression Model
x_min = min(IndependentVariable)  # Minimum x
x_max = max(IndependentVariable)  # Maximum x

pwr.siegel.test(sample_size = n, x_min = x_min, x_max = x_max, noise_level = sd, true_slope = m, n_simulations = 1000) #Post-Hoc
pwr.siegel.test(desired_power = 0.8, x_min = x_min, x_max = x_max, noise_level = sd, true_slope = m, n_simulations = 1000, starting_sample_size = 5) #Apriori


###############################################################################################
#Multivariable Linear Regression
###############################################################################################

#Assumptions

#Linearity Between Dependent Variable and Independent Variables & Multicolinearity
# Create multiple scatterplots
pairs(~DependentVariable + IndependentVariable1 + IndependentVariable2 + IndependentVariable3, data = DATA, main = 'DependentVariable scatterplots') #Keep adding variables as needed
# Create a correlation matrix from selected columns and round it to two decimals
round(cor(DATA[c("DependentVariable", "IndependentVariable1", "IndependentVariable2", "IndependentVariable3")]), 2)

#Should be roughly linear relationship between the dependent variable and each of the independent variables
#Should NOT be a correlation of 0.8 or higher for any pair of independent variables. If there is then remove one. 

#VIF: another measure of multicolinearity
library(car)
MulLinReg=lm(DependentVariable ~ IndependentVariable1 + IndependentVariable2 + IndependentVariable3, data=DATA)
vif(MulLinReg)

#Shouldn't have values > 10

#Independent Residuals: Should have Durbin Watson test statistic between 1.5 and 2.5
set.seed(13)
library(car)
durbinWatsonTest(MulLinReg)

#Homoscedasticity and Linearity of Residuals: Breusch-Pagan Test
set.seed(13)
library(car)
durbinWatsonTest(MulLinReg)

#Outliers
plot(MulLinReg, which = 4)
plot(MulLinReg, which = 5)

#Normality of Residuals
hist(MulLinReg$residuals, freq=FALSE, xlab = "Residuals", main="",cex.axis=1.5,cex.lab=1.5)
curve(dnorm(x,mean(MulLinReg$residuals), sd(MulLinReg$residuals)), -1.6, 2.1, add=TRUE,col=c('red'))
qqnorm(MulLinReg$residuals, main="")
qqline(MulLinReg$residuals)
shapiro.test(MulLinReg$residuals)

#Perform Multivariable Linear Regression
summary(MulLinReg)

#How to perform prediction if dataset is given
predicted_values = predict(MulLinReg, newdata = newdata) # Predict values using the regression model MulLinReg
actual_values = newdata$DependentVariable # Actual values from the test dataset

SS_res = sum((actual_values - predicted_values)^2) # Residual sum of squares
SS_tot = sum((actual_values - mean(actual_values))^2) # Total sum of squares
R_squared = 1 - (SS_res / SS_tot) # Calculate R-squared

print("Predicted Values:") # Print the predictions and evaluation metrics
print(predicted_values)
print(paste("R-squared:", R_squared))

#How to perform prediction on one point
newdata=data.frame(IndependentVariables)
predict(MulLinReg, newdata, interval="confidence")  # Our prediction, and its confidence interval

#Power Analysis
library(pwr)
r2  = summary(MulLinReg)$r.squared          # R-squared for our linear model
my.f2 = r2 / (1-r2)                         # Effect Size
my.u = length(MulLinReg$coefficients)-1     # number of independent variables
N=nrow(DATA)                                # Total sample size is N
my.v = N-my.u-1                             # v = N - u - 1

# POSTHOC POWER ANALYSIS: POWER OF THE REGRESSION ANALYSIS
pwr.f2.test(u = my.u, v = my.v, f2 = my.f2, sig.level = 0.05, power = NULL) 

# A-PRIORI POWER ANALYSIS: SAMPLE SIZE
ceiling(pwr.f2.test(u = my.u, v = NULL, f2 = my.f2, sig.level = 0.05, power = 0.8)$v + my.u + 1)

#How to perform stepwise regression
library(MASS)
step.model = stepAIC(MulLinReg, direction = "both", trace = FALSE) # Stepwise regression model
summary(step.model)


###############################################################################################
# Multivariable Binary Logistic Regression
###############################################################################################

#Shit ton of libraries
library(mlbench)     # For PimaIndiansDiabetes2 dataset
library(dplyr)       # For data manipulation (dplyr) 
library(broom)       # For making model summary tidy
library(visreg)      # For plotting logodds and probability 
library(rcompanion)  # To calculate pseudo R-squared
library(MASS)        # For stepwise model selection
library(ROCR)        # To find the probability threshold for best accuracy
library(car)         # For multicollinearity function vif()

#Fit the Model (entire dataset)
model_logi = glm(DependentVariable~., data=DATA, family = "binomial")      # Fitting a binary logistic regression
summary(model_logi)                   # Model summary

#Assumptions

#No multicolinearity
attach(DATA) # COULD ACCESS INDIVIDUAL COLUMNS WITHOUT WRITING: Diabetes$pedigree
round(cor(cbind(independentvariables)),2) # Note that we are not checking correlation for categorical variables.
vif(model) # VIF assesses the relationships between each independent variable and all the other variables.

#Outliers
plot(model_logi, which=4, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Cut offs we use are 0.5, and 1. n: number of samples.
plot(model_logi, which=5, cex.lab=1.3,cex.axis=1.3, cex=1.3, cex.id = 1.3, cex.caption = 1.3) # Points beyond dashed red curves are problematic

#Linear relationship between odds and independent variable
pred = predict(model_logi, DATA, type="response") 
logodds=logit(pred)
pairs(~logodds + independentvariables, data=DATA)



#How to perform stepwise regression
step.model = stepAIC(model_logi, direction = "both", trace = FALSE)   # Stepwise regression model
summary(step.model)

#Odds ratio
tidy(model, exponentiate = TRUE, conf.level = 0.95) # Odds ratio table

#NagelKerke Model Fitness
nagelkerke(model) # Pseudo R_squared values and Likelyhood ratio test

#Baseline test of model accuracy on new data
pred = predict(model, newdata, type="response")  # Predict the test dataset
predicted = round(pred)  # Round of the value; >0.5 will convert to 1 else 0
predicted
tab = table(Predicted = predicted, Reference = newdata$dependentvariable) # Creating a contingency table
tab
AllTestSamples = tab[1,1] + tab[1,2] + tab[2,1] + tab[2,2]
Accuracy = (tab[1,1]+tab[2,2]) / AllTestSamples     # Accuracy=(Correct Predictions)/(All Test Samples)
Accuracy         
BaselineAccuracy = max(tab[1,1] + tab[2,1], tab[1,2] + tab[2,2]) / AllTestSamples
BaselineAccuracy # BaselineAccuracy = (Number of Samples in the Larger Category)/(All Test Samples)



#Power Analysis

#preparation
null_model = glm(dependentvariable ~ 1, data = DATA, family = binomial) # Define the null model (intercept-only logistic regression model)
pseudo_r2 = nagelkerke(model_logi)$Pseudo.R.squared[1]  # Calculate McFadden's pseudo-R^2
effect_size = pseudo_r2 / (1 - pseudo_r2) # Calculate effect size (f^2)
num_predictors = length(model_logi$coefficients) - 1  # Number of predictors (u)
total_sample_size = nrow(DATA) # Total sample size
residual_dof = total_sample_size - num_predictors - 1 # Degrees of freedom for residuals (v)

#Post Hoc
posthoc_power = pwr.f2.test(
  u = num_predictors, 
  v = residual_dof, 
  f2 = effect_size, 
  sig.level = 0.05, 
  power = NULL
)
print(posthoc_power)

#Apriori
required_sample_size = ceiling(pwr.f2.test(
  u = num_predictors, 
  v = NULL, 
  f2 = effect_size, 
  sig.level = 0.05, 
  power = 0.8
)$v + num_predictors + 1)
print(paste("Required Sample Size:", required_sample_size))


#########################################################################################################
# Non-linear Regression
#########################################################################################################

# Create a plot of the data
plot(independentvariable, DependentVariable, 
     xlab = "independent variable",        # Label for x-axis
     ylab = "dependent variable",  # Label for y-axis
     cex.lab = 1.5,   # Increase size of axis labels
     cex.axis = 1.5,  # Increase size of axis tick marks
     cex = 1.5)       # Increase size of plot points


#Transform the data and plot it
transx = #transform independent variable 
transy = #transform dependent variable if needed
TransformedDF = data.frame(transx, transy)

# Perform linear regression using the transformed variables
line = lm(transy ~ transx, data = TransformedDF)

# Plot the transformed data
plot(transx, transy, 
     pch = 1,                  # Point character
     cex.lab = 1.5,            # Size of axis labels
     cex.axis = 1.5,           # Size of axis tick marks
     cex = 1.5)                # Size of plot points
# Add the regression line to the plot
abline(line, col = 'red', lwd = 2, lty = 1)

# Fit a linear regression model using the transformed dependent variable as the response variable and the transformed independent variable as the predictor variable
line = lm(transy ~ transx, data = TransformedDF)
# Display a summary of the linear regression model, including coefficients, R-squared value, and other statistics
summary(line)
# Calculate and display the 95% confidence interval for the slope coefficient 
confint(line, 'transx', level = 0.95)

#Plot the original data with the new regression line
c = line$coefficients[1]      # Intercept of the regression line
d = line$coefficients[2]      # Slope of the regression line
a= #transform back from earlier transformation
b #transform back from earlier transformation
plot(independentvariable, DependentVariable, xlab='independent variable',ylab='dependent variable', cex.axis=1.5,cex.lab=1.5, cex=2.0)
curve(formula, from=min(IndependentVariable), to=max(IndependentVariable), col="red", lwd=2, add=TRUE)
legend('bottomright', legend='formula', lwd=3, col=c("red"), lty=1, cex=1.5)

#How to make prediction with new data

# Transform the value you want to predict
newdata = data.frame(TransSubsConc = transformedvalue)

#calculate prediction and confidence interval
predicted_confidence_interval = predict(line, newdata, interval = "confidence")

# Extract the fitted value from the prediction and back-transform it
# Applying the back-transformation y = 1/ytilde
LinearFit = predicted_confidence_interval[,"fit"]
Fit = #transform back

# Extract and back-transform the lower bound of the confidence interval
LinearLB = predicted_confidence_interval[,"lwr"]
UB = #transform back

# Extract and back-transform the upper bound of the confidence interval
LinearUB = predicted_confidence_interval[,"upr"]
LB = #transform back

# Output the back-transformed fitted value and confidence interval bounds
cat("Fitted Value:", Fit, "\n")
cat("Lower Bound of Confidence Interval:", LB, "\n")
cat("Upper Bound of Confidence Interval:", UB, "\n")



#########################################################################################
# Linear Mixed Effects Model
#########################################################################################

#Random Intercepts

# Fit the Linear Mixed Effects Model
model = lmer(dependentvariable ~ independentvariables, data = DATA)
summary(model)

#Model Significance
model.null = lmer(dependentvariable ~ random_independent_variables, data = DATA)
anova(model, mnodel.null)

#Assumptions

#Linearity: Residuals vs Fitted Plot
ggplot(data = NULL, aes(x = fitted(model), y = resid(model))) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

#Collinearity
performance::check_collinearity(model) #Should be < 10


# Residual Normality
qqnorm(resid(model))
qqline(resid(model), col = "blue")
shapiro.test(resid(model))

#Homoscedasticity:
performance::check_heteroscedasticity(model)

#Independence of residuals
# Residual Autocorrelation Plot
acf(resid(model), main = "Autocorrelation of Residuals")
library(randtests)
runs.test(residuals(model))

#Outliers should be < 1
influence_results = influence(model, group = "subject")
plot(influence_results, which = "cook")

#Power analysis
library(simr)
# Linear Mixed Effects Model
model = lmer(dependentvariable ~ independentvariables, data = DATA)

# The number of subjects
extended_model = extend(model, along = "subject", n = 15) 

# Power analysis for fixed effect
powerSim(extended_model, fixed("fixed effect"), nsim = 100)


##############################################################################################
#Independence & Goodness of Fit Tests (Fisher's Exact, Chi-Squared, etc)
###############################################################################################

#Binomial: Only works for 2x2
binom.test(observations,trials,probability) # THIS IS ONE SIDED. NOTE: A binomial test is always 1-sided unless p = 0.5


#Chi-Square Test for goodness of it
obs.counts = c()
exp.probs = c()                             # Probability should sum to 1. 
exp.counts = exp.probs*sum(obs.counts)       # Expected counts
exp.counts

#Assumptions: 
# (i) Levels are exclusive and observations are independent.
# (ii) The sample size is large. (You can do Williams correction if the sample size is small.)
# (iii) Expected frequencies for each cell are at least 1
# (iv) More than 80% of cells have an expected frequency 5 or more.

#Run the test
chisq.test(obs.counts, p = exp.probs)                  # Chi-squared test of goodness-of-fit

#Follow up pairwise comparison
source("PairwiseChiSquareGoodnessOfFit.R")
PairwiseChiSquareGoodnessOfFit(obs.counts, exp.probs)


#Power Analysis
library(pwr)
obs.probs = c()/sum(c())
exp.probs = c()
effect.size = ES.w1(obs.probs, exp.probs)
degrees = length(obs.probs) - 1

# Posthoc Power Analysis: Number of samples needed to achieve a power of 0.8
pwr.chisq.test(w=effect.size,
               N=sum(obs.counts),   # Total number of observations
               df=degrees,          # Degrees of Freedom
               power=NULL,          # 1 minus Type II probability
               sig.level=0.05)      # Type I probability

# Apriori Power Analysis: Number of samples needed to achieve a power of 0.8
pwr.chisq.test(w=effect.size,
               N=NULL,   # Total number of observations
               df=degrees,          # Degrees of Freedom
               power=0.8,          # 1 minus Type II probability
               sig.level=0.05)      # Type I probability


#Fisher's Exact test of independence for 2x2 tables
my.table = matrix(c(),nrow=2)
row.names(my.table) = c()
colnames(my.table) = c()
my.table

library(DescTools)
fisher.test(my.table, workspace = 2e8)  # Fisher Exact Test is preferred as long as you can calculate it.

#Chi-Squared test for independence
my.table = matrix(c(), nrow=)
row.names(my.table) = c()
colnames(my.table) = c()
my.table

my.chisq = chisq.test(my.table)  # We give a table to Chisq function
my.chisq
my.chisq$expected                # Chi-square test assumptions pass
# If it fails, Williams correction is possible.


#Follow up pairwise comparison
source("PairwiseChiSquareTestOfIndependence.R")
PairwiseChiSquareTestOfIndependence(my.table, comparison_type = "row")
PairwiseChiSquareTestOfIndependence(my.table, comparison_type = "column")

#Power analysis
library(pwr)
# Calculate observed and expected probabilities
obs.probs = prop.table(my.table)  # Row-wise proportions
exp.probs = prop.table(my.chisq$expected)

# Calculate effect size
effect.size = ES.w1(as.vector(obs.probs), as.vector(exp.probs))

# Degrees of freedom for a chi-square test of independence
degrees = (nrow(my.table) - 1) * (ncol(my.table) - 1)

# Perform power analysis for the chi-square test
pwr.chisq.test(w = effect.size,
               N = sum(my.table),  # Total number of observations
               df = degrees,         # Degrees of Freedom
               power = NULL,         # 1 minus Type II probability
               sig.level = 0.05)     # Type I probability

# Perform power analysis for the chi-square test
pwr.chisq.test(w = effect.size,
               N = NULL,  # Total number of observations
               df = degrees,         # Degrees of Freedom
               power = 0.8,         # 1 minus Type II probability
               sig.level = 0.05)     # Type I probability


#######################################################################################
# Odds Ratio and Relative Risk
#######################################################################################

#Relative Risk
factor1 = c("","")
factor2 = c("", "")
data = matrix(c(,,,), nrow = 2)
dimnames(data) = list("factor1" = factor1, "factor2" = factor2)
RelativeRisk = riskratio(data)
RelativeRisk


#Odds Ratio
library(epitools)
factor1 = c("","")
factor2 = c("", "")
data = matrix(c(,,,), nrow = 2)
dimnames(data) = list("factor1" = factor1, "factor2" = factor2)
OddsRatio = oddsratio(data)
OddsRatio





