# load package
library(lme4) # for linear mixed models (not survival)
library(sjPlot)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(mediation)
library(ggplot2)
library(lmeresampler)


## load data
#setwd("/export04/data/mlevin/Uniformity_Illusion/bhv_data")
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
df <- read.csv('./minrates.csv',sep = ",",stringsAsFactors = FALSE)
df$FT = df$RT
df <- df %>%
  convert_as_factor(subject)

## ------- mixed effect models on all data ------- ##

## fit minimum rate
# BEST: lm6b: (1 + baselinerate | subject) + contrast + scaled_eccentricity + baselinerate
lm0 = lmer(msratemin ~ (1 | subject), df)
lm1 = lmer(msratemin ~ (1 | subject) + contrast, df)
lm2 = lmer(msratemin ~ (1 | subject) + eccentricity, df)
lm2b = lmer(msratemin ~ (1 | subject) + scaled_eccentricity, df)
lm3 = lmer(msratemin ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(msratemin ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(msratemin ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(msratemin ~ (1 | subject) + contrast*scaled_eccentricity, df)
lm5 = lmer(msratemin ~ (1 | subject) + contrast + eccentricity + baselinerate, df)
lm5b = lmer(msratemin ~ (1 | subject) + contrast + scaled_eccentricity + baselinerate, df)
lm5c = lmer(msratemin ~ (1 + baselinerate | subject) + contrast + baselinerate, df)
lm6 = lmer(msratemin ~ (1 + baselinerate | subject) + contrast + eccentricity + baselinerate, df)
lm6b = lmer(msratemin ~ (1 + baselinerate | subject) + contrast + scaled_eccentricity + baselinerate, df)
lm7 = lmer(msratemin ~ (1 + baselinerate | subject) + contrast*eccentricity + baselinerate, df)
lm7b = lmer(msratemin ~ (1 + baselinerate | subject) + contrast*scaled_eccentricity + baselinerate , df)
lm8 = lmer(msratemin ~ (1 | subject) + baselinerate*contrast + eccentricity, df)
lm8b = lmer(msratemin ~ (1 | subject) + baselinerate*contrast + scaled_eccentricity, df)
lm9 = lmer(msratemin ~ (1 | subject) + baselinerate*eccentricity + contrast, df)
lm9b = lmer(msratemin ~ (1 | subject) + baselinerate*scaled_eccentricity + contrast, df)
lm10 = lmer(msratemin ~ (1 + baselinerate | subject) + baselinerate*eccentricity + contrast, df)
lm10b = lmer(msratemin ~ (1 + baselinerate | subject) + baselinerate*scaled_eccentricity + contrast, df)
lm11 = lmer(msratemin ~ (1 + baselinerate | subject) + baselinerate*contrast + eccentricity, df)
lm11b = lmer(msratemin ~ (1 + baselinerate | subject) + baselinerate*contrast + scaled_eccentricity, df)
lm12 = lmer(msratemin ~ (1 | subject) + baselinerate*contrast*eccentricity, df)
lm12b = lmer(msratemin ~ (1 | subject) + baselinerate*contrast*scaled_eccentricity, df)
lm20 = lmer(msratemin ~ (1 + baselinerate | subject) + baselinerate*contrast*eccentricity, df)
lm21 = lmer(msratemin ~ (1 + baselinerate+contrast | subject) + baselinerate*contrast*eccentricity, df)
anova(lm0, lm1, lm2, lm2b, lm3, lm3b, lm4, lm4b, lm5, lm5b, lm6, lm6b, lm7, lm7b, lm8, lm8b, lm9, lm9b, lm10, lm10b, lm11, lm11b, test = "Chisq")
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11) # only reg eccentricity
anova(lm0, lm1, lm2b, lm3b, lm4b, lm5b, lm6b, lm7b, lm8b, lm9b, lm10b, lm11b) # only scaled_ecc


# calculate Bayes factor (evidence for 2nd)
exp((BIC(lm6) - BIC(lm5c))/2)


## fit same mixed models to amplitude
# pretty much best: just baselineampl predicts
lm0 = lmer(amplmin ~ (1 | subject), df)
lm0b = lmer(amplmin ~ (1 | subject) + baselineampl, df)
lm1 = lmer(amplmin ~ (1 | subject) + contrast, df)
lm2 = lmer(amplmin ~ (1 | subject) + eccentricity, df)
lm2b = lmer(amplmin ~ (1 | subject) + scaled_eccentricity, df)
lm3 = lmer(amplmin ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(amplmin ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(amplmin ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(amplmin ~ (1 | subject) + contrast*scaled_eccentricity, df)
lm5 = lmer(amplmin ~ (1 | subject) + contrast + eccentricity + baselineampl, df)
lm5b = lmer(amplmin ~ (1 | subject) + contrast + scaled_eccentricity + baselineampl, df)
lm6 = lmer(amplmin ~ (1 + baselineampl | subject) + contrast + eccentricity + baselineampl, df)
lm6b = lmer(amplmin ~ (1 + baselineampl | subject) + contrast + scaled_eccentricity + baselineampl, df)
lm7 = lmer(amplmin ~ (1 + baselineampl | subject) + contrast*eccentricity + baselineampl, df)
lm7b = lmer(amplmin ~ (1 + baselineampl | subject) + contrast*scaled_eccentricity + baselineampl , df)
lm8 = lmer(amplmin ~ (1 | subject) + baselineampl*contrast + eccentricity, df)
lm8b = lmer(amplmin ~ (1 | subject) + baselineampl*contrast + scaled_eccentricity, df)
lm9 = lmer(amplmin ~ (1 | subject) + baselineampl*eccentricity + contrast, df)
lm9b = lmer(amplmin ~ (1 | subject) + baselineampl*scaled_eccentricity + contrast, df)
lm10 = lmer(amplmin ~ (1 + baselineampl | subject) + baselineampl*eccentricity + contrast, df)
lm10b = lmer(amplmin ~ (1 + baselineampl | subject) + baselineampl*scaled_eccentricity + contrast, df)
lm11 = lmer(amplmin ~ (1 | subject) + baselineampl*contrast*eccentricity, df)
lm11b = lmer(amplmin ~ (1 | subject) + baselineampl*contrast*scaled_eccentricity, df)


anova(lm0, lm0b, lm1, lm2, lm2b, lm3, lm3b, lm4, lm4b, lm5, lm5b, lm6, lm6b, lm7, lm7b, lm8, lm8b, lm9, lm9b, lm10, lm10b, lm11, lm11b, test = "Chisq")
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11) # only reg eccentricity
anova(lm0, lm1, lm2b, lm3b, lm4b, lm5b, lm6b, lm7b, lm8b, lm9b, lm10b, lm11b) # only scaled_ecc



## model baseline rate
# BEST: just subject intercept (lm0)
# actually if you control for FT, eccentricity is significant now.
# also include eccentricity/contrast as random subject effect, because not all people
# might do this (only if model converges!)
#     -> eccentricity as random converges, contrast does not.
#     so using better models: lm2 is best (same as single trial analyses below)
lm0 = lmer(baselinerate ~ (1 | subject), df)
lm1 = lmer(baselinerate ~ (1 | subject) + contrast, df)
lm2 = lmer(baselinerate ~ (1 | subject) + eccentricity, df)
lm2b = lmer(baselinerate ~ (1 | subject) + scaled_eccentricity, df)
lm3 = lmer(baselinerate ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(baselinerate ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(baselinerate ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(baselinerate ~ (1 | subject) + contrast*scaled_eccentricity, df)
lm5 = lmer(baselinerate ~ (1 | subject) + contrast + eccentricity + msratemin, df) # not appropriate I think

anova(lm0, lm1, lm2, lm2b, lm3, lm3b, lm4, lm4b)

lm0 = lmer(baselinerate ~ (1 + FT | subject) + FT, df)
lm1 = lmer(baselinerate ~ (1 + FT | subject) + contrast + FT, df)
lm2 = lmer(baselinerate ~ (1 + FT + eccentricity | subject) + eccentricity + FT, df)
lm3 = lmer(baselinerate ~ (1 + FT | subject) + contrast + eccentricity + FT, df)
lm4 = lmer(baselinerate ~ (1 + FT | subject) + contrast*eccentricity + FT, df)

anova(lm0, lm1, lm2, lm3, lm4)


## model baseline blink rate - best: no predictors (lm0).
lm0 = lmer(baselineblinks ~ (1 | subject) + FT, df)
lm1 = lmer(baselineblinks ~ (1 | subject) + contrast, df)
lm2 = lmer(baselineblinks ~ (1 | subject) + eccentricity, df)
lm2b = lmer(baselineblinks ~ (1 | subject) + scaled_eccentricity, df)
lm3 = lmer(baselineblinks ~ (1 + FT | subject) + FT + contrast + eccentricity, df)
lm3b = lmer(baselineblinks ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(baselineblinks ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(baselineblinks ~ (1 | subject) + contrast*scaled_eccentricity, df)

anova(lm0, lm1, lm2, lm2b, lm3, lm3b, lm4, lm4b)

##  model mean FT
# BEST: lm5 (1 | subject) + contrast + eccentricity + baselinerate
# eccentricity*baselinerate is significant when included (lm10), but that model might have higher BIC than without the interaction
# seems almost equivalent BIC...
lm0 = lmer(FT ~ (1 | subject), df)
lm0b = lmer(FT ~ (1 | subject) + baselinerate, df)
lm1 = lmer(FT ~ (1 | subject) + contrast, df)
lm2 = lmer(FT ~ (1 | subject) + eccentricity, df)
lm2b = lmer(FT ~ (1 | subject) + scaled_eccentricity, df)
lm3 = lmer(FT ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(FT ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(FT ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(FT ~ (1 | subject) + contrast*scaled_eccentricity, df)
lm5 = lmer(FT ~ (1 | subject) + contrast + eccentricity + baselinerate, df)
lm5b = lmer(FT ~ (1 | subject) + contrast + scaled_eccentricity + baselinerate, df)
lm5c = lmer(FT ~ (1 | subject) + contrast + scaled_eccentricity + baselinerate + baselineblinks, df)
lm6 = lmer(FT ~ (1 | subject) + contrast*eccentricity + baselinerate, df)
lm6b = lmer(FT ~ (1 | subject) + contrast*scaled_eccentricity + baselinerate , df)
lm7 = lmer(FT ~ (1 + baselinerate | subject) + contrast + eccentricity + baselinerate, df)
lm7b = lmer(FT ~ (1 + baselinerate | subject) + contrast + scaled_eccentricity + baselinerate, df)
lm8 = lmer(FT ~ (1 + baselinerate | subject) + contrast*eccentricity + baselinerate, df)
lm8b = lmer(FT ~ (1 + baselinerate | subject) + contrast*scaled_eccentricity + baselinerate , df)
lm9 = lmer(FT ~ (1 | subject) + baselinerate*contrast + eccentricity, df)
lm9b = lmer(FT ~ (1 | subject) + baselinerate*contrast + scaled_eccentricity, df)
lm10 = lmer(FT ~ (1 | subject) + baselinerate*eccentricity + contrast, df)
lm10b = lmer(FT ~ (1 | subject) + baselinerate*scaled_eccentricity + contrast, df)
lm11 = lmer(FT ~ (1 + baselinerate | subject) + baselinerate*contrast + eccentricity, df)
lm11b = lmer(FT ~ (1 + baselinerate | subject) + baselinerate*contrast + scaled_eccentricity, df)
lm12 = lmer(FT ~ (1 + baselinerate | subject) + baselinerate*eccentricity + contrast, df)
lm12b = lmer(FT ~ (1 + baselinerate | subject) + baselinerate*scaled_eccentricity + contrast, df)
lm13 = lmer(FT ~ (1 | subject) + baselinerate*contrast*eccentricity, df)
lm13b = lmer(FT ~ (1 | subject) + baselinerate*contrast*scaled_eccentricity, df)
lm14 = lmer(FT ~ (1 + baselinerate | subject) + baselinerate*contrast*eccentricity, df)
lm14b = lmer(FT ~ (1 + baselinerate | subject) + baselinerate*contrast*scaled_eccentricity, df)

anova(lm0, lm0b, lm1, lm2, lm2b, lm3, lm3b, lm4, lm4b, lm5, lm5b, lm6, lm6b, lm7, lm7b, lm8, lm8b, lm9, lm9b, lm10, lm10b, lm11, lm11b, lm12, lm12b, lm13, lm13b, lm14, lm14b)
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11) # only reg eccentricity
anova(lm0, lm1, lm2b, lm3b, lm4b, lm5b, lm6b, lm7b, lm8b, lm9b, lm10b, lm11b, lm12b, lm13b, lm14b) # only scaled_ecc

exp((BIC(lm5) - BIC(lm10))/2)


## anova and boxplots
df2 = df

# regress out baselinerate & subjecteffects?
lmrt = lmer(FT ~ (1 | subject) + baselinerate, df2) # remove subj mean, baselinerate from FT
lmminrate = lmer(msratemin ~ (1 + baselinerate | subject) + baselinerate, df2) # remove subj mean, baselinerate from msratemin
lmbaselinerate = lmer(baselinerate ~ (1 | subject), df2)

df2$FT = residuals(lmrt) + mean(df$FT) # center around the total population mean FT
df2$msratemin = residuals(lmminrate) + mean(df$msratemin)
df2$baselinerate = residuals(lmbaselinerate) + mean(df$baselinerate)

# regress out contrast or eccentricity?
lmrtnocon = lm(FT ~ 1 + contrast, df2) # remove contrast from FT
lmminratenocon = lm(msratemin ~ 1 + contrast, df2) # remove contrast from minrate
lmbaselineratenocon = lm(baselinerate ~ 1 + contrast, df2) # remove contrast from baselinerate
lmrtnoecc = lm(FT ~ 1 + eccentricity, df2) # remove contrast from FT
lmminratenoecc = lm(msratemin ~ 1 + eccentricity, df2) # remove contrast from minrate
lmbaselineratenoecc = lm(baselinerate ~ 1 + eccentricity, df2) # remove contrast from baselinerate

df2$FT = residuals(lmrtnocon) + mean(df$FT)
df2$msratemin = residuals(lmminratenocon) + mean(df$msratemin)
df2$baselinerate = residuals(lmbaselineratenocon) + mean(df$baselinerate)

# convert for plotting
df2 <- df2 %>%
  convert_as_factor(subject, contrast, eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))

# summary stats
df2 %>%
  group_by(eccentricity) %>%
  get_summary_stats(msratemin, type = "mean_sd")

# boxplot of mean FTs / msratemin

plottitle="all trials"
ydata = "FT"; ylabel='FT (s)'; ylimits=c(0, 15)
#ydata="msratemin"; ylabel='minimum rate (Hz)'; ylimits=c(0, 2)
#ydata="baselinerate"; ylabel='baseline rate (Hz)'; ylimits=c(0, 2)

bxp <- ggboxplot(
  df2, x = "contrast", y = ydata,
  color = "eccentricity", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab(ylabel) + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = ylimits)


# boxplot of mean FTs / msratemin, averaged within conditions
bxp <- ggboxplot(
  df2, x = "contrast", y = ydata,
  color = "contrast", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab(ylabel) + theme(text = element_text(size = 20))
bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) # for eccentricity
#bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = ylimits) # normalized to pop. mean







## -------- comparing no-microsaccade trials to microsaccade trials -------- ##

# noms has 853 (1152) trials, onlyms has 5314 (5015).
#    (parentheses = ignore t=0-750ms when selecting trials)
#    (we shouldn't do that because even the 1st microsaccade can produce an eccentricity effect)
# onlyms_subset_minrates.csv is a random sample of ~853 (~1152) trials from onlyms data.
## load data
#df <- read.csv('./noms_minrates.csv',sep = ",",stringsAsFactors = FALSE)
df <- read.csv('./onlyms_minrates.csv',sep = ",",stringsAsFactors = FALSE)
#df <- read.csv('./onlyms_subset_minrates.csv',sep = ",",stringsAsFactors = FALSE)
df$FT = df$RT
df <- df %>%
  convert_as_factor(subject)

## predict FT from just contrast and eccentricity (compare ms to noms)
# Result: noms: best model is 1: only contrast predicts.
# onlyms: best model is 3/3b: both contrast & eccentricity predict.
# even if we use small onlyms subset to match trial # in noms!!
lm0 = lmer(FT ~ (1 | subject), df)
lm1 = lmer(FT ~ (1 | subject) + contrast, df)
lm2 = lmer(FT ~ (1 | subject) + eccentricity, df)
lm2b = lmer(FT ~ (1 | subject) + scaled_eccentricity, df)
lm3 = lmer(FT ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(FT ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(FT ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(FT ~ (1 | subject) + contrast*scaled_eccentricity, df)

lm9 = lmer(FT ~ (1 | subject) + baselinerate + eccentricity+ contrast, df)
lm10 = lmer(FT ~ (1 | subject) + baselinerate*eccentricity + contrast, df)

anova(lm0, lm1, lm2, lm2b, lm3, lm3b, lm4, lm4b, lm10, test = "Chisq")


## anova and boxplots
df2 <- df %>%
  convert_as_factor(subject, contrast, eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))

# regress out subject effect?
lmrt = lmer(FT ~ (1 | subject), df2) # remove subj mean from FT
df2$FT = residuals(lmrt) + mean(df$FT) # center around the total population mean FT

# summary stats
df2 %>%
  group_by(contrast) %>%
  get_summary_stats(FT, type = "mean_sd")


# boxplot of mean FTs
#plottitle="at least one microsaccade"
plottitle="no microsaccades or blinks"
bxp <- ggboxplot(
  df2, x = "contrast", y = "FT",
  color = "eccentricity", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15))

# boxplot of mean FTs, averaged within conditions
bxp <- ggboxplot(
  df2, x = "contrast", y = "FT",
  color = "contrast", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15)) # normalized to pop. mean FT






## ---- individual trials ---- ##
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
#setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/Uniformity_Illusion/eye_analysis")
df <- read.csv('./trial_data.csv',sep = ",",stringsAsFactors = FALSE)
df$FT = df$RT
df$baselinemsrate = df$baselinenms / df$FT
# subset of subjects?
#df = df[(df$subject < 12),]
#df = df[(df$subject > 11),]
df <- df %>%
  convert_as_factor(subject, mspresent, earlymspresent, latemspresent, blinkpresent, session)
dforig = df
# add scaled_eccentricity
#A = 0.063
#l0 = 36.54 * A
#df$scaled_eccentricity = (log(df$eccentricity) - l0) / A
#df$scaled_eccentricity = 1 / (A * df$eccentricity)

# select only trials w/o ms or blinks:
df = dforig[(dforig$baselinenms == 0) & (dforig$nblinks == 0),]
df = dforig[(dforig$baselinenms >= 0) & (dforig$nblinks == 0),]

# convert baselinenms and nblinks to rate (divide by FT)
dforig = df
df$baselinenms = dforig$baselinenms / dforig$FT
df$nblinks = dforig$nblinks / dforig$FT
# and set baselineretinalslip minimum to 0
df$baselineretinalslip = dforig$baselineretinalslip - min(dforig$baselineretinalslip, na.rm = TRUE)

# FT: fit mixed effect models LM20 BEST
# baselineretinalslip is NOT significant and does not improve model fit?
lm0 = lmer(FT ~ (1 | subject), df)
lm1 = lmer(FT ~ (1 | subject) + contrast, df)
lm2 = lmer(FT ~ (1 | subject) + eccentricity, df)
lm3 = lmer(FT ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(FT ~ (1 | subject) + contrast + eccentricity + eccentricity*baselineretinalslip, df)
lm4 = lmer(FT ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(FT ~ (1 | subject) + contrast*eccentricity*baselineretinalslip, df)
lm5 = lmer(FT ~ (1 | subject) + contrast + eccentricity + baselinenms, df)
lm6 = lmer(FT ~ (1 | subject) + contrast + eccentricity*baselinenms, df)
lm7 = lmer(FT ~ (1 | subject) + contrast*baselinenms + eccentricity, df)
lm8 = lmer(FT ~ (1 | subject) + contrast*baselinenms + eccentricity*baselinenms, df)

lm7 = lmer(FT ~ (1 | subject) + contrast + eccentricity + baselinenms + contrast:baselinenms+ eccentricity:baselinenms, df)
lm7b = lmer(FT ~ (1 | subject) + contrast + eccentricity + baselinenms + nblinks + contrast:baselinenms + eccentricity:baselinenms + contrast:nblinks + eccentricity:nblinks, df)

lm8b = lmer(FT ~ (1 + baselinenms | subject) + contrast + eccentricity + baselinenms, df)
lm9 = lmer(FT ~ (1 | subject) + contrast*eccentricity*baselinenms, df)
lm10 = lmer(FT ~ (1 + baselinenms | subject) + contrast*baselinenms + eccentricity*baselinenms, df)
lm11 = lmer(FT ~ (1 | subject) + contrast*nblinks + eccentricity*nblinks + baselinenms, df)
lm12 = lmer(FT ~ (1 + baselinenms | subject) + contrast*nblinks + eccentricity*nblinks + baselinenms, df)
lm13 = lmer(FT ~ (1 | subject) + contrast + eccentricity + baselinenms + nblinks + baselineretinalslip, df)
lm19 = lmer(FT ~ (1 | subject) + contrast*baselinenms + eccentricity*baselinenms + contrast*nblinks + eccentricity*nblinks, df)
lm20 = lmer(FT ~ (1 + baselinenms | subject) + contrast*baselinenms + eccentricity*baselinenms + contrast*nblinks + eccentricity*nblinks, df)
lm21 = lmer(FT ~ (1 + nblinks | subject) + contrast*baselinenms + eccentricity*baselinenms + contrast*nblinks + eccentricity*nblinks, df)
lm22 = lmer(FT ~ (1 + baselinenms | subject) + contrast*baselinenms + eccentricity*baselinenms + contrast*nblinks + eccentricity*nblinks + contrast*baselineretinalslip + eccentricity*baselineretinalslip, df)

anova(lm0, lm1, lm2, lm3, lm3b, lm4, lm4b)
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm8b, lm9, lm10, lm11, lm19, lm20, lm21, test="Chisq")


# include out/in baseline nms? (doesn't seem to give a clear difference)
lm0 = lmer(FT ~ (1 | subject) + contrast + eccentricity + baselinenms, df)
lm1 = lmer(FT ~ (1 | subject) + contrast + eccentricity + outbaselinenms + inbaselinenms, df)
lm2 = lmer(FT ~ (1 | subject) + contrast*baselinenms + eccentricity*baselinenms, df)
lm3 = lmer(FT ~ (1 | subject) + contrast*inbaselinenms  + eccentricity*inbaselinenms, df)

anova(lm0, lm1, lm2, lm3)

# first regress out FT -> baseline nms
# same result: only eccentricity predicts baseline nms
lmres = lmer(baselinenms ~ (1 + FT | subject) + FT, df)
df2 = df
df2$baselinenms = residuals(lmres)
lm0 = lm(baselinenms ~ 1, df2)
lm1 = lm(baselinenms ~ 1 + contrast, df2)
lm2 = lm(baselinenms ~ 1 + eccentricity, df2)
lm3 = lm(baselinenms ~ 1 + contrast + eccentricity, df2)
lm4 = lm(baselinenms ~ 1 + contrast*eccentricity, df2)
anova(lm0, lm1, lm2, lm3, lm4)


# predict nms / nblinks by contrast/eccentricity
# best lm2: higher eccentricity, more microsaccades. No contrast effect.
# We should report lm3.
lm0 = lmer(baselinenms ~ (1 + FT | subject) + FT, df)
lm1 = lmer(baselinenms ~ (1 + FT | subject) + contrast + FT, df)
lm2 = lmer(baselinenms ~ (1 + FT | subject) + eccentricity + FT, df)
lm3 = lmer(baselinenms ~ (1 + FT | subject) + contrast + eccentricity + FT, df)
lm4 = lmer(baselinenms ~ (1 + FT | subject) + contrast*eccentricity + FT, df)
lm5 = lmer(baselinenms ~ (1 + FT | subject) + contrast*FT + eccentricity*FT, df)
lm6 = lmer(baselinenms ~ (1 + FT | subject) + contrast*eccentricity*FT, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5)

# for nblinks, FT predicts but some don't converge
lm0 = lmer(nblinks ~ (1 + FT | subject) + FT, df)
lm1 = lmer(nblinks ~ (1 + FT | subject) + contrast + FT, df)
lm2 = lmer(nblinks ~ (1 + FT | subject) + eccentricity + FT, df)
lm3 = lmer(nblinks ~ (1 + FT | subject) + contrast + eccentricity + FT, df)
lm4 = lmer(nblinks ~ (1 + FT | subject) + contrast*eccentricity + FT, df)
lm5 = lmer(nblinks ~ (1 + FT | subject) + contrast*eccentricity*FT, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5)

# for retinalslip, contrast predicts (higher contrast = lower slip), FT might also weakly.
lm0 = lmer(baselineretinalslip ~ (1 | subject) + FT, df)
lm1 = lmer(baselineretinalslip ~ (1 | subject) + contrast + FT + baselinenms, df)
lm2 = lmer(baselineretinalslip ~ (1 | subject) + eccentricity + FT, df)
lm3 = lmer(baselineretinalslip ~ (1 | subject) + contrast + eccentricity + FT, df)
lm4 = lmer(baselineretinalslip ~ (1 | subject) + contrast*eccentricity + FT, df)
lm5 = lmer(baselineretinalslip ~ (1 | subject) + contrast*eccentricity*FT, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5)

# predict baseline RATE, only from trials with at least 1 microsaccade.
df2 = dforig[(dforig$baselinenms > 0),]
lm0 = lmer(baselinemsrate ~ (1 + FT | subject) + FT, df2)
lm1 = lmer(baselinemsrate ~ (1 + FT | subject) + contrast + FT, df2)
lm2 = lmer(baselinemsrate ~ (1 + FT | subject) + eccentricity + FT, df2)
lm3 = lmer(baselinemsrate ~ (1 + FT | subject) + contrast + eccentricity + FT, df2)
lm3b = lmer(baselinemsrate ~ (1 + FT | subject) + contrast + eccentricity, df2)
lm3c = lmer(baselinemsrate ~ (1 | subject) + contrast + eccentricity, df2)
lm4 = lmer(baselinemsrate ~ (1 + FT | subject) + contrast*eccentricity + FT, df2)
lm5 = lmer(baselinemsrate ~ (1 + FT | subject) + contrast*FT + eccentricity, df2)
lm6 = lmer(baselinemsrate ~ (1 + FT | subject) + contrast*eccentricity*FT, df2)
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6)
boo3 <- bootstrap(model = lm3, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo3, level = 0.9875, type = "perc")

# predict ms presence as function of contrast / eccentricity.
lm0 = glmer(mspresent ~ (1 + FT | subject) + FT, family = binomial, df)
lm1 = glmer(mspresent ~ (1 + FT | subject) + contrast + FT, family = binomial, df)
lm2 = glmer(mspresent ~ (1 + FT | subject) + eccentricity + FT, family = binomial, df)
lm3 = glmer(mspresent ~ (1 + FT | subject) + contrast + eccentricity + FT, family = binomial, df)
lm4 = glmer(mspresent ~ (1 + FT | subject) + contrast*eccentricity + FT, family = binomial, df)
lm5 = glmer(mspresent ~ (1 + FT | subject) + contrast*FT + eccentricity, family = binomial, df)
lm6 = glmer(mspresent ~ (1 + FT | subject) + contrast*eccentricity*FT, family = binomial, df)
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6)


# select only certain contrast or eccentricity, then re-run model
# kind of unclear what's going on here, or what the appropriate method is.
df1 = df[df$contrast == 1,]
df2 = df[df$contrast == 2,]
df3 = df[df$contrast == 3,]
lm21 = lmer(FT ~ (1 + baselinenms | subject) + eccentricity*baselinenms + eccentricity*nblinks + eccentricity*baselineretinalslip, df1)
lm22 = lmer(FT ~ (1 + baselinenms | subject) + eccentricity*baselinenms + eccentricity*nblinks + eccentricity*baselineretinalslip, df2)
lm23 = lmer(FT ~ (1 + baselinenms | subject) + eccentricity*baselinenms + eccentricity*nblinks + eccentricity*baselineretinalslip, df3)

df1 = df[df$eccentricity == 6,]
df2 = df[df$eccentricity == 4,]
df3 = df[df$eccentricity == 2,]
lm21 = lmer(FT ~ (1 + baselinenms | subject) + contrast*baselinenms + contrast*nblinks + contrast*baselineretinalslip, df1)
lm22 = lmer(FT ~ (1 + baselinenms | subject) + contrast*baselinenms + contrast*nblinks + contrast*baselineretinalslip, df2)
lm23 = lmer(FT ~ (1 + baselinenms | subject) + contrast*baselinenms + contrast*nblinks + contrast*baselineretinalslip, df3)

tab_model(lm21, lm22, lm23)

# above but in trials with no microsaccades
lm21 = lmer(FT ~ (1 | subject) + baselineretinalslip, df1)
lm22 = lmer(FT ~ (1 | subject) + baselineretinalslip, df2)
lm23 = lmer(FT ~ (1 | subject) + baselineretinalslip, df3)
tab_model(lm21, lm22, lm23)


## plot data

p <- ggplot(df, aes(x=baselineretinalslip, y=FT, colour=eccentricity)) + 
  geom_point(size=3) +
  geom_line(aes(y=predict(lm20), group=subject))
print(p)

df2 <- df %>%
  convert_as_factor(eccentricity, contrast)
df2$FT = residuals(lmer(FT ~ (1 | subject) + eccentricity + eccentricity:baselinenms + contrast:nblinks + nblinks + eccentricity:nblinks, df2))
df2$FT = residuals(lmer(FT ~ (1 | subject) + contrast + contrast:baselinenms + eccentricity:nblinks + nblinks + contrast:nblinks, df2))
lmres = lmer(FT ~ (1 + baselinenms | subject) + baselinenms + contrast*nblinks + eccentricity*nblinks, df2)
lmres = lmer(FT ~ (1 | subject) + nblinks, df)
df2$FT = residuals(lmres)

df2 <- df2 %>%
  convert_as_factor(eccentricity, contrast)
p <- ggplot(df2, aes(x=baselineretinalslip, y=FT, colour=contrast)) +
  geom_smooth(method = 'lm', se = F)
p


# look into "immobilization duration"
# confusing but seems like eccentricity is really the predictor.
# smaller eccentricity = last microsaccade has to happen earlier, or else it'll be too strong
# and prevent filling.
# best is lm3.
# when ignoring immobilization durations < 300 ms, best model still lm3
#, but can use lm5 or lm7: no sig contrast or lastmsamplitude.
df = dforig
df[df$immobduration < 0.3 & !is.na(df$immobduration),]$immobduration = NaN # ignore immobdurations less than 300ms (within motor RT)
df[df$immobduration > 0.3 & !is.na(df$immobduration),]$immobduration = NaN # ignore immobdurations longer than 300ms (outside motor RT)

# replace immobduration/amplitude with lastmstime/amplitude
df$immobduration = df$lastmstime
df$immobmsamplitude = df$lastmsamplitude

lm0 = lmer(immobduration ~ (1 | subject), df)
lm1 = lmer(immobduration ~ (1 | subject) + contrast, df)
lm2 = lmer(immobduration ~ (1 | subject) + eccentricity, df)
lm3 = lmer(immobduration ~ (1 | subject) + contrast + eccentricity, df)
#lm3a = lmer(immobduration ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(immobduration ~ (1 | subject) + contrast*eccentricity, df)
lm5 = lmer(immobduration ~ (1 | subject) + contrast + eccentricity + immobmsamplitude, df)
lm5star = lmer(immobduration ~ (1 + contrast | subject) + contrast + eccentricity + immobmsamplitude, df)
lm6 = lmer(immobduration ~ (1 | subject) + contrast*immobmsamplitude + eccentricity*immobmsamplitude, df)
lm7 = lmer(immobduration ~ (1 + immobmsamplitude | subject) + contrast + eccentricity + immobmsamplitude, df)
lm8 = lmer(immobduration ~ (1 + immobmsamplitude | subject) + contrast*immobmsamplitude + eccentricity*immobmsamplitude, df)
lm9 = lmer(immobduration ~ (1 + immobmsamplitude | subject) + contrast*immobmsamplitude*eccentricity, df)

#anova(lm0, lm1, lm2, lm3, lm3a, lm4, lm4b, lm5, lm5b, test = 'Chisq')
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9)

boo3 <- bootstrap(model = lm3, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo3, level = 0.983, type = "perc")
save(boo3, file = "immob_bootstrap.Rdata")
boo3b <- bootstrap(model = lm3b, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo3b, level = 0.983, type = "perc")
boo4 <- bootstrap(model = lm4, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo4, level = 0.9875, type = "perc")
boo5 <- bootstrap(model = lm5, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo5, level = 0.9875, type = "perc")
boo6 <- bootstrap(model = lm6, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo6, level = 0.9875, type = "perc")
boo7 <- bootstrap(model = lm7, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo7, level = 0.9875, type = "perc")
save(boo7, file = "immob_bootstrap_withamp.Rdata")
boo8 <- bootstrap(model = lm8, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo8, level = 0.95, type = "perc")


df2 <- df %>%
  convert_as_factor(subject, contrast, eccentricity)#, scaled_eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))
# average across subjects
df3 = df2 %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(immobduration, type = "mean")
# normalize by subject mean, where pop. mean = still pop. mean
pop_mean = mean(na.omit(dforig$immobduration))
for (s in unique(df$subject)) {
  subj_mean = mean(na.omit(dforig[(dforig$subject == s),]$immobduration))
  df3[(df3$subject == s),]$mean = pop_mean * df3[(df3$subject == s),]$mean / subj_mean
}

plottitle="immobilization time"
bxp <- ggboxplot(
  df3, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('immobilization time (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 3))

# boxplot of mean immob durations, averaged within conditions
df4 = df3 %>%
  group_by(subject, eccentricity) %>%
  get_summary_stats(mean, type = "mean")

bxp <- ggboxplot(
  df3, x = "eccentricity", y = "mean",
  color = "eccentricity", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + stat_summary(fun = mean, geom="point", shape=20, size=5, aes(color=paste(eccentricity)))
bxp = bxp + ylab('immobilization time (s)') + theme(text = element_text(size = 20))
bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) # for eccentricity
#bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 3)) # normalized to pop. mean FT

# summary stats
df2 %>%
  group_by(eccentricity) %>%
  get_summary_stats(immobduration, type = "mean_se")

# summary stats FT
df2 %>%
  group_by(eccentricity) %>%
  get_summary_stats(FT, type = "mean_sd")

# boxplots raw
bxp <- ggboxplot(
  df2, x = "eccentricity", y = "immobduration",
  #color = "contrast", palette = "jco",
  #outlier.shape = NA # hide outliers
)
bxp
# don't look at outliers
ylim1 = boxplot.stats(df2$immobduration)$stats[c(1,5)]
bxp1 = bxp + coord_cartesian(ylim = ylim1*1.05)
bxp1

exp((BIC(lm5b) - BIC(lm5))/2)






#### final model comparisons
# random slopes capture the variation of each subject from the population mean slope.
# if you don't also include the fixed effect, it assumes population mean = 0.
# with maximal design you have: estimate of true population effect (fixed), accounting for random subject variation in slope (random)
# final model is 2b.
# 7: mspresent and blinkpresent instead of baselinenms / nblinks. Good but better fit would account for cross-subject variability in baseline ms rate.
# 8: baselinemsRATE and nblinks instead of mspresent / blinkpresent. Bad because msrate not expected to correlate (1 Hz is 1 Hz...)
# 9: both rate and mspresent.
# Xb: include contrast:eccentricity
# Xstar = no baselineretinalslip random slope.
# models with b: include contrast:eccentricity interaction.
# 22: like 2, but replace baselinenms with baselinemsrate.
# I think 22bstar is new best. no baselinemsrate effect (expected), but ms present effect.
    # the difference from 9bstar is that 9bstar also includes interactions with baselinemsrate. I can argue they are irrelevant.
# current understanding is: it's possible that baselinemsrate and/or mspresent
# could predict FT and interact with ecc. So it's an empirical question which
# shows up.
# nblinks is better than blinkpresent.
lmallmax0 = lmer(FT ~ (nblinks + baselinenms + baselineretinalslip || subject) + contrast*mspresent + baselinenms + contrast*nblinks + contrast*baselineretinalslip + eccentricity*mspresent + eccentricity*nblinks + eccentricity*baselineretinalslip, df) # no baselinenms interactions
lmallmax0b = lmer(FT ~ (nblinks + baselinenms + baselineretinalslip || subject) + contrast*mspresent + baselinenms + contrast*nblinks + contrast*baselineretinalslip + eccentricity*mspresent + eccentricity*nblinks + eccentricity*baselineretinalslip + contrast:eccentricity, df) # include retinal slip, no baselinenms interactions
lmallmax1 = lmer(FT ~ (nblinks + baselinenms + baselineretinalslip || subject) + contrast*mspresent + contrast*baselinenms + contrast*nblinks + eccentricity*mspresent + eccentricity*baselinenms + eccentricity*nblinks, df) # no fixed retinal slip, but include baselinenms interactions
lmallmax1b = lmer(FT ~ (baselineretinalslip + baselinenms + nblinks || subject) + contrast*mspresent + contrast*baselinenms + contrast*nblinks + eccentricity*mspresent + eccentricity*baselinenms + eccentricity*nblinks + contrast:eccentricity, df) # no fixed retinal slip, but include baselinenms interactions
lmallmax2 = lmer(FT ~ (baselineretinalslip + baselinenms + nblinks || subject) + contrast*mspresent + baselinenms + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks, df) # no baselinenms interactions.
#lmallmax3 = lmer(FT ~ (baselinenms + baselineretinalslip + nblinks || subject) + contrast*mspresent + baselinenms + contrast*nblinks + scaled_eccentricity*mspresent + scaled_eccentricity*nblinks, df) # like 2, using scaled_eccentricity instead.
lmallmax2b = lmer(FT ~ (baselineretinalslip + baselinenms + nblinks || subject) + contrast*mspresent + baselinenms + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks + contrast:eccentricity, df)
lmallmax2c = lmer(FT ~ (baselineretinalslip + baselinenms + nblinks || subject) + contrast*mspresent + baselinenms + baselineretinalslip + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks + contrast:eccentricity, df) # like 2b but also include fixed baselineretinalslip (no interactions though)
lmallmaxfull = lmallmax0b # no baselinenms interactions
lmallmax4b = lmer(FT ~ (baselineretinalslip + baselinenms + nblinks || subject) + contrast*baselinenms + contrast*nblinks + eccentricity*baselinenms + eccentricity*nblinks + contrast:eccentricity, df) # no mspresent and no fixed baselineretinalslip
lmallmax5b = lmer(FT ~ (baselineretinalslip + baselinenms + nblinks || subject) + contrast*mspresent + nblinks + baselinenms + contrast*blinkpresent + eccentricity*mspresent + eccentricity*blinkpresent + contrast:eccentricity, df) # also include blinkpresent, in the same manner as mspresent.
lmallmax6b = lmer(FT ~ (baselineretinalslip + nblinks || subject) + contrast*mspresent + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks + contrast:eccentricity, df) # no baselinenms at all.
lmallmax6star = lmer(FT ~ (nblinks || subject) + contrast*mspresent + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks, df)
lmallmax7 = lmer(FT ~ (baselineretinalslip || subject) + contrast*mspresent + contrast*blinkpresent + eccentricity*mspresent + eccentricity*blinkpresent, df)
lmallmax7star = lmer(FT ~ (mspresent + blinkpresent | subject) + contrast*mspresent + contrast*blinkpresent + eccentricity*mspresent + eccentricity*blinkpresent, df)
lmallmax7b = lmer(FT ~ (baselineretinalslip || subject) + contrast*mspresent + contrast*blinkpresent + eccentricity*mspresent + eccentricity*blinkpresent + contrast:eccentricity, df) # no baselinenms at all, and no nblinks: blinkpresent instead.
lmallmax8 = lmer(FT ~ (baselinemsrate + baselineretinalslip + nblinks || subject) + contrast*baselinemsrate + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*nblinks, df) # convert microsaccades to baseline rate.
lmallmax8star = lmer(FT ~ (baselinemsrate + nblinks || subject) + contrast*baselinemsrate + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*nblinks, df) # 8 with no random retinal slip
lmallmax8b = lmer(FT ~ (baselinemsrate + baselineretinalslip + nblinks || subject) + contrast*baselinemsrate + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*nblinks + contrast:eccentricity, df) # convert microsaccades to baseline rate.
lmallmax8c = lmer(FT ~ (contrast*baselinemsrate + eccentricity*baselinemsrate + baselineretinalslip + nblinks || subject) + contrast*baselinemsrate + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*nblinks, df) # and also include random slopes for C/E interactions with baselinerate.
lmallmax8cb = lmer(FT ~ (contrast*baselinemsrate + eccentricity*baselinemsrate + baselineretinalslip + nblinks || subject) + contrast*baselinemsrate + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*nblinks + contrast:eccentricity, df) # and also include random slopes for C/E interactions with baselinerate.
lmallmax9 = lmer(FT ~ (baselinemsrate + nblinks || subject) + contrast*baselinemsrate + contrast*nblinks + contrast*mspresent + eccentricity*baselinemsrate + eccentricity*nblinks + eccentricity*mspresent, df) # baseline ms RATE AND mspresent, no random retinalslip (bc doesn't converge with it)
lmallmax9star = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*baselinemsrate + contrast*nblinks + contrast*mspresent + eccentricity*baselinemsrate + eccentricity*nblinks + eccentricity*mspresent, df) # baseline ms RATE AND mspresent, no random retinalslip (bc doesn't converge with it)
lmallmax9bstar = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*baselinemsrate + contrast*nblinks + contrast*mspresent + eccentricity*baselinemsrate + eccentricity*nblinks + eccentricity*mspresent + contrast:eccentricity, df) # baseline ms RATE AND mspresent, no random retinalslip (bc doesn't converge with it)
lmallmax22 = lmer(FT ~ (baselineretinalslip + baselinemsrate + nblinks || subject) + contrast*mspresent + baselinemsrate + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks, df)
lmallmax22b = lmer(FT ~ (baselineretinalslip + baselinemsrate + nblinks || subject) + contrast*mspresent + baselinemsrate + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks + contrast:eccentricity, df)
lmallmax22star = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*mspresent + baselinemsrate + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks, df) # include mspresent random slope.
lmallmax22bstar = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*mspresent + baselinemsrate + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks + contrast:eccentricity, df) # include mspresent random slope.
lmallmax22full = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*mspresent + contrast*baselinemsrate + contrast*nblinks + contrast*baselineretinalslip + eccentricity*mspresent + eccentricity*baselinemsrate + eccentricity*nblinks + eccentricity*baselineretinalslip + contrast:eccentricity, df)
lmallmax23star = lmer(FT ~ (baselinemsrate + blinkpresent + mspresent | subject) + contrast*mspresent + baselinemsrate + contrast*blinkpresent + eccentricity*mspresent + eccentricity*blinkpresent, df) # include mspresent random slope.

lmallmax2bposthoc = lmallmax1b
lmallmax2cposthoc = lmer(FT ~ (baselinenms + baselineretinalslip + nblinks || subject) + contrast*mspresent + contrast*baselinenms + contrast*nblinks + baselineretinalslip + eccentricity*mspresent + eccentricity*baselinenms + eccentricity*nblinks + contrast:eccentricity, df)
lmallmaxfullposthoc = lmer(FT ~ (baselinenms + nblinks + baselineretinalslip || subject) + contrast*mspresent + contrast*baselinenms + contrast*nblinks + contrast*baselineretinalslip + eccentricity*mspresent + eccentricity*baselinenms + eccentricity*nblinks + contrast:eccentricity + eccentricity*baselineretinalslip, df)
# full posthoc has everything: retinalslip+interactions and baselinenms+interactions.
lmallmax4bposthoc = lmer(FT ~ (baselinenms + baselineretinalslip + nblinks || subject) + contrast*baselinenms + contrast*nblinks + contrast*baselineretinalslip + eccentricity*baselinenms + eccentricity*nblinks + eccentricity*baselineretinalslip + contrast:eccentricity, df) # no mspresent, but include baselineretinalslip+interactions.
lmallmax7bposthoc = lmer(FT ~ (baselineretinalslip || subject) + contrast*mspresent + contrast*blinkpresent + contrast*baselineretinalslip + eccentricity*mspresent + eccentricity*blinkpresent + eccentricity*baselineretinalslip + contrast:eccentricity, df) # no baselinenms at all, and no nblinks: blinkpresent instead. Posthoc also add retinalslip effects.

exp((BIC(lmallmax9bstar) - BIC(lmallmax22bstar))/2)

# make baselinemsrate variable. The rationale is: baselinenms trivially correlates with RT, but baselinsmrate doesn't
# however baselinemsrate should interact with stimulus properties if they do.
df$baselinemsrate = df$baselinenms / df$FT

# bootstrap
mySumm <- function(.) {
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta")))
}

set.seed(1234)
# I only want to bootstrap the fixed effects; those are what I'm interested in.
# For resampling, I resample subsets of subjects. If a subject is included
# in a bootstrap sample, we include all of their observations.
boo2b <- bootstrap(model = lmallmax2b, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
boo0b <- bootstrap(model = lmallmax0b, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo1b <- bootstrap(model = lmallmax1b, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo4b <- bootstrap(model = lmallmax4b, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo7 <- bootstrap(model = lmallmax7, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo7star <- bootstrap(model = lmallmax7star, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo7b <- bootstrap(model = lmallmax7b, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo8 <- bootstrap(model = lmallmax8, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo8star <- bootstrap(model = lmallmax8star, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo9 <- bootstrap(model = lmallmax9, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo9star <- bootstrap(model = lmallmax9star, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo9bstar <- bootstrap(model = lmallmax9bstar, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo22bstar <- bootstrap(model = lmallmax22bstar, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo22bfull <- bootstrap(model = lmallmax22full, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boo2bposthoc <- bootstrap(model = lmallmax2bposthoc, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
boofullposthoc <- bootstrap(model = lmallmaxfullposthoc, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
#boo2_pval <- bootstrap_pvals(model = lmallmax2, type = "case", B = 1000, resample = c(TRUE, FALSE))
# which confidence interval? not normal - that assumes normal distribution.
# basic = like percentile but somewhat corrected for bootstrap bias I think?
# percentile = self explanatory, the 2.5th & 97.5th percentiles from the bootstraps.
confint(boo2b, level = 0.9955, type = "perc")
confint(boo1b, level = 0.9962, type = "perc")
confint(boo0b, level = 0.9964, type = "perc")
confint(boo4b, level = 0.95, type = "perc")
confint(boo7, level = 0.9944, type = "perc")
confint(boo7star, level = 0.9944, type = "perc")
confint(boo7b, level = 0.995, type = "perc")
confint(boo8, level = 0.9944, type = "perc")
confint(boo8star, level = 0.9944, type = "perc")
confint(boo8c, level = 0.995, type = "perc")
confint(boo9, level = 0.9958, type = "perc")
confint(boo9star, level = 0.9964, type = "perc")
confint(boo9bstar, level = 0.9964, type = "perc")
confint(boo22bstar, level = 0.9964, type = "perc")
confint(boo22bfull, level = 0.997, type = "perc")
confint(boo2bposthoc, level = 0.9955, type = "perc")
confint(boofullposthoc, level = 0.9964, type = "perc")
# visualize
plot(boo2b)

# save anything?
save(boo2b, file = "final_bootstrap.Rdata")

# copy pvalues into a list called pvalues, then: (requires library stats)
padjusted = p.adjust(pvalues, method = "fdr")

# - # - # ALTERNATIVE METHOD - initial parameter selection # - # - #
# first model everything with no interactions and only keep sig. effects (liberal p = 0.05, uncorrected)
# bootstrap using extract_parameters: first vc = random intercept, last vc = sigma2 (standard deviation of the residuals)
lm0 = lmer(FT ~ (1 | subject) + baselinemsrate + mspresent + nblinks + blinkpresent + baselineretinalslip + contrast + eccentricity, df)

# contrast coding: I use the default (contr.treatment), which means factors (categorical variables)
# act where 0 = reference, 1 = treatment. Other main effects correspond to when that factor = 0,
# e.g., when mspresent = 0. Other option, contr.sum, sets factors to -1 and 1, so other main effects
# correspond to the MEAN of both factor options - i.e., avg of mspresent and msabsent.
#options(contrasts = c("contr.sum", "contr.poly"))
#options(contrasts = c("contr.treatment", "contr.poly"))

lmallmain = lmer(FT ~ (baselinemsrate + mspresent + nblinks + blinkpresent | subject) + baselinemsrate + mspresent + nblinks + blinkpresent + baselineretinalslip + contrast + eccentricity, df)
#lmallmain = lmer(FT ~ (mspresent + nblinks + baselinemsrate | subject) + baselinemsrate + mspresent + nblinks + blinkpresent + baselineretinalslip + contrast, df)
booallmain <- bootstrap(model = lmallmain, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
ciallmain = confint(booallmain, level = 0.99375, type = "perc")

lmallmainearlylate = lmer(FT ~ (baselinemsrate + earlymspresent + latemspresent| subject) + baselinemsrate + earlymspresent + latemspresent + nblinks + blinkpresent + contrast + eccentricity, df)
lmallmainearlylate2 = lmer(FT ~ (baselinemsrate + latemspresent + nblinks + blinkpresent | subject) + baselinemsrate + earlymspresent + latemspresent + nblinks + blinkpresent + contrast + eccentricity, df) # better.
lmallmainimmob = lmer(FT ~ (baselinemsrate + immobmsamplitude + nblinks | subject) + baselinemsrate + nblinks + blinkpresent + contrast + eccentricity + immobduration + immobmsamplitude, df)

booallmainearlylate <- bootstrap(model = lmallmainearlylate2, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
ciallmainearlylate = confint(booallmainearlylate, level = 0.9944, type = "perc")

# compare to model with non-sig effects removed.
lmtrimmain = lmer(FT ~ (baselinemsrate + mspresent + nblinks + blinkpresent | subject) + contrast + eccentricity + mspresent + nblinks, df)
anova(lmallmain, lmtrimmain)
# (note fit is tiniest bit better without random blinkpresent, but tiny bit worse without random baselinemsrate)
# (I think whatever, just leave them in.)
lmtrimmain2 = lmer(FT ~ (baselinemsrate + mspresent + nblinks | subject) + contrast + eccentricity + mspresent + nblinks, df)
lmtrimmain3 = lmer(FT ~ (mspresent + nblinks | subject) + baselinemsrate + contrast + eccentricity + mspresent + nblinks, df)


# then add interactions with C/E between those, and test with Bonferroni-correction.
# singular fit if include random blinkpresent.
lmallselectedtmp = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*eccentricity + contrast*mspresent + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks, df)
lmallselected = lmer(FT ~ (nblinks + mspresent | subject) + contrast*eccentricity + contrast*mspresent + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks, df)

lmallselected2 = lmer(FT ~ (baselinemsrate + mspresent + nblinks | subject) + contrast*eccentricity + contrast*baselinemsrate + contrast*mspresent + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*mspresent + eccentricity*nblinks, df)
lmallselectedposthoc = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*mspresent + contrast*nblinks + eccentricity*mspresent + eccentricity*nblinks, df) # no C*E interaction
lmallselectedposthoc2 = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast + eccentricity*mspresent + eccentricity*nblinks, df) # no insignificant interactions
lmallselectedposthoc3 = lmer(FT ~ (baselinemsrate + nblinks + mspresent | subject) + contrast*eccentricity + contrast*mspresent + contrast*nblinks + eccentricity*mspresent, df) # no ecc*nblinks interaction, fits worse.

lmallselectedimmob = lmer(FT ~ (baselinemsrate + nblinks | subject) + contrast*eccentricity + contrast*immobduration + contrast*nblinks + eccentricity*immobduration + eccentricity*nblinks, df)

exp((BIC(lmallselected) - BIC(lmallselectedb))/2)

# nothing sig different if splitting microsaccades into early vs late.
# specifically: maybe early microsaccades aren't really effective at all.
lmallselectedearly = lmer(FT ~ (baselinemsrate + nblinks + earlymspresent | subject) + contrast*eccentricity + contrast*earlymspresent + contrast*nblinks + eccentricity*earlymspresent + eccentricity*nblinks + latemspresent, df)
lmallselectedlate = lmer(FT ~ (baselinemsrate + nblinks + latemspresent | subject) + contrast*eccentricity + contrast*latemspresent + contrast*nblinks + eccentricity*latemspresent + eccentricity*nblinks + earlymspresent, df)

booallselected <- bootstrap(model = lmallselected, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(booallselected, level = 0.99, type = "perc")
booallselected2 <- bootstrap(model = lmallselected2, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(booallselected2, level = 0.95, type = "perc")

# compare to model without mspresent - it's worse and doesn't fully account for eccentricity effect.
lmallselectedstar = lmer(FT ~ (baselinemsrate + nblinks | subject) + contrast*eccentricity + contrast*baselinemsrate + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*nblinks, df)
anova(lmallselected, lmallselectedstar)

# save bootstraps
save(booallmain, booallselected, file = "final_mainselect_bootstrap.Rdata")

# predict nms / nblinks by contrast/eccentricity
# best lm2: higher eccentricity, more microsaccades. No contrast effect.
# We should report lm3.
lm0 = lmer(baselinenms ~ (1 + FT | subject) + FT, df)
lm1 = lmer(baselinenms ~ (1 + FT | subject) + contrast + FT, df)
lm2 = lmer(baselinenms ~ (1 + FT | subject) + eccentricity + FT, df)
lm3 = lmer(baselinenms ~ (1 + FT | subject) + contrast + eccentricity + FT, df)
lm3rate = lmer(baselinemsrate ~ (1 + FT | subject) + contrast + eccentricity + FT, df)
lm4 = lmer(baselinenms ~ (1 + FT | subject) + contrast*eccentricity + FT, df)
lm5 = lmer(baselinenms ~ (1 + FT | subject) + contrast*FT + eccentricity*FT, df)
lm6 = lmer(baselinenms ~ (1 + FT | subject) + contrast*eccentricity*FT, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5)
boo3 <- bootstrap(model = lm3, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo3, level = 0.9875, type = "perc")
boo3rate <- bootstrap(model = lm3rate, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo3rate, level = 0.9875, type = "perc")
boo5 <- bootstrap(model = lm5, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo5, level = 0.9875, type = "perc")

# save bootstraps
save(boo3, boo3rate, file = "final_msgeneration_bootstraps.Rdata")




## now only trials without microsaccades or blinks
# eccentricity-retinalslip interaction. No contrast effect, but it's because including the contrast*slip interaction messes it up.
# model without that interaction is much better fit and has a strong contrast effect.
# that kind of fix doesn't work for the main model above: the only insignificant main effect, baselineretinalslip, does not become
#  significant if you remove its interactions. So the full model is still valid.
df2 = df[(df$baselinenms == 0) & (df$nblinks == 0),]
#df2 = df[(df$latemspresent == 0) & (df$nblinks == 0) & (df$earlymspresent == 1),]
lmallmax2 = lmer(FT ~ (baselineretinalslip || subject) + contrast*baselineretinalslip + eccentricity*baselineretinalslip, df2)
lmallmax3 = lmer(FT ~ (baselineretinalslip || subject) + contrast*baselineretinalslip + eccentricity, df2)
lmallmax4 = lmer(FT ~ (baselineretinalslip || subject) + contrast + eccentricity*baselineretinalslip, df2)
lmallmax5 = lmer(FT ~ (baselineretinalslip || subject) + contrast + baselineretinalslip + eccentricity, df2)
lmallmax6 = lmer(FT ~ (baselineretinalslip | subject) + contrast + baselineretinalslip + eccentricity + earlymspresent, df2)
anova(lmallmax2, lmallmax3, lmallmax4, lmallmax5)

lm5 = lmer(FT ~ (1 | subject) + contrast + eccentricity, df2)

boo4 <- bootstrap(model = lmallmax4, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo4, level = 0.95, type = "perc")

boo5 <- bootstrap(model = lmallmax5, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo5, level = 0.9875, type = "perc")
save(boo5, boo5, file = "final_noms_bootstraps.Rdata")


# for only trials with at least one blink, no microsaccades (only 84 trials)
df2 = df[(df$baselinenms == 0) & (df$nblinks > 0),]
lmallmax2 = lmer(FT ~ (1 | subject) + contrast*nblinks + eccentricity*nblinks + baselineretinalslip, df2)
lmallmax3 = lmer(FT ~ (1 | subject) + contrast*nblinks + eccentricity + baselineretinalslip, df2)
lmallmax4 = lmer(FT ~ (1 | subject) + contrast + eccentricity*nblinks + baselineretinalslip, df2)
lmallmax5 = lmer(FT ~ (1 | subject) + contrast + baselineretinalslip + eccentricity + nblinks, df2)
anova(lmallmax2, lmallmax3, lmallmax4, lmallmax5)

boo5 <- bootstrap(model = lmallmax5, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo5, level = 0.95, type = "perc")

# for only trials with at least one microsaccade
# doesn't look like we can get an interaction bw microsaccades and eccentricity.
# but I think it's a matter of the study design...
df2 = df[(df$baselinenms > 0) & (df$nblinks >= 0),]
lmallmax2 = lmer(FT ~ (baselinenms + nblinks | subject) + contrast*baselinenms + contrast*nblinks + eccentricity*baselinenms + eccentricity*nblinks, df2)
lmallmax3 = lmer(FT ~ (baselinemsrate + nblinks | subject) + contrast*baselinemsrate + contrast*nblinks + eccentricity*baselinemsrate + eccentricity*nblinks, df2)

boo2 <- bootstrap(model = lmallmax2, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
confint(boo2, level = 0.95, type = "perc")
boo3 <- bootstrap(model = lmallmax3, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
confint(boo3, level = 0.95, type = "perc")

# take a subset of onlyms trials, to match the trial # for noms
# need to equally sample each subject:
dfms = dforig[(dforig$baselinenms > 0) & (dforig$nblinks == 0),]
ntrials = 769
ntrialspersubject = round(ntrials / length(unique(dforig$subject)))
df = data.frame()
for (s in unique(dforig$subject)) {
  subj_df = dfms[(dfms$subject == s),]
  if (nrow(subj_df) >= ntrialspersubject) {
    df = rbind(df, subj_df[sample(nrow(subj_df), ntrialspersubject),])
  } else {
    df = rbind(df, subj_df)
  }
}


# # plotting # #
df2 = df
# convert for plotting
df2 <- df2 %>%
  convert_as_factor(subject, contrast, eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))
# average within subjects
df3 = df2 %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(FT, type = "mean")
# normalize by subject mean, where pop. mean = still pop. mean
pop_mean = mean(na.omit(dforig$FT))
for (s in unique(df$subject)) {
  subj_mean = mean(na.omit(dforig[(dforig$subject == s),]$FT))
  df3[(df3$subject == s),]$mean = pop_mean * df3[(df3$subject == s),]$mean / subj_mean
}

# boxplot of mean FTs
#plottitle = "all trials"
#plottitle="at least one microsaccade"
#plottitle="at least one microsaccade, subset"
#plottitle ="at least one blink"
plottitle="no microsaccades or blinks"
bxp <- ggboxplot(
  df3, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15))

# boxplot of mean FTs, averaged within conditions
df4 = df3 %>%
  group_by(subject, contrast) %>%
  get_summary_stats(mean, type = "mean")

bxp <- ggboxplot(
  df4, x = "contrast", y = "mean",
  color = "contrast", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
bxp = bxp + stat_summary(fun = mean, geom="point", shape=20, size=5, aes(color=paste(contrast)))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15)) # normalized to pop. mean FT






### - ##### - fully all trials - ##### - ###
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
df <- read.csv('./full_behavior.csv',sep = ",",stringsAsFactors = FALSE)
#df$FT = df$RT
# subset of subjects?
df <- df %>%
  convert_as_factor(subject, trialtype, block, pressed)
dforig = df

# select any rows?
df = dforig[(dforig$trialtype == 3),]

# predict FT
# (no E effect in replay trials, but still C effect)
# (both effects in sharp trials)

lm0 = lmer(FT ~ (1 | subject), df)
lm1 = lmer(FT ~ (1 | subject) + contrast, df)
lm2 = lmer(FT ~ (1 | subject) + eccentricity, df)
lm3 = lmer(FT ~ (1 | subject) + contrast + eccentricity, df)
lm4 = lmer(FT ~ (1 | subject) + contrast*eccentricity, df)
lm5 = lmer(FT ~ (1 | subject) + contrast + eccentricity + trialtype, df)

lm10 = lmer(FT ~ (1 | subject) + contrast + eccentricity + medianfixdist, df)
lm11 = lmer(FT ~ (1 | subject) + contrast*medianfixdist + eccentricity*medianfixdist, df)

anova(lm0, lm1, lm2, lm3, lm4)

# plot across trial types
df2 = df
# convert for plotting
df2 <- df2 %>%
  convert_as_factor(subject, contrast, eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))
# average within subjects, separated by stim condition
df3 = df2 %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(FT, type = "median")
# normalize by subject mean, where pop. mean = still pop. mean
# doesn't really make sense if we are trying to plot the pop. mean within trialtypes.
pop_mean = mean(na.omit(dforig$FT))
for (s in unique(df$subject)) {
  subj_mean = mean(na.omit(dforig[(dforig$subject == s),]$FT))
  df3[(df3$subject == s),]$mean = pop_mean * df3[(df3$subject == s),]$mean / subj_mean
}
# or average within subjects and across stim conditions
df3 = df2 %>%
  group_by(subject, trialtype) %>%
  get_summary_stats(FT, type = "median")

# t tests or wilcoxon signed rank tests by trialtype
t.test(df3[df3$trialtype == 1,]$mean, df3[df3$trialtype == 3,]$mean, paired = TRUE)
wilcox.test(df3[df3$trialtype == 1,]$median, df3[df3$trialtype == 3,]$median, paired = TRUE)
wilcox.test(df3[df3$trialtype == 2,]$median, mu = 6) # test if replay FT is different from something

# boxplot of mean FTs by C&E (averaged across trialtypes or just 1 trialtype)
plottitle ="all data"
#plottitle = "main trials"
#plottitle = "replay trials"
#plottitle = "sharp trials"
bxp <- ggboxplot(
  df3, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15))

# boxplot of mean FTs, averaged within conditions and across trialtypes
bxp <- ggboxplot(
  df3, x = "contrast", y = "mean",
  color = "contrast", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15)) # normalized to pop. mean FT

# boxplot of mean FTs, averaged within trialtypes
bxp <- ggboxplot(
  df3, x = "trialtype", y = "median",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
bxp = bxp + scale_x_discrete(labels = c("main", "replay", "sharp")) # for contrast
ggpar(bxp, ylim = c(0, 15))







## NEW TESTING NOV 2023

# load packages
library(lme4)
library(sjPlot)
library(tidyverse)
library(ggpubr)
library(rstatix)
#library(mediation)
library(ggplot2)
library(lmeresampler)
library(performance)
library(qqplotr)
library(ggeffects)

## ---- load individual trials ---- ##
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
#setwd("/export04/data/mlevin/UI_eyelink/bhv_data")
df <- read.csv('./trial_data_withblinks.csv',sep = ",",stringsAsFactors = FALSE)
df$FT = df$RT
df$baselinemsrate = df$baselinenms / df$FT
df$blinkrate = df$nblinks / df$FT
df <- df %>%
  convert_as_factor(subject, mspresent, earlymspresent, latemspresent, blinkpresent, session)
dforig = df
# add scaled_eccentricity
A = 0.063
l0 = 36.54 * A
df$scaled_eccentricity = (log(df$eccentricity) - l0) / A
df$scaled_eccentricity = 1 / (A * df$eccentricity)

df[(df$baselinenms == 0),]$baselinemeanmsampl = 0
df$baselinemsampltotal = df$baselinemeanmsampl * df$baselinenms

# rescale and zero-mean all continuous predictors (but not FT)
numcols <- sapply(df, is.numeric)
dfs <- df
dfs[,numcols] <- scale(dfs[,numcols])
#dfs$FT = dforig$FT / 10
dfnocenter <- df
dfnocenter[,numcols] <- scale(dfnocenter[,numcols], center = FALSE, scale = apply(dfnocenter[,numcols], 2, sd, na.rm = TRUE))
dfs$FT = dfnocenter$FT
dfs$immobduration = dfnocenter$immobduration
#dfs = dfnocenter
logdfs = dfs
logdfs$FT = log(dfs$FT)

# log of amplitude too
logdfsampl = logdfs
logdfsampl$baselinemeanmsampl = log(dfnocenter$baselinemeanmsampl)

# predict baseline RATE, only from trials with at least 1 microsaccade.
# best is lm2 (higher ecc -> higher rate), but use lm3 to show contrast is not significant.
# actually with random C/E effects, no consistent population diff. (lm6)
# but later trial: more microsaccades?
# don't predict number of microsaccades, because it would trivially correlate with C & E (longer trials)
# use LOG baselinemsrate to aid model fit. (lm6log)
dfswithms = dfs
dfswithms$baselinemsrate = dfnocenter$baselinemsrate
dfswithms = dfswithms[(dforig$baselinenms > 0),]
lm0 = lmer(baselinemsrate ~ (1 | subject), dfswithms)
lm1 = lmer(baselinemsrate ~ (1 | subject) + contrast, dfswithms)
lm2 = lmer(baselinemsrate ~ (1 | subject) + eccentricity, dfswithms)
lm3 = lmer(baselinemsrate ~ (1 | subject) + contrast + eccentricity, dfswithms)
lm3b = lmer(baselinemsrate ~ (1 + contrast + eccentricity | subject) + contrast + eccentricity, dfswithms)
lm4 = lmer(baselinemsrate ~ (1 | subject) + contrast*eccentricity, dfswithms)
lm5 = lmer(baselinemsrate ~ (1 + trial | subject) + contrast + eccentricity + trial, dfswithms, control=lmerControl("bobyqa"))
lm6 = lmer(baselinemsrate ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, dfswithms, control=lmerControl("bobyqa"))
logdfswithms = dfswithms
logdfswithms$baselinemsrate = log(logdfswithms$baselinemsrate)
lm6log = lmer(baselinemsrate ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, logdfswithms, control=lmerControl("bobyqa"))
anova(lm0, lm1, lm2, lm3, lm4)
boo3 <- bootstrap(model = lm3, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boo3, level = 0.9833, type = "perc")
boo6 <- bootstrap(model = lm6log, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo6, level = 0.9875, type = "perc")
save(boo6, file = "baselinerate_bootstrap.Rdata")

# don't use glm, too complex here (hard to fit with trial parameter) and we can do log(baselinemsrate) I think.
glm3 = glmer(baselinemsrate ~ (1 | subject) + contrast + eccentricity, dfswithms, family=inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm3gamma = glmer(baselinemsrate ~ (1 | subject) + contrast + eccentricity, dfswithms, family=Gamma(link="identity"), control=glmerControl("bobyqa"))

glm6gamma = glmer(baselinemsrate ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, dfswithms, family=Gamma(link="log"), control=glmerControl("bobyqa"))
glm6invgauss = glmer(baselinemsrate ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, dfswithms, family=inverse.gaussian(link="log"), control=glmerControl("bobyqa"))


# look into "immobilization duration"
# best is lm2 (just eccentricity)
# smaller eccentricity = last microsaccade has to happen earlier, or else it'll be too strong
# and prevent filling.
# but use lm5 for paper: no sig contrast or lastmsamplitude.

# I think I can use log actually, immobduration doesn't have as obvious an interpretation?
# conclusions are the same as with glm, but better convergence with log.
# Though: can't directly interpret the parameter estimate (i.e., in milliseconds)
df = dfs
#dfs[dforig$immobduration < 0.3 & !is.na(dforig$immobduration),]$immobduration = NaN # ignore immobdurations less than 300ms (within motor RT)

# replace immobduration/amplitude with lastmstime/amplitude
df$immobduration = df$lastmstime
df$immobmsamplitude = df$lastmsamplitude

lm0 = lmer(immobduration ~ (1 | subject), df, control=lmerControl("bobyqa"))
lm1 = lmer(immobduration ~ (contrast | subject) + contrast, df, control=lmerControl("bobyqa"))
lm2 = lmer(immobduration ~ (eccentricity | subject) + eccentricity, df, control=lmerControl("bobyqa"))
lm3 = lmer(immobduration ~ (contrast + eccentricity | subject) + contrast + eccentricity, df, control=lmerControl("bobyqa"))
#lm3a = lmer(immobduration ~ (contrast + scaled_eccentricity | subject) + contrast + scaled_eccentricity, df, control=lmerControl("bobyqa"))
lm4 = lmer(immobduration ~ (contrast + eccentricity + immobmsamplitude | subject) + contrast + eccentricity + immobmsamplitude, df, control=lmerControl("bobyqa"))
lm5 = lmer(immobduration ~ (contrast + eccentricity + immobmsamplitude + trial | subject) + contrast + eccentricity + immobmsamplitude + trial, df, control=lmerControl("bobyqa"))
anova(lm0, lm1, lm2, lm3, lm4, lm5)

boo5 <- bootstrap(model = lm5, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo5, level = 0.9875, type = "perc")
save(boo5, file = "immob_bootstrap.Rdata")

logdf = df
logdf$immobduration = log(df$immobduration)
lm3log = lmer(immobduration ~ (contrast + eccentricity + immobmsamplitude | subject) + contrast + eccentricity + immobmsamplitude, logdf, control=lmerControl("bobyqa"))
lm3logstar = lmer(immobduration ~ (contrast + eccentricity | subject) + contrast + eccentricity, logdf, control=lmerControl("bobyqa"))
lm5log = lmer(immobduration ~ (contrast + eccentricity + immobmsamplitude + trial | subject) + contrast + eccentricity + immobmsamplitude + trial, logdf, control=lmerControl("bobyqa"))
boo5 <- bootstrap(model = lm5log, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo5, level = 0.9875, type = "perc")
save(boo5, file = "immob_bootstrap.Rdata")

glm0 = glmer(immobduration ~ (1 | subject) + immobmsamplitude, df, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm1 = glmer(immobduration ~ (contrast | subject) + contrast + immobmsamplitude, df, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm2 = glmer(immobduration ~ (eccentricity | subject) + eccentricity + immobmsamplitude, df, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm3 = glmer(immobduration ~ (contrast + eccentricity | subject) + contrast + eccentricity + immobmsamplitude, df, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm3star = glmer(immobduration ~ (contrast + eccentricity | subject) + contrast + eccentricity, df, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm4 = glmer(immobduration ~ (contrast + eccentricity | subject) + contrast + eccentricity + trial + immobmsamplitude, df, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm5 = glmer(immobduration ~ (contrast + eccentricity + immobmsamplitude | subject) + contrast + eccentricity + immobmsamplitude, df, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glm5gamma = glmer(immobduration ~ (contrast + eccentricity + immobmsamplitude | subject) + contrast + eccentricity + immobmsamplitude, df, family = Gamma(link="identity"), control=glmerControl("bobyqa"))
boo4 <- bootstrap(model = glm4, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
confint(boo4, level = 0.95, type = "perc") # significant eccentricity and trial effects.

# get summary stats
dfsummary <- dforig %>%
  convert_as_factor(eccentricity) %>%
  group_by(subject, eccentricity) %>%
  get_summary_stats(immobduration, type = "mean") %>%
  group_by(eccentricity) %>%
  get_summary_stats(mean, type = "mean_se")

# get subject averages across all stim conditions
df3 %>%
  group_by(subject) %>%
  get_summary_stats(mean, type = "mean")


## hypothesis driven model (microsaccade hypothesis)
glmhypothesisall = glmer(FT ~ (baselinemsrate*contrast + baselinemsrate*eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glmhypothesisall2 = glmer(FT ~ (baselinenms*contrast + baselinenms*eccentricity | subject) + baselinenms*contrast + baselinenms*eccentricity, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glmhypothesis = glmer(FT ~ (baselinemsrate + contrast + eccentricity | subject) + contrast*baselinemsrate + eccentricity*baselinemsrate, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glmhypothesistrim1 = glmer(FT ~ (baselinemsrate + contrast + eccentricity | subject) + contrast + eccentricity*baselinemsrate, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glmhypothesis2 = glmer(FT ~ (baselinemsrate + contrast + eccentricity + blinkrate| subject) + contrast*baselinemsrate + eccentricity*baselinemsrate + contrast*blinkrate + eccentricity*blinkrate, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glmhypothesis2trim1 = glmer(FT ~ (baselinemsrate + contrast + eccentricity + blinkrate| subject) + eccentricity*baselinemsrate + contrast*blinkrate + eccentricity*blinkrate, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glmhypothesis2trim2 = glmer(FT ~ (baselinemsrate + contrast + eccentricity + blinkrate| subject) + eccentricity*baselinemsrate + eccentricity*blinkrate, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
# hypothesis driven model with log FT
lmloghypothesisall = lmer(FT ~ (baselinemsrate*contrast + baselinemsrate*eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity, logdfs, control=lmerControl("bobyqa"))
lmloghypothesisall2 = lmer(FT ~ (baselinenms*contrast + baselinenms*eccentricity | subject) + baselinenms*contrast + baselinenms*eccentricity, logdfs, control=lmerControl("bobyqa"))
lmloghypothesisall2trim = lmer(FT ~ (baselinenms*contrast + eccentricity | subject) + baselinenms*contrast + baselinenms*eccentricity, logdfs, control=lmerControl("bobyqa"))
lmloghypothesis2 = lmer(FT ~ (baselinenms + contrast + eccentricity | subject) + baselinenms*contrast + baselinenms*eccentricity, logdfs, control=lmerControl("bobyqa"))
lmloghypothesisalltrim = lmer(FT ~ (baselinemsrate*contrast + eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity, logdfs, control=lmerControl("bobyqa"))
lmloghypothesis = lmer(FT ~ (baselinemsrate + contrast + eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity, logdfs, control=lmerControl("bobyqa"))
boolmloghyp <- bootstrap(model = lmloghypothesis, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
confint(boolmloghyp, level = 0.95, type = "perc") # liberal, take any effects with p < 0.05 uncorrected.
lmloghyptest = lmer(FT ~ (baselinemsrate | subject) + contrast*baselinemsrate + eccentricity*baselinemsrate, logdfs, control=lmerControl("bobyqa"))
lmloghyptesttrim = lmer(FT ~ (1 | subject) + contrast*baselinemsrate + eccentricity*baselinemsrate, logdfs, control=lmerControl("bobyqa"))

loglmhypothesis = lmer(FT ~ (mspresent + nblinks + baselineretinalslip | subject) + contrast*mspresent + eccentricity*mspresent + nblinks + baselineretinalslip + contrast*trial, logdfs, control=lmerControl("bobyqa"))
loglmhypothesis2 = lmer(FT ~ (mspresent + nblinks + baselineretinalslip +contrast + eccentricity | subject) + contrast*mspresent + eccentricity*mspresent + nblinks + baselineretinalslip + trial, logdfs, control=lmerControl("bobyqa"))

lmhypothesisall = lmer(FT ~ (baselinemsrate*contrast + baselinemsrate*eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity, dfs, control=lmerControl("bobyqa"))
lmhypothesis = lmer(FT ~ (baselinemsrate + contrast + eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity, dfs, control=lmerControl("bobyqa"))

lmloghypothesisnew = lmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip | subject) + contrast*baselinemsrate + eccentricity*baselinemsrate + contrast*nblinks + eccentricity*nblinks + contrast*baselineretinalslip + eccentricity*baselineretinalslip, logdfs, control=lmerControl("bobyqa"))
glmhypothesisnew = glmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip | subject) + contrast*baselinemsrate + eccentricity*baselinemsrate + contrast*nblinks + eccentricity*nblinks + contrast*baselineretinalslip + eccentricity*baselineretinalslip, dfs, family = inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
## 3 step model creation. For GLMER I need to: center & rescale all continuous variables
# except for FT, scale it but do NOT center
# Step 1: model with no eye movements.
lmnoeyes = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial + eccentricity*trial, dfs, control=lmerControl(optimizer="bobyqa"))
boonoeyes <- bootstrap(model = lmnoeyes, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
cinoeyes = confint(boonoeyes, level = 0.9917, type = "perc") # liberal, take any effects with p < 0.05 uncorrected.
lmnoeyestrim1 = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*trial + eccentricity*trial, dfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(lmnoeyes) - BIC(lmnoeyestrim1))/2) # BF = 872.55, so evidence for removing contrast*eccentricity
lmnoeyestrim2 = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial, dfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(lmnoeyes) - BIC(lmnoeyestrim2))/2) # BF = 3537.137, so evidence for removing eccentricity*trial
lmnoeyestrimfinal = lmer(FT ~ (trial + contrast + eccentricity | subject) + eccentricity + contrast*trial, dfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(lmnoeyes) - BIC(lmnoeyestrimfinal))/2) # BF = 3062407, so evidence for removing both insignificant parameters.

glmnoeyes = glmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial + eccentricity*trial, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
boonoeyes <- bootstrap(model = glmnoeyes, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
# check BICs here, but I don't think I will use them. Because some interactions might give better fits for now, until
# we later include the eye movement parameters that make everything more valid.
# So instead, use confidence intervals to determine. (even liberal: 95% CI)
glmnoeyestrim1 = glmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*trial + eccentricity*trial, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmnoeyes) - BIC(glmnoeyestrim1))/2) # BF = 0.0009, so evidence for keeping contrast*eccentricity???
glmnoeyestrim2 = glmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmnoeyes) - BIC(glmnoeyestrim2))/2) # BF = 2.08, so weak evidence for removing eccentricity*trial
glmnoeyestrimfinal = glmer(FT ~ (trial + contrast + eccentricity | subject) + eccentricity + contrast*trial, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmnoeyes) - BIC(glmnoeyestrimfinal))/2) # BF = 0.001, so more evidence for keeping both

loglmnoeyes = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial + eccentricity*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmnoeyestrim1 = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*trial + eccentricity*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmnoeyes) - BIC(loglmnoeyestrim1))/2) # BF = 5590, so evidence for removing contrast*eccentricity
loglmnoeyestrim2 = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmnoeyes) - BIC(loglmnoeyestrim2))/2) # BF = 5565, so evidence for removing eccentricity*trial
loglmnoeyestrimfinal = lmer(FT ~ (trial + contrast + eccentricity | subject) + eccentricity + contrast*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmnoeyes) - BIC(loglmnoeyestrimfinal))/2) # BF = 305582, so evidence for removing both insignificant parameters.


# Step 2: use selected parameters from 1, and add eye movement main effects.
lmmain = lmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent + blinkpresent + baselineretinalslip + trial*contrast + eccentricity, dfs, control=lmerControl(optimizer="bobyqa"))
boomain <- bootstrap(model = lmmain, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
cimain = confint(boomain, level = 0.9929, type = "perc")
lmmaintrim1 = lmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent + blinkpresent + trial*contrast + eccentricity, dfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(lmmain) - BIC(lmmaintrim1))/2) # BF = 166.32, so evidence for removing the fixed baselineretinalslip effect
lmmaintrim2 = lmer(FT ~ (mspresent + blinkpresent + trial + contrast + eccentricity | subject) + mspresent + blinkpresent + baselineretinalslip + trial*contrast + eccentricity, dfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(lmmain) - BIC(lmmaintrim2))/2) # BF = 35.25, so evidence for removing the random baselineretinalslip effect
lmmaintrimfinal = lmer(FT ~ (mspresent + blinkpresent + trial + contrast + eccentricity | subject) + mspresent + blinkpresent + trial*contrast + eccentricity, dfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(lmmain) - BIC(lmmaintrimfinal))/2) # BF = 0.000253, so evidence for keeping at least 1 of the removed effects.
# lmmaintrim1 has lowest BIC, so we will use that one: keep random baselineretinalslip but not fixed.

glmmain = glmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent + blinkpresent + baselineretinalslip + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
# glm doesn't converge, so drop 1 random effect at a time. Dropping blinkpresent or baselineretinalslip don't converge either
# lowest BIC is glmmaina: drop random trial effect.
glmmaina = glmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + contrast + eccentricity | subject) + mspresent + blinkpresent + baselineretinalslip + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
glmmainb = glmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + eccentricity | subject) + mspresent + blinkpresent + baselineretinalslip + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
glmmainc = glmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + contrast | subject) + mspresent + blinkpresent + baselineretinalslip + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
glmmainatrim1 = glmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + contrast + eccentricity | subject) + mspresent + blinkpresent + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmmaina) - BIC(glmmainatrim1))/2) # BF = 47.51, so evidence for removing the fixed baselineretinalslip effect
glmmainatrim2 = glmer(FT ~ (mspresent + blinkpresent + contrast + eccentricity | subject) + mspresent + blinkpresent + baselineretinalslip + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmmaina) - BIC(glmmainatrim2))/2) # BF ~0 , so evidence for keeping the random baselineretinalslip effect
# glmmaintrim1 has lowest BIC, so we will use that one: keep random baselineretinalslip but not fixed.
# instead use ms rate / n blinks
glmmain2 = glmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + baselinemsrate + nblinks + baselineretinalslip + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
glmmain2trim1 = glmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + baselinemsrate + nblinks + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmmain2) - BIC(glmmain2trim1))/2) # BF = 77.59, so evidence for removing the fixed baselineretinalslip effect
glmmain2trim2 = glmer(FT ~ (baselinemsrate + nblinks + trial + contrast + eccentricity | subject) + baselinemsrate + nblinks + baselineretinalslip + trial*contrast + eccentricity, dfs, family = inverse.gaussian(link = "identity"), control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmmain2) - BIC(glmmain2trim2))/2) # BF ~0, so evidence for keeping the random baselineretinalslip effect

# log: use loglmmain2 (include mean amplitude too)
loglmmain = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity  + baselineretinalslip | subject) + mspresent + nblinks + baselineretinalslip + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmmain2 = lmer(FT ~ (mspresent + baselinemeanmsampl + nblinks + trial + contrast + eccentricity  + baselineretinalslip | subject) + mspresent + baselinemeanmsampl + nblinks + baselineretinalslip + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmmain3 = lmer(FT ~ (mspresent + baselinemeanmsampl + baselinenms + nblinks + trial + contrast + eccentricity  + baselineretinalslip | subject) + mspresent + baselinemeanmsampl + baselinenms + nblinks + baselineretinalslip + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmmaintrim1 = lmer(FT ~ (mspresent + baselinemeanmsampl + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent + baselinemeanmsampl + nblinks + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmmain2) - BIC(loglmmaintrim1))/2) # BF = 691, so evidence for removing the fixed baselineretinalslip effect
loglmmaintrim2 = lmer(FT ~ (mspresent + baselinemeanmsampl + nblinks + trial + contrast + eccentricity | subject) + mspresent + baselinemeanmsampl + nblinks + baselineretinalslip + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmmain2) - BIC(loglmmaintrim2))/2) # BF ~0, so no evidence for removing the random baselineretinalslip effect
loglmmaintrimfinal = lmer(FT ~ (mspresent + baselinemeanmsampl + nblinks + trial + contrast + eccentricity | subject) + mspresent + baselinemeanmsampl + nblinks + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmmain2) - BIC(loglmmaintrimfinal))/2) # BF ~0, so evidence for keeping at least 1 of the removed effects.
# lmmaintrim1 has lowest BIC, so we will use that one: keep random baselineretinalslip but not fixed.
# Actually: because no fixed, no reason to keep a random effect. So get rid of retinal slip entirely.


# Step 3: use selected parameters from 2, and add interactions bw C/E and eye movements. Bonferroni-corrected interaction selection.
# full model did not converge, so I dropped one random effect at a time and took the one with lowest BIC: lmselectedb (no blinkpresent)
lmselected = lmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + blinkpresent*contrast + blinkpresent*eccentricity + trial*contrast, dfs, control=lmerControl(optimizer="bobyqa"))
lmselecteda = lmer(FT ~ (blinkpresent + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + blinkpresent*contrast + blinkpresent*eccentricity  + trial*contrast, dfs, control=lmerControl(optimizer="bobyqa"))
lmselectedb = lmer(FT ~ (mspresent + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + blinkpresent*contrast + blinkpresent*eccentricity + trial*contrast, dfs, control=lmerControl(optimizer="bobyqa"))
lmselectedc = lmer(FT ~ (mspresent + blinkpresent + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + blinkpresent*contrast + blinkpresent*eccentricity + trial*contrast, dfs, control=lmerControl(optimizer="bobyqa"))
lmselectedd = lmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + blinkpresent*contrast + blinkpresent*eccentricity + trial*contrast, dfs, control=lmerControl(optimizer="bobyqa"))
lmselectede = lmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + blinkpresent*contrast + blinkpresent*eccentricity  + trial*contrast, dfs, control=lmerControl(optimizer="bobyqa"))
lmselectedf = lmer(FT ~ (mspresent + blinkpresent + baselineretinalslip + trial + contrast | subject) + mspresent*contrast + mspresent*eccentricity + blinkpresent*contrast + blinkpresent*eccentricity + trial*contrast, dfs, control=lmerControl(optimizer="bobyqa"))
BIC(lmselected, lmselecteda, lmselectedb, lmselectedc, lmselectedd, lmselectede, lmselectedf)
booselected <- bootstrap(model = lmselectedb, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
ciselected = confint(booselected, level = 0.995, type = "perc")

glmselected = glmer(FT ~ (mspresent + nblinks + trial + baselineretinalslip + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, dfs, family=inverse.gaussian(link="identity"),control=glmerControl(optimizer="bobyqa"))
booselected <- bootstrap(model = glmselected, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
ciselected = confint(booselected, level = 0.95, type = "perc")
# instead use ms rate / n blinks
glmselected2 = glmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, dfs, family=inverse.gaussian(link="identity"),control=glmerControl(optimizer="bobyqa"))
glmselected2trim1 = glmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + baselinemsrate*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, dfs, family=inverse.gaussian(link="identity"),control=glmerControl(optimizer="bobyqa"))
glmselected2trim2 = glmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + baselinemsrate*contrast + baselinemsrate*eccentricity + nblinks*eccentricity + trial*contrast, dfs, family=inverse.gaussian(link="identity"),control=glmerControl(optimizer="bobyqa"))
glmselected2trimfinal = glmer(FT ~ (baselinemsrate + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + baselinemsrate*eccentricity + nblinks*eccentricity + trial*contrast, dfs, family=inverse.gaussian(link="identity"),control=glmerControl(optimizer="bobyqa"))
exp((BIC(glmselected2) - BIC(glmselected2trimfinal))/2) # BF = 120
# trimfinal has lowest BIC and converges: evidence for no baselinerate x contrast OR blink x contrast. Double check with permutations
booselected <- bootstrap(model = glmselected2, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
ciselected = confint(booselected, level = 0.95, type = "perc")
glmselected3 = glmer(FT ~ (baselinenms + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + baselinenms*contrast + baselinenms*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, dfs, family=inverse.gaussian(link="identity"),control=glmerControl(optimizer="bobyqa"))

# use loglmselected2: include baselinemeanmsampl.
loglmselected = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmselected2 = lmer(FT ~ (mspresent + baselinemeanmsampl + nblinks + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + baselinemeanmsampl*contrast + baselinemeanmsampl*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, logdfs, control=lmerControl(optimizer="bobyqa"))

loglmselected2trim1 = lmer(FT ~ (mspresent + baselinemeanmsampl + nblinks + trial + contrast + eccentricity | subject) + mspresent*eccentricity + baselinemeanmsampl*contrast + baselinemeanmsampl*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmselected2) - BIC(loglmselected2trim1))/2) # BF = 1661, so no contrast*mspresent
loglmselected2trim2 = lmer(FT ~ (mspresent + baselinemeanmsampl + nblinks + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + baselinemeanmsampl + nblinks*contrast + nblinks*eccentricity + trial*contrast, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmselected2) - BIC(loglmselected2trim2))/2) # BF = 6481731, so no interactions with baselinemsampl

boologselected <- bootstrap(model = loglmselected2, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
cilogselected = confint(boologselected, level = 0.99, type = "perc")


## now only trials with or without microsaccades/blinks
# make sure each subject has same number of trials in both models (with ms or without ms)
dfswithmsnoblinks = dfs[(dforig$baselinenms > 0) & (dforig$nblinks == 0),]
dfsnomsnoblinks = dfs[(dforig$baselinenms == 0) & (dforig$nblinks == 0),]
dforigwithmsnoblinks = dforig[(dforig$baselinenms > 0) & (dforig$nblinks == 0),]
dforignomsnoblinks = dforig[(dforig$baselinenms == 0) & (dforig$nblinks == 0),]
dfswithmssubset = data.frame()
dforigwithmssubset = data.frame()
dfsnomssubset = data.frame()
dforignomssubset = data.frame()
for (s in unique(dforig$subject)) {
  ntrials_subj_noms = nrow(dfsnomsnoblinks[(dfsnomsnoblinks$subject == s),])
  ntrials_subj_withms = nrow(dfswithmsnoblinks[(dfswithmsnoblinks$subject == s),])
  subj_dfswithms = dfswithmsnoblinks[(dfswithmsnoblinks$subject == s),]
  subj_dforigwithms = dforigwithmsnoblinks[(dforigwithmsnoblinks$subject == s),]
  if (ntrials_subj_noms > 0) {
    subj_dfsnoms = dfsnomsnoblinks[(dfsnomsnoblinks$subject == s),]
    subj_dforignoms = dforignomsnoblinks[(dforignomsnoblinks$subject == s),]
    if (ntrials_subj_withms >= ntrials_subj_noms) {
      dfswithmssubset = rbind(dfswithmssubset, subj_dfswithms[sample(nrow(subj_dfswithms), ntrials_subj_noms),])
      dfsnomssubset = rbind(dfsnomssubset, subj_dfsnoms)
      dforigwithmssubset = rbind(dforigwithmssubset, subj_dforigwithms[sample(nrow(subj_dforigwithms), ntrials_subj_noms),])
      dforignomssubset = rbind(dforignomssubset, subj_dforignoms)
    } else {
      dfswithmssubset = rbind(dfswithmssubset, subj_dfswithms)
      dfsnomssubset = rbind(dfsnomssubset, subj_dfsnoms[sample(nrow(subj_dfsnoms), ntrials_subj_withms),])
      dforigwithmssubset = rbind(dforigwithmssubset, subj_dforigwithms)
      dforignomssubset = rbind(dforignomssubset, subj_dforignoms[sample(nrow(subj_dforignoms), ntrials_subj_withms),])
    }
  }
}
# First: no microsaccades or blinks.
lmnomsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + trial:contrast + baselineretinalslip, dfsnomssubset, control=lmerControl("bobyqa"))
boonomsnoblinks <- bootstrap(model = lmnomsnoblinks, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
confint(boonomsnoblinks, level = 0.95, type = "perc")
save(boonomsnoblinks, file = "final_noms_bootstraps.Rdata")
lmnomsnoblinkstrim = lmer(FT ~ (1 | subject) + contrast + trial + trial:contrast + baselineretinalslip, dfsnomssubset, control=lmerControl("bobyqa"))
exp((BIC(lmnomsnoblinks) - BIC(lmnomsnoblinkstrim))/2) # BF = 389.7, so evidence for no eccentricity effect.
# confirm worse fit without contrast effect
lmnomsnoblinksnocontrast = lmer(FT ~ (1 | subject) + eccentricity + trial + baselineretinalslip, dfsnomssubset, control=lmerControl("bobyqa"))
exp((BIC(lmnomsnoblinks) - BIC(lmnomsnoblinksnocontrast))/2) # BF = ~0, so evidence for including contrast effect.

glmnomsnoblinks = glmer(FT ~ (1 | subject) + contrast + eccentricity + trial + trial:contrast + baselineretinalslip, dfsnomssubset, family=inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))

logdfsnomssubset = dfsnomssubset
logdfsnomssubset$FT = log(dfsnomssubset$FT)
loglmnomsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfsnomssubset, control=lmerControl("bobyqa"))


# Next: trial subset with at least 1 microsaccade, no blinks.
lmwithmsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + trial:contrast + baselineretinalslip, dfswithmssubset, control=lmerControl("bobyqa"))
boowithmsnoblinks <- bootstrap(model = lmwithmsnoblinks, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
confint(boowithmsnoblinks, level = 0.95, type = "perc")
save(boowithmsnoblinks, file = "final_withms_bootstraps.Rdata")
lmwithmsnoblinkstrim = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, dfswithmssubset, control=lmerControl("bobyqa"))
exp((BIC(lmwithmsnoblinks) - BIC(lmwithmsnoblinkstrim))/2) # BF = 96.77, so evidence for no contrast:trial effect here.
# confirm worse fit without eccentricity effect
lmwithmsnoblinksnoeccentricity = lmer(FT ~ (1 | subject) + contrast + trial + trial:contrast + baselineretinalslip, dfswithmssubset, control=lmerControl("bobyqa"))
exp((BIC(lmwithmsnoblinks) - BIC(lmwithmsnoblinksnoeccentricity))/2) # BF = ~0, so evidence for including eccentricity effect here.

glmwithmsnoblinks = glmer(FT ~ (1 | subject) + contrast + eccentricity + trial + trial:contrast + baselineretinalslip, dfswithmssubset, family=inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))

logdfswithmssubset = dfswithmssubset
logdfswithmssubset$FT = log(dfswithmssubset$FT)
loglmwithmsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfswithmssubset, control=lmerControl("bobyqa"))


# Next: trials with at least one blink, no microsaccades (only 118 trials)
# no interactions done because there are so, so few trials. And just C + E parameters.
dfsnomswithblinks = dfs[(dforig$baselinenms == 0) & (dforig$nblinks > 0),]
lmnomswithblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + trial:contrast + baselineretinalslip, dfsnomswithblinks, control=lmerControl("bobyqa"))
boonomswithblinks <- bootstrap(model = lmnomswithblinks, .f = fixef, type = "case", B = 100, resample = c(TRUE, FALSE))
confint(boonomswithblinks, level = 0.95, type = "perc")
save(boonomswithblinks, file = "final_nomswithblinks_bootstraps.Rdata")
lmnomswithblinkstrim1 = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, dfsnomswithblinks, control=lmerControl("bobyqa"))
exp((BIC(lmnomswithblinks) - BIC(lmnomswithblinkstrim1))/2) # BF = 55.26, so evidence for no contrast:trial effect.
lmnomswithblinkstrim2 = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + trial:contrast, dfsnomswithblinks, control=lmerControl("bobyqa"))
exp((BIC(lmnomswithblinks) - BIC(lmnomswithblinkstrim2))/2) # BF = 39.85, so evidence for no ocular drift effect.
lmnomswithblinkstrimfinal = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial, dfsnomswithblinks, control=lmerControl("bobyqa"))
exp((BIC(lmnomswithblinks) - BIC(lmnomswithblinkstrimfinal))/2) # BF = 2076.73, so evidence neither insig effect.
# confirm worse fit without eccentricity effect
lmnomswithblinksnoeccentricity = lmer(FT ~ (1 | subject) + contrast + trial + trial:contrast + baselineretinalslip, dfsnomswithblinks, control=lmerControl("bobyqa"))
exp((BIC(lmnomswithblinks) - BIC(lmnomswithblinksnoeccentricity))/2) # BF = 0.02, so evidence for including eccentricity effect here.

glmnomswithblinks = glmer(FT ~ (1 | subject) + contrast + eccentricity + trial + trial:contrast + baselineretinalslip, dfsnomswithblinks, family=inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
glmnomswithblinksnoeccentricity = glmer(FT ~ (1 | subject) + contrast + trial + trial:contrast + baselineretinalslip, dfsnomswithblinks, family=inverse.gaussian(link="identity"), control=glmerControl("bobyqa"))
exp((BIC(glmnomswithblinks) - BIC(glmnomswithblinksnoeccentricity))/2) # BF = 0.0007, so evidence for including eccentricity effect here.

logdfsnomswithblinks = dfsnomswithblinks
logdfsnomswithblinks$FT = log(dfsnomswithblinks$FT)
loglmnomswithblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfsnomswithblinks, control=lmerControl("bobyqa"))

# boxplot of FTs by C or E, in trials with / without microsaccades
combined = bind_rows(list("no ms" = dforignomssubset, ">1 ms" = dforigwithmssubset), .id = "source")
dfplot <- combined %>%
  convert_as_factor(subject, contrast, eccentricity)
dfplot <- dfplot %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))
dfplot <- dfplot %>%
  group_by(eccentricity, subject, source) %>%
  get_summary_stats(FT, type = "mean_sd")
plottitle ="all data"
bxp <- ggboxplot(
  dfplot, x = "eccentricity", y = "mean",
  color = "source",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15))


# for all the trials with at least one microsaccade, to test baseline ms rate or count
# But I don't think I'm going to use this analysis.
dfswithms = dfs[(dforig$baselinenms > 0),]
logdfswithms = logdfs[(dforig$baselinenms > 0),]
loglmallwithms = lmer(FT ~ (baselinemsrate + baselinemeanmsampl + contrast + eccentricity + trial + nblinks | subject) + contrast*baselinemsrate + baselinemeanmsampl + contrast*trial + eccentricity*baselinemsrate + nblinks, logdfswithms, control=lmerControl(optimizer="bobyqa"))
loglmallwithms2 = lmer(FT ~ (baselinenms + baselinemeanmsampl + contrast + eccentricity + trial + nblinks | subject) + contrast*baselinenms + contrast*trial + baselinemeanmsampl*contrast + eccentricity*baselinenms + nblinks, logdfswithms, control=lmerControl(optimizer="bobyqa"))
loglmallwithms4 = lmer(FT ~ (baselinemsrate + baselinemeanmsampl + contrast + eccentricity + trial + nblinks | subject) + contrast+ baselinemsrate + contrast*trial + baselinemeanmsampl*contrast + eccentricity*baselinemeanmsampl + nblinks, logdfswithms, control=lmerControl(optimizer="bobyqa"))
loglmallwithms5 = lmer(FT ~ (baselinemsampltotal + contrast + eccentricity + trial + nblinks | subject) + contrast*trial + baselinemsampltotal*contrast + eccentricity*baselinemsampltotal + nblinks, logdfswithms, control=lmerControl(optimizer="bobyqa"))
boowithms <- bootstrap(model = lmallwithmse, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
confint(boowithms, level = 0.975, type = "perc")

glmallwithms = glmer(FT ~ (baselinemsrate + contrast + eccentricity + trial + nblinks | subject) + contrast*baselinemsrate + contrast*trial + eccentricity*baselinemsrate + nblinks, dfswithms, family=inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))
glmallwithms2 = glmer(FT ~ (baselinenms + contrast + eccentricity + trial + nblinks | subject) + contrast*baselinenms + contrast*trial + eccentricity*baselinenms + nblinks, dfswithms, family=inverse.gaussian(link="identity"), control=glmerControl(optimizer="bobyqa"))



## test some trial subsets with numbers of microsaccades
logdfs2 = logdfs[(dforig$baselinenms > 2),]
loglmtest = lmer(FT ~ (contrast + eccentricity + baselinenms + baselinemeanmsampl + trial | subject) + contrast + eccentricity + baselinenms + baselinemeanmsampl + trial, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmtesttrim = lmer(FT ~ (contrast + eccentricity + baselinenms + baselinemeanmsampl + trial | subject) + contrast + eccentricity + baselinenms + baselinemeanmsampl, logdfs, control=lmerControl(optimizer="bobyqa"))


### - ##### - fully all trials - ##### - ###
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
df <- read.csv('./full_behavior.csv',sep = ",",stringsAsFactors = FALSE)
df <- df %>%
  convert_as_factor(subject, trialtype, block, pressed)
dforig = df

# select any rows?
df = dforig[(dforig$trialtype == 3),]

# remove trials with median fix distance > 1 deg
df = df[(df$medianfixdist < 1) & (df$maxfixdist < 1),]

# predict FT
# (no E effect in replay trials, but still C effect)
# (both effects in sharp trials)
# negative trial effect in all: later in block, FT is faster
# also interaction bw trial and contrast: later in block, contrast less effective
# I think this could relate to normalization.
# no equivalent normalization for eccentricity, because they are already normalized - 
# it's only microsaccades that have an eccentricity effect.
# Use lm7 I think, at least for main: tests both interactions with trial, and C*E.

lm0 = lmer(FT ~ (1 | subject), df)
lm1 = lmer(FT ~ (contrast | subject) + contrast, df, control=lmerControl("bobyqa"))
lm2 = lmer(FT ~ (eccentricity | subject) + eccentricity, df, control=lmerControl("bobyqa"))
lm3 = lmer(FT ~ (contrast + eccentricity | subject) + contrast + eccentricity, df, control=lmerControl("bobyqa"))
lm4 = lmer(FT ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, df, control=lmerControl("bobyqa"))
lm5 = lmer(FT ~ (contrast + eccentricity + trial | subject) + contrast*trial + eccentricity*trial, df, control=lmerControl("bobyqa"))
lm6 = lmer(FT ~ (contrast + eccentricity + trial | subject) + contrast*eccentricity + trial, df, control=lmerControl("bobyqa"))
lm7 = lmer(FT ~ (contrast + eccentricity + trial | subject) + contrast*eccentricity + contrast*trial + eccentricity*trial, df, control=lmerControl("bobyqa"))
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7)

lm3b = lmer(FT ~ (1 | subject) + contrast + eccentricity, df, control=lmerControl("bobyqa"))


exp((BIC(lm3b) - BIC(lm3d))/2)

booallnoeyes <- bootstrap(model = lm3c, .f = fixef, type = "case", B = 1000, resample = c(TRUE, FALSE))
ciallnoeyes = confint(booallnoeyes, level = 0.9917, type = "perc")
save(booallnoeyes, file = "final_allnoeyedata_bootstraps.Rdata")

# huge effect of median fixation distance in all trials
# negatively interacts with C/E effect directions in main, meaning: as median fix
# increases, C/E effects are less relevant...
lm10 = lmer(FT ~ (1 | subject) + contrast + eccentricity + medianfixdist, df)
lm11 = lmer(FT ~ (1 | subject) + contrast*medianfixdist + eccentricity*medianfixdist, df)

# predict median fixation distance from contrast/eccentricity
# higher contrast or smaller eccentricity = more eye movement?
lm = lmer(medianfixdist ~ (1 | subject) + trial + contrast + eccentricity, df)

# predict FT from stim parameters + trial type (main or sharp)
# every subject has positive sharp random slope
df = dforig
df$sharp = NaN
df[(df$trialtype == 1),]$sharp = 0
df[(df$trialtype == 3),]$sharp = 1
lmsharp = lmer(FT ~ (contrast + eccentricity + sharp | subject) + contrast + eccentricity + sharp, df, control=lmerControl("bobyqa"))
lmsharpint = lmer(FT ~ (contrast + eccentricity + sharp | subject) + contrast*sharp + eccentricity*sharp, df, control=lmerControl("bobyqa"))


# get summary stats of each subject (useful for replay trials)
df2 <- df %>%
  convert_as_factor(subject, contrast, eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))
df3 = df2 %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(FT, type = "mean")
# get subject averages across all stim conditions
df3 %>%
  group_by(subject) %>%
  get_summary_stats(mean, type = "mean")


# predict likelihood to fill-in a trial
# note: higher contrast = MORE likely! (& higher eccentricity = more likely)
# later trial = more likely.
lm0 = glmer(pressed ~ (1 | subject), df, family="binomial", control=glmerControl("bobyqa"))
lm1 = glmer(pressed ~ (contrast | subject) + contrast, df, family="binomial", control=glmerControl("bobyqa"))
lm2 = glmer(pressed ~ (eccentricity | subject) + eccentricity, df, family="binomial", control=glmerControl("bobyqa"))
lm3 = glmer(pressed ~ (contrast + eccentricity | subject) + contrast + eccentricity, df, family="binomial", control=glmerControl("bobyqa"))
lm3b = glmer(pressed ~ (contrast + eccentricity | subject) + contrast*eccentricity, df, family="binomial", control=glmerControl("bobyqa"))
lm4 = glmer(pressed ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, df, family="binomial", control=glmerControl("bobyqa"))
lm5 = glmer(pressed ~ (contrast + eccentricity + trial | subject) + contrast*trial + eccentricity*trial, df, family="binomial", control=glmerControl("bobyqa"))
lm6 = glmer(pressed ~ (contrast + eccentricity + trial | subject) + contrast*eccentricity + trial, df, family="binomial", control=glmerControl("bobyqa"))
anova(lm0, lm1, lm2, lm3, lm3b, lm4, lm5)


# plot across trial types
df2 = df
# convert for plotting
df2 <- df2 %>%
  convert_as_factor(subject, contrast, eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))
# average within subjects, separated by stim condition
df3 = df2 %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(FT, type = "median")
# normalize by subject mean, where pop. mean = still pop. mean
# doesn't really make sense if we are trying to plot the pop. mean within trialtypes.
pop_mean = mean(na.omit(dforig$FT))
for (s in unique(df$subject)) {
  subj_mean = mean(na.omit(dforig[(dforig$subject == s),]$FT))
  df3[(df3$subject == s),]$mean = pop_mean * df3[(df3$subject == s),]$mean / subj_mean
}
# or average within subjects and across stim conditions
df3 = df2 %>%
  group_by(subject, trialtype) %>%
  get_summary_stats(FT, type = "median")

# t tests or wilcoxon signed rank tests, comparing trial types
t.test(df3[df3$trialtype == 1,]$mean, df3[df3$trialtype == 3,]$mean, paired = TRUE)
wilcox.test(df3[df3$trialtype == 1,]$median, df3[df3$trialtype == 3,]$median, paired = TRUE)
wilcox.test(df3[df3$trialtype == 2,]$median, mu = 6) # test if replay FT is different from some time (e.g. 6 seconds)

# boxplot of mean FTs by C&E (averaged across trialtypes or just 1 trialtype)
plottitle ="all data"
#plottitle = "main trials"
#plottitle = "replay trials"
#plottitle = "sharp trials"
bxp <- ggboxplot(
  df3, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15))

# boxplot of mean FTs, averaged within conditions and across trialtypes
bxp <- ggboxplot(
  df3, x = "contrast", y = "mean",
  color = "contrast", palette = "rgb",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) # for contrast
#ggpar(bxp, ylim = c(-6, 5)) # normalized to 0
ggpar(bxp, ylim = c(0, 15)) # normalized to pop. mean FT

# boxplot of mean FTs, averaged within trialtypes
bxp <- ggboxplot(
  df3, x = "trialtype", y = "median",
  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
bxp = bxp + scale_x_discrete(labels = c("main", "replay", "sharp")) # for contrast
ggpar(bxp, ylim = c(0, 15))
