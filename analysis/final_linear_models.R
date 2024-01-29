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
logdfs$immobduration = log(dfs$immobduration)

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
logdfswithms = dfswithms
logdfswithms$baselinemsrate = log(logdfswithms$baselinemsrate)
lm6log = lmer(baselinemsrate ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, logdfswithms, control=lmerControl("bobyqa"))
boo6 <- bootstrap(model = lm6log, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo6, level = 0.983333, type = "perc")
save(boo6, file = "baselinerate_bootstrap.Rdata")

# look into "immobilization duration"
# smaller eccentricity = last microsaccade has to happen earlier, or else it'll be too strong
# and prevent filling.

# use log here too: immobduration doesn't have as obvious an interpretation
# conclusions are the same as with glm, but better convergence with log.
# Though: can't directly interpret the parameter estimate (i.e., in milliseconds)
lm5log = lmer(immobduration ~ (contrast + eccentricity + immobmsamplitude + trial | subject) + contrast + eccentricity + immobmsamplitude + trial, logdfs, control=lmerControl("bobyqa"))
boo5 <- bootstrap(model = lm5log, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boo5, level = 0.983333, type = "perc")
save(boo5, file = "immob_bootstrap.Rdata")
# note below shows: immob duration does not depend with FT (one prediction that microsaccades are more or less effective with adaptation)
lm5logstar = lmer(immobduration ~ (contrast + eccentricity + immobmsamplitude + trial | subject) + contrast + eccentricity + immobmsamplitude + trial + FT, logdfs, control=lmerControl("bobyqa"))

# permutation tests by randomizing contrast / eccentricity
nperm = 10


# get summary stats
dfimmobsummary <- dforig %>%
  convert_as_factor(eccentricity) %>%
  group_by(subject, eccentricity) %>%
  get_summary_stats(immobduration, type = "mean") %>%
  group_by(eccentricity) %>%
  get_summary_stats(mean, type = "mean_se")

# get summary stats of immobilizing microsaccade amplitude
dfamplsummary <- dforig %>%
  convert_as_factor(eccentricity) %>%
  group_by(subject, eccentricity) %>%
  get_summary_stats(immobmsamplitude, type = "mean") %>%
  group_by(eccentricity) %>%
  get_summary_stats(mean, type = "mean_se")

## 3 step model creation. Fcenter & rescaled all continuous variables
# except for FT, scale it but do NOT center
# Step 1: model with no eye movements.
loglmnoeyes = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial + eccentricity*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmnoeyestrim1 = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*trial + eccentricity*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmnoeyes) - BIC(loglmnoeyestrim1))/2) # BF = 1230, so evidence for removing contrast*eccentricity
loglmnoeyestrim2 = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmnoeyes) - BIC(loglmnoeyestrim2))/2) # BF = 218, so evidence for removing eccentricity*trial
loglmnoeyestrimfinal = lmer(FT ~ (trial + contrast + eccentricity | subject) + eccentricity + contrast*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmnoeyes) - BIC(loglmnoeyestrimfinal))/2) # BF = 305582, so evidence for removing both insignificant parameters.
loglmnoeyestrim3 = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + eccentricity*trial, logdfs, control=lmerControl(optimizer="bobyqa"))

boonoeyes <- bootstrap(model = loglmnoeyes, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boonoeyes, level = 0.9917, type = "perc")
save(boonoeyes, file = "noeyes_bootstrap.Rdata")

# Step 2: use selected parameters from 1, and add eye movement main effects.
loglmmain = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity  + baselineretinalslip | subject) + mspresent + nblinks + baselineretinalslip + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmmaintrim1 = lmer(FT ~ (mspresent + nblinks + baselineretinalslip + trial + contrast + eccentricity | subject) + mspresent + nblinks + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmmain) - BIC(loglmmaintrim1))/2) # BF = 676, so evidence for removing the fixed baselineretinalslip effect
loglmmaintrim2 = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity | subject) + mspresent + nblinks + baselineretinalslip + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmmain) - BIC(loglmmaintrim2))/2) # BF ~0, so no evidence for removing the random baselineretinalslip effect
loglmmaintrimfinal = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity | subject) + mspresent + nblinks + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmmain) - BIC(loglmmaintrimfinal))/2) # BF ~0, so evidence for keeping at least 1 of the removed effects.
# because better fit without fixed effect, no reason to keep a random effect. So get rid of retinal slip entirely.
boologmain <- bootstrap(model = loglmmain, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
cilogmain = confint(boologselected, level = 0.9929, type = "perc")
save(boologmain, file = "main_bootstrap.Rdata")

# Step 3: use selected parameters from 2, and add interactions bw C/E and eye movements. Bonferroni-corrected interaction selection.
loglmselected = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, logdfs, control=lmerControl(optimizer="bobyqa"))
loglmselectedtrim1 = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity | subject) + mspresent*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, logdfs, control=lmerControl(optimizer="bobyqa"))
exp((BIC(loglmselected) - BIC(loglmselectedtrim1))/2) # BF = 2012, so evidence for removing contrast*mspresent
boologselected <- bootstrap(model = loglmselected, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
cilogselected = confint(boologselected, level = 0.995, type = "perc")
save(boologselected, file = "selected_bootstrap.Rdata")


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
logdfsnomssubset = dfsnomssubset
logdfsnomssubset$FT = log(dfsnomssubset$FT)
loglmnomsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfsnomssubset, control=lmerControl("bobyqa"))
# confirm better fit without eccentricity effect, but worse without contrast effect
loglmnomsnoblinksnoeccentricity = lmer(FT ~ (1 | subject) + contrast + trial + baselineretinalslip, logdfsnomssubset, control=lmerControl("bobyqa"))
exp((BIC(loglmnomsnoblinks) - BIC(loglmnomsnoblinksnoeccentricity))/2) # BF = 792, so evidence for removing eccentricity effect
loglmnomsnoblinksnocontrast = lmer(FT ~ (1 | subject) + eccentricity + trial + baselineretinalslip, logdfsnomssubset, control=lmerControl("bobyqa"))
exp((BIC(loglmnomsnoblinks) - BIC(loglmnomsnoblinksnocontrast))/2) # BF = ~0, so evidence for including contrast effect

boonomsnoblinks <- bootstrap(model = loglmnomsnoblinks, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boonomsnoblinks, level = 0.95, type = "perc")
save(boonomsnoblinks, file = "nomssubset_bootstraps.Rdata")

# Next: trial subset with at least 1 microsaccade, no blinks.
logdfswithmssubset = dfswithmssubset
logdfswithmssubset$FT = log(dfswithmssubset$FT)
loglmwithmsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfswithmssubset, control=lmerControl("bobyqa"))
# confirm worse fit without eccentricity effect OR contrast effect
loglmwithmsnoblinksnoeccentricity = lmer(FT ~ (1 | subject) + contrast + trial + baselineretinalslip, logdfswithmssubset, control=lmerControl("bobyqa"))
exp((BIC(loglmwithmsnoblinks) - BIC(loglmwithmsnoblinksnoeccentricity))/2) # BF = ~0, so evidence for including eccentricity effect
loglmwithmsnoblinksnocontrast = lmer(FT ~ (1 | subject) + eccentricity + trial + baselineretinalslip, logdfswithmssubset, control=lmerControl("bobyqa"))
exp((BIC(loglmwithmsnoblinks) - BIC(loglmwithmsnoblinksnocontrast))/2) # BF = ~0, so evidence for including contrast effect

boowithmsnoblinks <- bootstrap(model = loglmwithmsnoblinks, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boowithmsnoblinks, level = 0.95, type = "perc")
save(boowithmsnoblinks, file = "withmssubset_bootstraps.Rdata")

# Next: trials with at least one blink, no microsaccades (only 118 trials)
# no interactions done because there are so, so few trials. And just C + E parameters.
dfsnomswithblinks = dfs[(dforig$baselinenms == 0) & (dforig$nblinks > 0),]
dforignomswithblinks = dforig[(dforig$baselinenms == 0) & (dforig$nblinks > 0),]
logdfsnomswithblinks = dfsnomswithblinks
logdfsnomswithblinks$FT = log(dfsnomswithblinks$FT)
loglmnomswithblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfsnomswithblinks, control=lmerControl("bobyqa"))
# confirm worse fit without eccentricity effect OR contrast effect
loglmnomswithblinksnoeccentricity = lmer(FT ~ (1 | subject) + contrast + trial + baselineretinalslip, logdfsnomswithblinks, control=lmerControl("bobyqa"))
exp((BIC(loglmnomswithblinks) - BIC(loglmnomswithblinksnoeccentricity))/2) # BF = 0.02, so evidence for including eccentricity effect
loglmnomswithblinksnocontrast = lmer(FT ~ (1 | subject) + eccentricity + trial + baselineretinalslip, logdfsnomswithblinks, control=lmerControl("bobyqa"))
exp((BIC(loglmnomswithblinks) - BIC(loglmnomswithblinksnocontrast))/2) # BF = 68, so evidence for including contrast effect

boonomswithblinks <- bootstrap(model = loglmnomswithblinks, .f = fixef, type = "case", B = 10000, resample = c(TRUE, FALSE))
confint(boonomswithblinks, level = 0.95, type = "perc")
save(boonomswithblinks, file = "nomswithblinks_bootstraps.Rdata")


# boxplot of FTs by C or E, in trials with / without microsaccades
#combined = bind_rows(list("no ms" = dforignomssubset, ">1 ms" = dforigwithmssubset), .id = "source")
combined = bind_rows(list("no ms" = dforignomsnoblinks, ">1 ms" = dforigwithmsnoblinks), .id = "source")
dfplot <- combined %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(eccentricity, subject, source) %>%
  get_summary_stats(FT, type = "mean_sd")

bxp <- ggboxplot(
  dfplot, x = "eccentricity", y = "mean",
  color = "source",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(name = "contrast (deg)", labels = c("10", "20", "30")) + scale_colour_discrete(name = "subset") # for eccentricity
bxp = bxp + scale_x_discrete(name="eccentricity (deg)", labels = c("6", "4", "2")) + scale_colour_discrete(name = "subset") # for contrast
ggpar(bxp, ylim = c(0, 15))


# box plot of FTs by E, in trials without microsaccades but with blinks
dfplot <- dforignomswithblinks %>%
  convert_as_factor(subject, contrast) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(eccentricity, subject) %>%
  get_summary_stats(FT, type = "mean")

bxp <- ggboxplot(
  dfplot, x = "eccentricity", y = "mean",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
bxp = bxp + scale_x_discrete(name="eccentricity (deg)", labels = c("6", "4", "2"))# for contrast
ggpar(bxp, ylim = c(0, 15))




### ---- load control trials ---- ###
# replay trials
dfreplay <- read.csv('./trial_data_replay_withblinks.csv',sep = ",",stringsAsFactors = FALSE)
dfreplay$FT = dfreplay$RT
dfreplay$baselinemsrate = dfreplay$baselinenms / dfreplay$FT
dfreplay$blinkrate = dfreplay$nblinks / dfreplay$FT
dfreplay <- dfreplay %>%
  convert_as_factor(subject, mspresent, earlymspresent, latemspresent, blinkpresent, session)
dfreplayorig = dfreplay
# add scaled_eccentricity
A = 0.063
l0 = 36.54 * A
dfreplay$scaled_eccentricity = (log(dfreplay$eccentricity) - l0) / A
dfreplay$scaled_eccentricity = 1 / (A * dfreplay$eccentricity)

# rescale and zero-mean all continuous predictors (but not FT)
numcols <- sapply(dfreplay, is.numeric)
dfsreplay <- dfreplay
dfsreplay[,numcols] <- scale(dfsreplay[,numcols])
dfreplaynocenter <- dfreplay
dfreplaynocenter[,numcols] <- scale(dfreplaynocenter[,numcols], center = FALSE, scale = apply(dfreplaynocenter[,numcols], 2, sd, na.rm = TRUE))
dfsreplay$FT = dfreplaynocenter$FT
dfsreplay$immobduration = dfreplaynocenter$immobduration
#dfsreplay = dfreplaynocenter
logdfsreplay = dfsreplay
logdfsreplay$FT = log(dfsreplay$FT)


# sharp trials
dfsharp <- read.csv('./trial_data_sharp_withblinks.csv',sep = ",",stringsAsFactors = FALSE)
dfsharp$FT = dfsharp$RT
dfsharp$baselinemsrate = dfsharp$baselinenms / dfsharp$FT
dfsharp$blinkrate = dfsharp$nblinks / dfsharp$FT
dfsharp <- dfsharp %>%
  convert_as_factor(subject, mspresent, earlymspresent, latemspresent, blinkpresent, session)
dfsharporig = dfsharp
# add scaled_eccentricity
A = 0.063
l0 = 36.54 * A
dfsharp$scaled_eccentricity = (log(dfsharp$eccentricity) - l0) / A
dfsharp$scaled_eccentricity = 1 / (A * dfsharp$eccentricity)

# rescale and zero-mean all continuous predictors (but not FT)
numcols <- sapply(dfsharp, is.numeric)
dfssharp <- dfsharp
dfssharp[,numcols] <- scale(dfssharp[,numcols])
dfsharpnocenter <- dfsharp
dfsharpnocenter[,numcols] <- scale(dfsharpnocenter[,numcols], center = FALSE, scale = apply(dfsharpnocenter[,numcols], 2, sd, na.rm = TRUE))
dfssharp$FT = dfsharpnocenter$FT
dfssharp$immobduration = dfsharpnocenter$immobduration
#dfssharp = dfsharpnocenter
logdfssharp = dfssharp
logdfssharp$FT = log(dfssharp$FT)


## linear model (don't use)
loglmsharp = lmer(FT ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, logdfssharp, control=lmerControl(optimizer="bobyqa"))
loglmreplay = lmer(FT ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, logdfsreplay, control=lmerControl(optimizer="bobyqa"))



### ---- summary stats and plotting ---- ###

#### plot mean FT for every stimulus condition, 1 dot per subject
curr_df = dfsharp

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(FT, type = "mean")

# get subject mean FT, averaged across all stim conditions
dfplot %>%
  group_by(subject) %>%
  get_summary_stats(mean, type = "mean")

# optional: normalize by subject mean, where pop. mean = still pop. mean?
pop_mean = mean(na.omit(curr_df$FT))
for (s in unique(dfplot$subject)) {
  subj_mean = mean(na.omit(curr_df[(curr_df$subject == s),]$FT))
  dfplot[(dfplot$subject == s),]$mean = pop_mean * dfplot[(dfplot$subject == s),]$mean / subj_mean
}

# boxplot of mean FTs by C&E 
#plottitle = "main trials"
#plottitle = "replay trials"
#plottitle = "sharp trials"
bxp <- ggboxplot(
  dfplot, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
#  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('norm FT (s)') + xlab('color contrast (deg)') + theme(text = element_text(size = 20, family="Helvetica"))
#bxp = bxp + scale_x_discrete(labels = c("6", "4", "2")) + scale_colour_discrete(labels=c("10", "20", "30")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="eccentricity (deg)", labels=c("6", "4", "2")) # for contrast
ggpar(bxp, ylim = c(0, 15))


#### plot mean baseline microsaccade rate for every stimulus condition, 1 dot per subject
curr_df = dforig

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(baselinemsrate, type = "mean")

# get subject mean msrate, averaged across all stim conditions
dfplot %>%
  group_by(subject) %>%
  get_summary_stats(mean, type = "mean")

# optional: normalize by subject mean, where pop. mean = still pop. mean?
pop_mean = mean(na.omit(curr_df$baselinemsrate))
for (s in unique(dfplot$subject)) {
  subj_mean = mean(na.omit(curr_df[(curr_df$subject == s),]$baselinemsrate))
  dfplot[(dfplot$subject == s),]$mean = pop_mean * dfplot[(dfplot$subject == s),]$mean / subj_mean
}

# boxplot of mean baseline ms rate by C&E
bxp <- ggboxplot(
  dfplot, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  #  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('norm microsaccade rate (Hz)') + xlab('color contrast (deg)') + theme(text = element_text(size = 20, family="Helvetica"))
#bxp = bxp + scale_x_discrete(labels = c("6", "4", "2")) + scale_colour_discrete(labels=c("10", "20", "30")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="eccentricity (deg)", labels=c("6", "4", "2")) # for contrast
ggpar(bxp, ylim = c(0, 2))


#### plot mean immobilization duration for every stimulus condition, 1 dot per subject
curr_df = dforig[!is.na(dforig$immobduration),]

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(immobduration, type = "mean")

# get subject mean msrate, averaged across all stim conditions
dfplot %>%
  group_by(subject) %>%
  get_summary_stats(mean, type = "mean")

# optional: normalize by subject mean, where pop. mean = still pop. mean?
pop_mean = mean(na.omit(curr_df$immobduration))
for (s in unique(dfplot$subject)) {
  subj_mean = mean(na.omit(curr_df[(curr_df$subject == s),]$immobduration))
  dfplot[(dfplot$subject == s),]$mean = pop_mean * dfplot[(dfplot$subject == s),]$mean / subj_mean
}

# boxplot of mean baseline ms rate by C&E
bxp <- ggboxplot(
  dfplot, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  #  title = plottitle,
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('norm immobilization (s)') + xlab('color contrast (deg)') + theme(text = element_text(size = 20, family="Helvetica"))
#bxp = bxp + scale_x_discrete(labels = c("6", "4", "2")) + scale_colour_discrete(labels=c("10", "20", "30")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="eccentricity (deg)", labels=c("6", "4", "2")) # for contrast
ggpar(bxp, ylim = c(0, 4))




# # plot FT as function of trial, after regressing out everything else?
#loglmmost = lmer(FT ~ (mspresent + nblinks + contrast + eccentricity | subject) + contrast + mspresent*eccentricity + nblinks*eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
#dftrialres = dforig
#dftrialres$FT = residuals(loglmmost)
#pp <- ggplot(dftrialres, aes(x=trial, y=FT)) + geom_point(size=2, shape=23)



# t tests or wilcoxon signed rank tests, comparing trial types
t.test(df3[df3$trialtype == 1,]$mean, df3[df3$trialtype == 3,]$mean, paired = TRUE)
wilcox.test(df3[df3$trialtype == 1,]$median, df3[df3$trialtype == 3,]$median, paired = TRUE)
wilcox.test(df3[df3$trialtype == 2,]$median, mu = 6) # test if replay FT is different from some time (e.g. 6 seconds)



# boxplot of mean FTs, averaged within trialtypes
dfall = bind_rows(list("main" = dforig, "replay" = dfreplay, "sharp" = dfsharp), .id = "trialtype")
dfplot <- dfall %>%
  group_by(subject, trialtype) %>%
  get_summary_stats(FT, type = "mean")

bxp <- ggboxplot(
  dfplot, x = "trialtype", y = "mean",
  color = "variable", palette = "black", # just placeholder to make plot shape the same
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + xlab('trial type') + theme(text = element_text(size = 20, family = "Helvetica"))
bxp = bxp + scale_x_discrete(labels = c("main", "replay", "sharp"))
ggpar(bxp, ylim = c(0, 15))


## normalize FTs by subject, then plot distributions across all trials, by trialtype
#dfallnorm <- dfall
#pop_mean = mean(na.omit(dfall$FT))
#for (s in unique(dfall$subject)) {
#  subj_mean = mean(na.omit(dfall[(dfall$subject == s),]$FT))
#  dfallnorm[(dfall$subject == s),]$FT = pop_mean * dfall[(dfall$subject == s),]$FT / subj_mean
#}

# density plots by trialtype
bxp <- ggplot(dfall, aes(x=FT, linetype=trialtype)) + stat_density(geom="line", linewidth=1, position="identity") + theme_classic() + xlab('FT (s)') +
  theme(text = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 20, colour="black"), axis.text.y = element_text(size = 20, colour="black")) + 
  scale_linetype_manual(name="trial type", values=c("solid", "dashed", "dotted")) + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5)))
ggpar(bxp, xlim = c(0, 20))






###### Summary data. Load all trials
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
df <- read.csv('./full_behavior.csv',sep = ",",stringsAsFactors = FALSE)
df <- df %>%
  convert_as_factor(subject, trialtype, block, pressed)
levels(df$trialtype) <- c("main", "replay", "sharp")
logdf = df
logdf$FT = log(df$FT)

dfsubject <- df %>%
  group_by(subject, trialtype) %>%
  get_summary_stats(FT, type = "mean")

dfsum <- dfsubject %>%
  group_by(trialtype) %>%
  get_summary_stats(mean, type = "mean_se")

dfsubjectplot <- dfsubject[(dfsubject$trialtype == "main") | (dfsubject$trialtype == "sharp"),]
# plot subject averages in 2 conditions
bxp <- ggpaired(dfsubjectplot, x = "trialtype", y = "mean",
                id = "subject", line.color = "gray",
                ylab = "FT (s)", xlab = "trial type")
bxp# + stat_compare_means(paired=TRUE)

# wilcoxon tests
wilcox.test(dfsubject[dfsubject$trialtype == "main",]$mean, dfsubject[dfsubject$trialtype == "replay",]$mean, paired = TRUE)
wilcox.test(dfsubject[dfsubject$trialtype == "main",]$mean, dfsubject[dfsubject$trialtype == "sharp",]$mean, paired = TRUE)
# sharp: V = 0, meaning the lowest sum of ranks is 0 because all ranks are negative. There are no main values larger than sharp value.
# test statistic W is be the lowest of W- (sum of negative ranks) or W+ (sum of positive ranks)
# while V is ONLY W+ 

# COMPUTE W:
a <- dfsubject[dfsubject$trialtype == "main",]$mean
b <- dfsubject[dfsubject$trialtype == "replay",]$mean
diff <- c(a - b) # vector of paired differences
diff <- diff[ diff!=0 ] # delete all differences equal to zero
diff.rank <- rank(abs(diff)) # get the absolute value of the ranks of the differences
diff.rank.sign <- diff.rank * sign(diff) # now add sign back to the absolute ranks, based on which group (a or b) is larger
ranks.pos <- sum(diff.rank.sign[diff.rank.sign > 0]) # sum of ranks assigned to positive differences (W+), a > b, i.e., V statistic
ranks.neg <- -sum(diff.rank.sign[diff.rank.sign < 0]) # sum of ranks assigned to negative differences (W-), a < b
W <- min(ranks.pos, ranks.neg)
W

wilcox.test(dfsubject[dfsubject$trialtype == 2,]$mean, mu = 6) # test if replay FT is different from some time (e.g. 6 seconds)

# plot all trialtypes
bxp <- ggboxplot(
  dfsubject, x = "trialtype", y = "mean",
  order = c("replay", "main", "sharp"),
  color = "variable", palette = "black" # just placeholder to make plot shape the same
  #add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + geom_point(aes(group=subject),shape=21,color="black",fill="gray") + geom_line(aes(group=subject), color="gray")
bxp = bxp + ylab('FT (s)') + xlab('trial type') + theme(text = element_text(size = 20, family = "Helvetica"))
#bxp = bxp + scale_x_discrete(labels = c("main", "replay", "sharp"))
ggpar(bxp, ylim = c(0, 15))
