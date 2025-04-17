# load packages
library(lme4)
library(sjPlot)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(lmeresampler)
library(performance)
library(qqplotr)
library(ggeffects)

## ---- load individual trials ---- ##
setwd("/export03/data/mlevin/UI_microsaccades/bhv_data")
df <- read.csv('./trial_data_withblinks.csv',sep = ",",stringsAsFactors = FALSE)
df$FT = df$RT
df$baselinemsrate = df$baselinenms / df$FT
df <- df %>%
  convert_as_factor(subject, mspresent, earlymspresent, latemspresent, blinkpresent, session)
# add scaled_eccentricity
A = 0.063
l0 = 36.54 * A
df$scaled_eccentricity = (log(df$eccentricity) - l0) / A
df$scaled_eccentricity = 1 / (A * df$eccentricity)
dforig = df


### create sub datasets for plots

# trials with or without microsaccades/blinks
# make sure each subject has same number of trials in both models (with ms or without ms)
dforigwithmsnoblinks = dforig[(dforig$baselinenms > 0) & (dforig$nblinks == 0),]
dforignomsnoblinks = dforig[(dforig$baselinenms == 0) & (dforig$nblinks == 0),]

# load control trials
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

dfall = bind_rows(list("main" = dforig, "replay" = dfreplay, "sharp" = dfsharp), .id = "trialtype")



## Figures 1c, S1a, S1b: plot mean FT for every stimulus condition, 1 dot per subject
curr_df = dforig # Fig. 1c
#curr_df = dfreplay # Fig. S1a
#curr_df = dfsharp # Fig. S1b

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(FT, type = "mean")

# get subject mean FT, averaged across all stim conditions
dfplot %>%
  group_by(subject) %>%
  get_summary_stats(mean, type = "mean")

# normalize by subject mean, where pop. mean = still pop. mean
pop_mean = mean(na.omit(curr_df$FT))
for (s in unique(dfplot$subject)) {
  subj_mean = mean(na.omit(curr_df[(curr_df$subject == s),]$FT))
  dfplot[(dfplot$subject == s),]$mean = pop_mean * dfplot[(dfplot$subject == s),]$mean / subj_mean
}

# boxplot of mean FTs by C&E 
bxp <- ggboxplot(
  dfplot, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('norm FT (s)') + xlab('color contrast (deg)') + theme(text = element_text(size = 20, family="Helvetica"))
#bxp = bxp + scale_x_discrete(labels = c("6", "4", "2")) + scale_colour_discrete(labels=c("10", "20", "30")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="eccentricity (deg)", labels=c("6", "4", "2")) # for contrast
ggpar(bxp, ylim = c(0, 15))


## Figure 1d: density plots by trialtype
bxp <- ggplot(dfall, aes(x=FT, linetype=trialtype)) + stat_density(geom="line", linewidth=1, position="identity") + theme_classic() + xlab('FT (s)') +
  theme(text = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 20, colour="black"), axis.text.y = element_text(size = 20, colour="black")) + 
  scale_linetype_manual(name="trial type", values=c("solid", "dashed", "dotted")) + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5)))
ggpar(bxp, xlim = c(0, 20))


## Figure 2 c-d: boxplot of FTs by C or E, in trials with / without microsaccades
combined = bind_rows(list("no ms" = dforignomsnoblinks, ">1 ms" = dforigwithmsnoblinks), .id = "source")
dfplot <- combined %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(eccentricity, subject, source) %>%
  get_summary_stats(FT, type = "mean_sd")

bxp <- ggboxplot(
  dfplot, x = "eccentricity", y = "mean", # set to eccentricity or contrast
  color = "source",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('FT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(name = "contrast (deg)", labels = c("10", "20", "30")) + scale_colour_discrete(name = "subset") # for eccentricity
bxp = bxp + scale_x_discrete(name="eccentricity (deg)", labels = c("6", "4", "2")) + scale_colour_discrete(name = "subset") # for contrast
ggpar(bxp, ylim = c(0, 15))




## Figure 3a: mean baseline microsaccade rate for every stimulus condition, 1 dot per subject
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

# normalize by subject mean, where pop. mean = pop. mean
pop_mean = mean(na.omit(curr_df$baselinemsrate))
for (s in unique(dfplot$subject)) {
  subj_mean = mean(na.omit(curr_df[(curr_df$subject == s),]$baselinemsrate))
  dfplot[(dfplot$subject == s),]$mean = pop_mean * dfplot[(dfplot$subject == s),]$mean / subj_mean
}

# boxplot of mean baseline ms rate by C&E
bxp <- ggboxplot(
  dfplot, x = "contrast", y = "mean",
  color = "eccentricity", palette = "rgb",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('norm microsaccade rate (Hz)') + xlab('color contrast (deg)') + theme(text = element_text(size = 20, family="Helvetica"))
#bxp = bxp + scale_x_discrete(labels = c("6", "4", "2")) + scale_colour_discrete(labels=c("10", "20", "30")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="eccentricity (deg)", labels=c("6", "4", "2")) # for contrast
ggpar(bxp, ylim = c(0, 2))


## Figure 3b: mean immobilization duration for every stimulus condition, 1 dot per subject
curr_df = dforig[!is.na(dforig$immobduration),]

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(immobduration, type = "mean")

# get subject mean, averaged across all stim conditions
dfplot %>%
  group_by(subject) %>%
  get_summary_stats(mean, type = "mean")

# normalize by subject mean, where pop. mean = pop. mean
pop_mean = mean(na.omit(curr_df$immobduration))
for (s in unique(dfplot$subject)) {
  subj_mean = mean(na.omit(curr_df[(curr_df$subject == s),]$immobduration))
  dfplot[(dfplot$subject == s),]$mean = pop_mean * dfplot[(dfplot$subject == s),]$mean / subj_mean
}

# boxplot of mean immobilization by C&E
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


###### Summary data. Load all trials, without excluding by eye movement criteria
df <- read.csv('./full_behavior.csv',sep = ",",stringsAsFactors = FALSE)
df <- df %>%
  convert_as_factor(subject, trialtype, block, pressed)
levels(df$trialtype) <- c("main", "replay", "sharp")

dfsubject <- df %>%
  group_by(subject, trialtype) %>%
  get_summary_stats(FT, type = "mean")

dfsum <- dfsubject %>%
  group_by(trialtype) %>%
  get_summary_stats(mean, type = "mean_se")

# wilcoxon test
wilcox.test(dfsubject[dfsubject$trialtype == "main",]$mean, dfsubject[dfsubject$trialtype == "sharp",]$mean, paired = TRUE)
# sharp: V = 0, meaning the lowest sum of ranks is 0 because all ranks are negative. There are no main values larger than sharp value.
# test statistic W is be the lowest of W- (sum of negative ranks) or W+ (sum of positive ranks)
# while V is ONLY W+ 

# COMPUTE W:
a <- dfsubject[dfsubject$trialtype == "main",]$mean
b <- dfsubject[dfsubject$trialtype == "sharp",]$mean
diff <- c(a - b) # vector of paired differences
diff <- diff[ diff!=0 ] # delete all differences equal to zero
diff.rank <- rank(abs(diff)) # get the absolute value of the ranks of the differences
diff.rank.sign <- diff.rank * sign(diff) # now add sign back to the absolute ranks, based on which group (a or b) is larger
ranks.pos <- sum(diff.rank.sign[diff.rank.sign > 0]) # sum of ranks assigned to positive differences (W+), a > b, i.e., V statistic
ranks.neg <- -sum(diff.rank.sign[diff.rank.sign < 0]) # sum of ranks assigned to negative differences (W-), a < b
W <- min(ranks.pos, ranks.neg)
W