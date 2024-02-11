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
setwd("./") # set current working directory
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

dfall = bind_rows(list("Main" = dforig, "Replay" = dfreplay, "Sharp" = dfsharp), .id = "trialtype")



## Figures 1c, S1a, S1b: plot mean FT for every stimulus condition, 1 dot per subject
curr_df = dforig # Fig. 1c
#curr_df = dfreplay # Fig. S1a
#curr_df = dfsharp # Fig. S1b

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(FT, type = "mean")

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
  add = c("jitter"), add.params = list(alpha=.5),
)
bxp = bxp + ylab('Normalized Fading Time (s)') + xlab('Color Contrast (deg)') + ggtitle('Fading Times by Boundary Strength')
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="Eccentricity (deg)", labels=c("6", "4", "2"), guide = guide_legend(nrow=1)) # for contrast
bxp = bxp + theme(text = element_text(size = 9, family="Helvetica"), title = element_text(size = 10),
                  legend.text = element_text(size = 9), legend.title = element_text(size = 9), plot.title = element_text(size = 12, hjust=0.5))
bxp = bxp + theme(legend.justification=c(0,1), legend.position=c(0.01,1), legend.background = element_rect(fill = NA))
bxp = ggpar(bxp, ylim = c(0, 15))
#bxp = ggpar(bxp, ylim = c(0, 20)) + ggtitle('Replay Trials') # uncomment to change title to 'replay', or 'sharp'
bxp
ggsave(plot = bxp, width = 4, height = 3, dpi = 600, filename = "fig1c.pdf")
#ggsave(plot = bxp, width = 4, height = 3, dpi = 600, filename = "figs1a.pdf")
#ggsave(plot = bxp, width = 4, height = 3, dpi = 600, filename = "figs1b.pdf")



## Figure 1d: density plots by trialtype
dp <- ggplot(dfall, aes(x=FT, linetype=trialtype)) + stat_density(geom="line", linewidth=1, position="identity") + theme_classic() + 
  xlab('Fading Time (s)') + ylab('Density') + ggtitle('Fading Time Distributions') +
  scale_linetype_manual(name="Trial Type", values=c("solid", "dashed", "dotted")) + guides(linetype = guide_legend(override.aes = list(linewidth = 0.5)))
dp = dp + theme(text = element_text(size = 9, family="Helvetica"), title = element_text(size = 10), plot.title = element_text(size = 12, hjust=0.5),
                  legend.text = element_text(size = 9), legend.title = element_text(size = 9), axis.text = element_text(size = 9, colour="black"))
dp = dp + theme(legend.justification=c(1,1), legend.position=c(1,1))
dp
ggsave(plot = dp, width = 4, height = 3, dpi = 600, filename = "fig1d.pdf")



## Figure 2 c-d: boxplot of FTs by C or E, in trials with / without microsaccades
combined = bind_rows(list("0 m.s." = dforignomsnoblinks, ">1 m.s." = dforigwithmsnoblinks), .id = "source")
dfplot <- combined %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(eccentricity, subject, source) %>% # set to contrast (2c) or eccentricity (2d)
  get_summary_stats(FT, type = "mean_sd")

bxp <- ggboxplot(
  dfplot, x = "eccentricity", y = "mean", # set to contrast (2c) or eccentricity (2d)
  color = "source",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('Fading Time (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(name = "Color Contrast (deg)", labels = c("10", "20", "30")) + scale_colour_discrete(name = "Subset") # for contrast (Fig. 2c)
bxp = bxp + scale_x_discrete(name="Eccentricity (deg)", labels = c("6", "4", "2")) + scale_colour_discrete(name = "Subset") # for eccentricity (Fig. 2d)
bxp = bxp + theme(text = element_text(size = 9, family="Helvetica"), title = element_text(size = 10), legend.text = element_text(size = 9), legend.title = element_text(size = 9),
                  axis.text = element_text(size = 9, colour="black"))
bxp = bxp + theme(legend.justification=c(0,1), legend.position=c(0.01,1), legend.direction="horizontal",legend.background = element_rect(fill = NA))
bxp = ggpar(bxp, ylim = c(0, 15))
ggsave(plot = bxp, width = 4, height = 3, dpi = 600, filename = "fig2d.pdf")



## Figure 3a: mean baseline microsaccade rate for every stimulus condition, 1 dot per subject
curr_df = dforig

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(baselinemsrate, type = "mean")

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
bxp = bxp + ylab('Normalized Microsaccade Rate (Hz)') + xlab('Color Contrast (deg)') + ggtitle('Microsaccade Rate by Boundary Strength')
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="Eccentricity (deg)", labels=c("6", "4", "2"), guide = guide_legend(nrow=1)) # for contrast
bxp = bxp + theme(text = element_text(size = 9, family="Helvetica"), title = element_text(size = 10),
                  legend.text = element_text(size = 9), legend.title = element_text(size = 9), plot.title = element_text(size = 12, hjust=0.5))
bxp = bxp + theme(legend.justification=c(1,1), legend.position=c(1,1), legend.background = element_rect(fill = NA))
bxp = ggpar(bxp, ylim = c(0, 2.2))
ggsave(plot = bxp, width = 4, height = 3, dpi = 600, filename = "fig3a.pdf")




## Figure 3b: mean immobilization duration for every stimulus condition, 1 dot per subject
curr_df = dforig[!is.na(dforig$immobduration),]

dfplot <- curr_df %>%
  convert_as_factor(subject, contrast, eccentricity) %>%
  reorder_levels("eccentricity", order=c("6", "4", "2")) %>%
  group_by(subject, contrast, eccentricity) %>%
  get_summary_stats(immobduration, type = "mean")

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
bxp = bxp + ylab('Normalized Immobilization Time (s)') + xlab('Color Contrast (deg)') + ggtitle('Immobilization Time by Boundary Strength')
bxp = bxp + scale_x_discrete(labels = c("10", "20", "30")) + scale_colour_discrete(name="Eccentricity (deg)", labels=c("6", "4", "2"), guide = guide_legend(nrow=1)) # for contrast
bxp = bxp + theme(text = element_text(size = 9, family="Helvetica"), title = element_text(size = 10),
                  legend.text = element_text(size = 9), legend.title = element_text(size = 9), plot.title = element_text(size = 12, hjust=0.5))
bxp = bxp + theme(legend.justification=c(0,1), legend.position=c(0.01,1), legend.background = element_rect(fill = NA))
bxp = ggpar(bxp, ylim = c(0, 4.2))
ggsave(plot = bxp, width = 4, height = 3, dpi = 600, filename = "fig3b.pdf")



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