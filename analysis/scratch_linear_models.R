# load package
library(lme4) # for linear mixed models (not survival)
library(sjPlot)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(mediation)


# load data as df_full
#setwd("/export04/data/mlevin/Uniformity_Illusion/bhv_data")
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
df <- read.csv('./minrates.csv',sep = ",",stringsAsFactors = FALSE)
#df <- read.csv('./noms_minrates.csv',sep = ",",stringsAsFactors = FALSE)
#df <- read.csv('./onlyms_minrates.csv',sep = ",",stringsAsFactors = FALSE)
df <- df %>%
  convert_as_factor(subject)

# fit mixed effect models
# LM5B IS BEST
lm0 = lmer(msratemin ~ (1 | subject), df)
lm1 = lmer(msratemin ~ (1 | subject) + contrast, df)
lm2 = lmer(msratemin ~ (1 | subject) + eccentricity, df)
lm3 = lmer(msratemin ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(msratemin ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(msratemin ~ (1 | subject) + contrast*eccentricity, df)
lm5 = lmer(msratemin ~ (1 | subject) + contrast + eccentricity + baselinerate, df)
lm5a = lmer(msratemin ~ (1 + baselinerate | subject) + contrast + eccentricity + baselinerate, df)
lm5b = lmer(msratemin ~ (1 + baselinerate | subject) + contrast + scaled_eccentricity + baselinerate, df)
lm6 = lmer(msratemin ~ (1 + baselinerate | subject) + contrast*scaled_eccentricity + baselinerate , df)
lm7 = lmer(msratemin ~ (1 | subject) + baselinerate*contrast + eccentricity, df)
lm8 = lmer(msratemin ~ (1 | subject) + baselinerate*contrast*eccentricity, df)

anova(lm0, lm1, lm2, lm3, lm3b, lm4, test = "Chisq")
anova(lm0, lm1, lm2, lm3, lm3b, lm4, lm5, lm5a, lm5b, lm6, lm7, lm8, test = "Chisq")
anova(lm0, lm1, lm2, lm3, lm4, lm5, test = "Chisq")
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, test = "Chisq")


# calculate Bayes factor
exp((BIC(lm5a) - BIC(lm5b))/2)

# explore predicting RT from msratemin.
lm0 = lmer(RT ~ (1 | subject) + msratemin, df)
lm0b = lmer(RT ~ (1 + msratemin | subject) + msratemin, df)
lm1 = lmer(RT ~ (1 | subject) + msratemin + baselinerate, df)
lm1b = lmer(RT ~ (1 + baselinerate | subject) + msratemin + baselinerate, df)
lm1c = lmer(RT ~ (1 + baselinerate + msratemin | subject) + msratemin + baselinerate, df)
lm2 = lmer(RT ~ (1 | subject) + msratemin + contrast + eccentricity, df)
lm2c = lmer(RT ~ (1 + msratemin | subject) + msratemin + contrast + eccentricity, df)
lm3 = lmer(RT ~ (1 | subject) + msratemin + baselinerate + contrast + eccentricity, df)
lm3b = lmer(RT ~ (1 + baselinerate | subject) + msratemin + baselinerate + contrast + eccentricity, df)
lm3c = lmer(RT ~ (1 + baselinerate + msratemin | subject) + msratemin + baselinerate + contrast + eccentricity, df)

anova(lm0, lm0b, lm1, lm1b, lm1c, lm2, lm2c, lm3, lm3b, lm3c, test = "Chisq")

# predict RT from just contrast and eccentricity (compare ms to noms)
# Result: noms: best model is 1: only contrast predicts.
# onlyms: best model is 3/3b: both contrast & eccentricity predict.
lm0 = lmer(RT ~ (1 | subject), df)
lm1 = lmer(RT ~ (1 | subject) + contrast, df)
lm2 = lmer(RT ~ (1 | subject) + eccentricity, df)
lm2b = lmer(RT ~ (1 | subject) + scaled_eccentricity, df)
lm3 = lmer(RT ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(RT ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(RT ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(RT ~ (1 | subject) + contrast*scaled_eccentricity, df)

anova(lm0, lm1, lm2, lm2b, lm3, lm3b, lm4, lm4b, test = "Chisq")

# Does minimum rate (baseline-corrected / relative) mediate contrast or eccentricity > RT?
# using best-fitting min rate model (lm5b above)
# NO because minrate only weakly correlates with RT, only if you don't include con/ecc and baselinerate.
# because RT and minrate both reflect microsaccade influence, and minrate is noisy
# only correlates with model: lmer(RT ~ (1 + baselinerate | subject) + msratemin + baselinerate, df)
# the second you add either contrast or eccentricity, corrected-msratemin no longer predicts RT.
lmmed = lmer(msratemin ~ (1 + baselinerate | subject) + contrast + scaled_eccentricity + baselinerate, df) # mediation model
lmout = lmer(RT ~ (1 + baselinerate | subject) + contrast + baselinerate + scaled_eccentricity + msratemin, df) # outcome model
med_out_con = mediate(lmmed, lmout, treat = "contrast", mediator = "msratemin", control.value = 1, treat.value = 3, sims = 1000)
med_out_ecc = mediate(lmmed, lmout, treat = "scaled_eccentricity", mediator = "msratemin", control.value = 2, treat.value = 6, sims = 100)



# fit mixed models to amplitude
lm0 = lmer(amplmin ~ (1 | subject), df)
lm1 = lmer(amplmin ~ (1 | subject) + contrast, df)
lm2 = lmer(amplmin ~ (1 | subject) + eccentricity, df)
lm3 = lmer(amplmin ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(amplmin ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(amplmin ~ (1 | subject) + contrast*eccentricity, df)
lm5 = lmer(amplmin ~ (1 | subject) + contrast + eccentricity + baselineampl, df)

# fit mixed models to slip
lm0 = lmer(slipmin ~ (1 | subject), df)
lm1 = lmer(slipmin ~ (1 | subject) + contrast, df)
lm2 = lmer(slipmin ~ (1 | subject) + eccentricity, df)
lm3 = lmer(slipmin ~ (1 | subject) + contrast + eccentricity, df)
lm4 = lmer(slipmin ~ (1 | subject) + contrast*eccentricity, df)
lm5 = lmer(slipmin ~ (1 | subject) + contrast + eccentricity + baselineslip, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5, test = 'Chisq')


# same but fixed effects
lm0 = lm(msratemin ~ 1, df)
lm1 = lm(msratemin ~ contrast, df)
lm2 = lm(msratemin ~ eccentricity, df)
lm3 = lm(msratemin ~ contrast + eccentricity, df)
lm4 = lm(msratemin ~ contrast*eccentricity, df)
lm5 = lm(msratemin ~ contrast + eccentricity + baselinerate, df)
lm5b = lm(msratemin ~ contrast + scaled_eccentricity + baselinerate, df)
lm6 = lm(msratemin ~ contrast*eccentricity + baselinerate, df)
lm7 = lm(msratemin ~ baselinerate*contrast + eccentricity, df)
lm8 = lm(msratemin ~ baselinerate*contrast*eccentricity, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, test = 'Chisq')

## mediation modeling mixed
lmmedbase = lmer(baselinerate ~ (1 | subject) + contrast + eccentricity, df) # mediation model
lmmedmin = lmer(msratemin ~ (1 | subject) + contrast + eccentricity + baselinerate, df)
lmoutbase = lmer(RT ~ (1 | subject) + contrast + baselinerate + eccentricity, df) # outcome model
lmoutmin = lmer(RT ~ (1 | subject) + contrast + msratemin + eccentricity + baselinerate, df) # outcome model
med_out_base = mediate(lmmedbase, lmoutbase, treat = "contrast", mediator = "baselinerate", control.value = 1, treat.value = 3, sims = 1000)
med_out_min = mediate(lmmedmin, lmoutmin, treat = "eccentricity", mediator = "msratemin", control.value = 2, treat.value = 6, sims = 100)

## mixed mediation of baseline on relative minimum (does contrast -> baseline drive contrast -> relative min)
df2 = df
df2$msratemin = df2$msratemin - df2$baselinerate # get relative rate
lmmed = lmer(baselinerate ~ (1 | subject) + contrast + eccentricity, df2)
lmout = lmer(msratemin ~ (1 | subject) + contrast + eccentricity + baselinerate, df2)
med_out = mediate(lmmed, lmout, treat = "contrast", mediator = "baselinerate", control.value = 1, treat.value = 3, sims = 1000)


## mediation fixed

# first regress out subject random effects.
df2 = df
#tab_model(lmer(RT~ (1 + contrast + eccentricity | subject), df2))
df2$baselinerate = residuals(lmer(baselinerate ~ (1| subject), df))
df2$RT = residuals(lmer(RT ~ (1 | subject), df))
df2$msratemin = residuals(lmer(msratemin ~ (1 | subject), df))

# RESULTS: contrast - RT/baseline: all effects: contrast -> RT partially (seems weakly!) mediated by baseline change.
#          contrast - RT/min rate: direct & total only, means contrast -> RT with no mediation by min rate.
#          eccentricity - RT/baseline: direct & total only, means eccentricity -> RT with no mediation by baseline rate.
#          eccentricity: RT/min rate: direct & total only: means eccentricity -> RT with no mediation by min rate.
#          NOTE MIN RATE DOES NOT SIGNIFICANTLY MODULATE RT INDEPENDENTLY OF BASELINE RATE ANYWAY!!
#             ^^^ explained below in other mediations.
#          The last confusing observation is: minimum rate doesn't affect RT at all?
#             the weird thing is that (ONLY) if I model absolute min without including baseline, eccentricity -> RT is mediated by minimum.
#             but ecc->RT is not otherwise mediated by baseline alone or by relative rate. What's up???
#             ie putting it all together it's as if: eccentricity -> min rate -> RT. The min rate -> RT is fully mediated by baseline -> RT.
#               that's because relative min rate doesn't -> RT. But eccentricity also doesn't -> baseline rate. So how could E -> min -> RT if E doesn't -> base and min doesn't -> RT?
#          so.... the baseline rate just happens to also influence minimum rate. *As long as the baseline difference is not driven by contrast* ??
#             aka there are random fluctuations in OVERALL microsaccade rate - fluctuates baseline + minimum equally.
#             so you REALLY have to include absolute baseline in the model if you include absolute minimum, because the raw absolute minimum tracks this random fluctuation in overall microsaccade rate.
#               if you don't: WHAT HAPPENS - STILL CONFUSED ABOUT WHY IT SHOWS ECC -> MIN RATE -> RT.
#                   it almost sounds like ECC -> the random fluc -> RT. but then why doesn't it show ECC -> baseline -> RT?
#          ok, I think: eccentricity -> minimum rate. random fluctuations in overall (baseline) -> minimum rate.
#            if you don't take into account the random fluctuations, the mediation model gets *confused* and introduces false causal relationship.
#             this happens because the model is *incorrectly specified*. The correct model, with baseline included, shows no mediation effect.
lmmedbase = lm(baselinerate ~ contrast + eccentricity, df2) # mediation model
lmmedmin = lm(msratemin ~ contrast + eccentricity + baselinerate, df2)
lmoutbase = lm(RT ~ contrast + baselinerate + eccentricity, df2) # outcome model
lmoutmin = lm(RT ~ contrast + msratemin + eccentricity + baselinerate, df2) # outcome model
med_out_base = mediate(lmmedbase, lmoutbase, treat = "contrast", mediator = "baselinerate", control.value = 1, treat.value = 3, boot = TRUE, sims = 1000)
med_out_min = mediate(lmmedmin, lmoutmin, treat = "eccentricity", mediator = "msratemin", control.value = 2, treat.value = 6, boot = TRUE, sims = 1000)

senss = medsens(med_out_min, rho.by = 0.1, effect.type = "indirect", sims = 1000)
plot(senss, sens.par = 'rho', main = 'baseline', ylim = c(-0.3, 0.3))
plot(senss, sens.par = "R2", r.type = "total", sign.prod = "positive")

## mediation of baseline on relative minimum (does contrast -> baseline drive contrast -> relative min)
# RESULT: contrast - ACME & total but no ADE, means contrast -> relative min fully mediated by baseline change.
#         eccentricity - ADE & total but no ACME, means eccentricity -> relative min with no mediation by baseline.
df3 = df2
df3$msratemin = df2$msratemin - df2$baselinerate # get relative rate
lmmed = lm(baselinerate ~ contrast + eccentricity, df3)
lmout = lm(msratemin ~ contrast + eccentricity + baselinerate, df3)
med_out = mediate(lmmed, lmout, treat = "contrast", mediator = "baselinerate", control.value = 1, treat.value = 3, boot = TRUE, sims = 1000)
senss = medsens(med_out, rho.by = 0.1, effect.type = "indirect", sims = 1000)

## mediation of baseline on absolute minimum  (this isn't needed because the indiv models already show things aren't significant.)
# RESULT: contrast - ACME but no total effect, means contrast -> baseline but no contrast -> absolute minimum anyway.
#         eccentricity - ADE & total but no ACME, means eccentricity -> absolute min with no mediation by baseline.
lmmed = lm(baselinerate ~ contrast + eccentricity, df)
lmout = lm(msratemin ~ contrast + eccentricity + baselinerate, df)
med_out = mediate(lmmed, lmout, treat = "contrast", mediator = "baselinerate", control.value = 1, treat.value = 3, boot = TRUE, sims = 1000)
med_out = mediate(lmmed, lmout, treat = "eccentricity", mediator = "baselinerate", control.value = 2, treat.value = 6, boot = TRUE, sims = 1000)

## does baseline mediate the effect of minimum rate on RT?
# I'm not sure this is appropriate, obviously - minimum rate happens after baseline.
#  but perhaps it can mean baseline mediates the relative minimum rate effect.
# RESULTS: ACME but no ADE, meaning minimum rate influences RT fully via baseline.
# (if I first convert msratemin to relative rate by subtracting, you get ACME but no total. So relative min has no effect on RT.)
lmmed = lm(baselinerate ~ contrast + eccentricity + msratemin, df2)
lmout = lm(RT ~ contrast + eccentricity + baselinerate + msratemin, df2)
med_out = mediate(lmmed, lmout, treat = "msratemin", mediator = "baselinerate", control.value = -1, treat.value = 1, boot = TRUE, sims = 1000)

## does minimum mediate the effect of baseline rate on RT?
# RESULTS: baseline directly drives RT with no mediation by minimum rate (BECAUSE MIN RATE DOESN'T DRIVE RT INDEPENDENTLY OF BASELINE.)
#          same thing if I model absolute min rate or relative min rate.
lmmed = lm(msratemin ~ contrast + eccentricity + baselinerate, df)
lmout = lm(RT ~ contrast + eccentricity + baselinerate + msratemin, df)
med_out = mediate(lmmed, lmout, treat = "baselinerate", mediator = "msratemin", control.value = .06, treat.value = 1, boot = TRUE, sims = 1000)


## model baseline rate

lm0 = lmer(baselinerate ~ (1 | subject), df)
lm1 = lmer(baselinerate ~ (1 | subject) + contrast, df)
lm2 = lmer(baselinerate ~ (1 | subject) + eccentricity, df)
lm3 = lmer(baselinerate ~ (1 | subject) + contrast + eccentricity, df)
lm3b = lmer(baselinerate ~ (1 | subject) + contrast + scaled_eccentricity, df)
lm4 = lmer(baselinerate ~ (1 | subject) + contrast*eccentricity, df)
lm5 = lmer(baselinerate ~ (1 | subject) + contrast + eccentricity + msratemin, df)

lm0 = lm(baselinerate ~ 1, df)
lm1 = lm(baselinerate ~ 1 + contrast, df)
lm2 = lm(baselinerate ~ 1 + eccentricity, df)
lm3 = lm(baselinerate ~ 1 + contrast + eccentricity, df)
lm3b = lm(baselinerate ~ 1 + contrast + scaled_eccentricity, df)
lm4 = lm(baselinerate ~ 1 + contrast*eccentricity, df)
lm5 = lm(baselinerate ~ 1 + contrast + eccentricity + msratemin, df) # IT DOESN'T MAKE SENSE TO DO THIS. MS RATE MIN CANNOT PREDICT BASELINE .
anova(lm0, lm1, lm2, lm3, lm3b, lm4, lm5, test = "Chisq")

# calculate Bayes factor
exp((BIC(lm0) - BIC(lm1))/2)


# model mean rt
# eccentricity*baselinerate is significant in onlyms, but that model might have higher BIC than without the interaction
lm0 = lmer(RT ~ (1 | subject) + contrast + eccentricity + baselinerate, df) # best
lm1 = lmer(RT ~ (1 + baselinerate | subject) + contrast + eccentricity + baselinerate, df) # maybe also best (not when using onlyms)
lm1b = lmer(RT ~ (1 + baselinerate | subject) + contrast + scaled_eccentricity + baselinerate, df)
lm2 = lmer(RT ~ (1 + baselinerate | subject) + contrast*eccentricity + baselinerate, df)
lm3 = lmer(RT ~ (1 | subject) + contrast + eccentricity + baselinerate + contrast*baselinerate + eccentricity*baselinerate, df)
lm3b = lmer(RT ~ (1 + baselinerate | subject) + contrast + eccentricity + baselinerate + contrast*baselinerate + eccentricity*baselinerate, df)
lm4 = lmer(RT ~ (1 + baselinerate | subject) + contrast + eccentricity + baselinerate + msratemin, df)
lm5 = lmer(RT ~ (1 + msratemin | subject) + contrast + eccentricity + msratemin, df)
lm6 = lmer(RT ~ (1 | subject) + contrast + eccentricity + baselinerate + eccentricity*baselinerate, df)
lm7 = lmer(RT ~ (1 | subject) + contrast*eccentricity*baselinerate, df)
anova(lm0, lm1, lm1b, lm2, lm3, lm3b, lm4, lm5, test = "Chisq")

lm8 = lmer(msratemin ~ (1 | subject) + scaled_eccentricity, df)

# model mean RT - only random intercept
lm0 = lmer(RT ~ (1 | subject), df)
lm1 = lmer(RT ~ (1 | subject) + contrast, df)
lm2 = lmer(RT ~ (1 | subject) + eccentricity, df)
lm3 = lmer(RT ~ (1 | subject) + contrast + eccentricity, df)
lm4 = lmer(RT ~ (1 | subject) + contrast*eccentricity, df) # not sig in any scenario. No interaction.
lm5 = lmer(RT ~ (1 | subject) + contrast + eccentricity + baselinerate, df) # BEST.
lm5b = lmer(RT ~ (1 | subject) + contrast + scaled_eccentricity + baselinerate, df)
lm6 = lmer(RT ~ (1 | subject) + contrast*eccentricity + baselinerate, df)
lm10 = lmer(RT ~ (1 | subject) + contrast + eccentricity + baselinerate + msratemin, df)
lm11 = lmer(RT ~ (1 | subject) + contrast + baselinerate + msratemin, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5, lm5b, lm6, lm10, test = 'Chisq')





# model slope
df2 <- df
df$msratemin = df$msratemin - df$msratestart

# model selection
anova(lm0, lm1, lm2, lm3, lm4, test = "Chisq")
anova(lm0, lm1, lm2, lm3, lm4, lm5, test = "Chisq")

# regress out contrast&eccentricity&random effects, then plot residuals
lm0 = lmer(msratemin ~ (1 | subject) + contrast*eccentricity, df)
df3 = df
df3$msratemin = residuals(lm0)
plot(df3$msratemin, df3$baselinerate)

# calculate Bayes factor
exp((BIC(lm0) - BIC(lm2))/2)

# view results of a model
summary(lm3)
tab_model(lm3)

## anova and boxplots
df2 <- df %>%
  convert_as_factor(subject, contrast, eccentricity)
df2 <- df2 %>%
  reorder_levels("eccentricity", order=c("6", "4", "2"))

# summary stats
df2 %>%
  group_by(contrast) %>%
  get_summary_stats(baselinerate, type = "mean_sd")

# boxplots
bxp <- ggboxplot(
  df2, x = "eccentricity", y = "msratemin",
  color = "contrast", palette = "jco"
)
bxp

bxp <- ggboxplot(
  df2, x = "eccentricity", y = "baselinerate",
  color = "contrast", palette = "jco"
)
bxp

# boxplot of mean RTs
bxp <- ggboxplot(
  df2, x = "contrast", y = "RT",
  color = "eccentricity", palette = "rgb",
  title = "SIMULATED",#"all trials",#only trials with at least one microsaccade",#"only trials with no microsaccades or blinks",#"normalized RT by stimulus condition",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('RT (s)') + theme(text = element_text(size = 20))
#bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) + scale_colour_discrete(labels=c("low", "medium", "high")) # for eccentricity
bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) + scale_colour_discrete(labels=c("large", "medium", "small")) # for contrast
ggpar(bxp, ylim = c(3, 15))

# boxplot of mean RTs, averaged within conditions
bxp <- ggboxplot(
  df2, x = "eccentricity", y = "RT",
  color = "eccentricity", palette = "rgb",
  title = "only trials with at least one microsaccade",#"only trials with no microsaccades or blinks",#"SIMULATED",#"all trials",#"normalized RT by stimulus condition",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp = bxp + ylab('RT (s)') + theme(text = element_text(size = 20))
bxp = bxp + scale_x_discrete(labels = c("large", "medium", "small")) # for eccentricity
#bxp = bxp + scale_x_discrete(labels = c("low", "medium", "high")) # for contrast
ggpar(bxp, ylim = c(-6, 5))
ggpar(bxp, ylim = c(2, 14))



# changes
# group =eccentricity; session = contrast; index = subject, value = RT
df <- data.frame(group, session, value, index, U = interaction(session, group))
p <- ggplot(df, aes(x = U, y = value, fill = session)) + 
  scale_x_discrete(labels = rep(unique(group), each = 2))
p <- p + geom_line(aes(group = index), alpha = 0.6, colour = "black", data = df) 


# remove some columns from dataframe
keeps <- c("contrast","eccentricity","subject","msratemin","baselinerate","RT")
df2 = df2[keeps]
# check assumptions
df2 %>%
  group_by(contrast, eccentricity) %>%
  identify_outliers(RT)

df2 %>%
  group_by(contrast, eccentricity) %>%
  shapiro_test(RT)

ggqqplot(df2, "RT", ggtheme = theme_bw()) +
  facet_grid(contrast ~ eccentricity, labeller = "label_both")

# run anova
res.aov <- anova_test(
  data = df2, dv = RT, wid = subject,
  within = c(contrast, eccentricity)
)
get_anova_table(res.aov)

# pair wise post hoc
df2 %>%
  pairwise_t_test(
    msratemin ~ contrast, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )
df2 %>%
  pairwise_t_test(
    msratemin ~ eccentricity, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )


## model post-button press stuff
lm0 = lmer(msratepost ~ (1 | subject) + slippost, df)
lm0 = lmer(msratepost ~ (1 | subject) + amplpost + contrast + eccentricity, df)


### loaded everything raw: regress out stuff then test residuals.
lmrt = lmer(RT ~ (1 | subject) + contrast+eccentricity, df) # remove subj mean, contrast, eccentricity from RT
lmbr = lmer(baselinerate ~ (1 | subject) + contrast+eccentricity, df) # ^ from baseline rate
lmmr = lmer(msratemin ~ (1 | subject) + contrast+eccentricity, df) # ^ from min rate

dfres = df
dfres$RT = residuals(lmrt)
dfres$baselinerate = residuals(lmbr)
dfres$msratemin = residuals(lmmr)
lmres = lm(RT ~ baselinerate + msratemin, dfres)
plot(dfres$baselinerate, dfres$RT,
     xlab="baseline rate residuals",
     ylab = "RT residuals")
abline(lmres)

# also regress out baselinerate from msratemin then test it on RT
lmmr = lmer(msratemin ~ (1 | subject) + contrast + eccentricity + baselinerate, df)
dfres$msratemin = residuals(lmmr)
lmres = lm(RT ~ msratemin, dfres)

# also regress out all rate changes, then test for contrast&eccentricity effect.
lmrt = lmer(RT ~ (1 | subject) + baselinerate + msratemin, df)
dfres = df
dfres$RT = residuals(lmrt)
lmres = lm(RT ~ contrast + eccentricity, df)

# plot it
dfres2 <- dfres %>%
  convert_as_factor(subject, contrast, eccentricity)
bxp <- ggboxplot(
  dfres2, x = "eccentricity", y = "RT",
  color = "contrast", palette = "rgb",
  title = "regressed out microsaccade rates",
  ylab = "RT residuals",
  add = c("jitter"), add.params = list(alpha=.5)
)
bxp








######### RT
setwd("~/Documents/McGill/OneDrive - McGill University/neurospeed/UI_eyelink/bhv_data")
df <- read.csv('./trial_data.csv',sep = ",",stringsAsFactors = FALSE)
df <- df %>%
  convert_as_factor(subject, mspresent, blinkpresent, session)

# IT'S MAYBE NOT APPROPRIATE TO INCLUDE BASELINE RATE. BC MICROSACCADES ARE SO SPARSE, YOU'LL HAVE 1 OR 2 PER TRIAL. BUT A LONGER RT AUTOMATICALLY LOWERS THE HZ IF YOU ONLY HAVE 1 PER TRIAL.
# fit mixed effect models
lm0 = lmer(RT ~ (1 | subject), df)
lm1 = lmer(RT ~ (1 | subject) + contrast, df)
lm2 = lmer(RT ~ (1 | subject) + eccentricity, df)
lm3 = lmer(RT ~ (1 | subject) + contrast + eccentricity, df)
lm4 = lmer(RT ~ (1 | subject) + contrast + eccentricity + contrast:eccentricity, df)
lm5 = lmer(RT ~ (1 | subject) + contrast + eccentricity + mspresent, df)
lm5b = lmer(RT ~ (1 + mspresent | subject) + contrast + eccentricity + mspresent, df)
lm6 = lmer(RT ~ (1 | subject) + contrast*eccentricity*mspresent, df)
lm6b = lmer(RT ~ (1 + mspresent | subject) + contrast*eccentricity*mspresent, df)
lm7 = lmer(RT ~ (1 | subject) + contrast + eccentricity  + baselinenms, df)
lm8 = lmer(RT ~ (1 + baselinenms | subject) + contrast + eccentricity+ baselinenms + minnms, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5, test="Chisq")


# mediate ecc -> RT by mspresent
lmmed = lmer(msratemin ~ contrast + eccentricity + baselinerate, df)
lmout = lm(RT ~ contrast + eccentricity + baselinerate + msratemin, df)
med_out = mediate(lmmed, lmout, treat = "baselinerate", mediator = "msratemin", control.value = .06, treat.value = 1, boot = TRUE, sims = 1000)



# lm11 is good I think. accounts for all expected random subject effects.
lm10 = lmer(RT ~ (1 | subject) + contrast + eccentricity + minnms + baselinenms + baselineretinalslip, df)
lm11 = lmer(RT ~ (1 + baselineretinalslip + baselinenms | subject) + contrast + eccentricity + contrast:eccentricity + minnms + baselinenms + baselineretinalslip, df)
lm10 = lmer(RT ~ (1 | subject) + contrast + eccentricity + minnms + minretinalslip, df)
lm20 = lmer(RT ~ (1 | subject) + contrast*eccentricity*baselinenms*baselineretinalslip*minnms*minretinalslip, df)
lm21 = lmer(RT ~ (1 + contrast*eccentricity*baselinenms*baselineretinalslip*minnms*minretinalslip | subject) + contrast*eccentricity*baselinenms*baselineretinalslip*minnms*minretinalslip, df)

# model selection
anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, test = "Chisq")
anova(lm0, lm1, lm2, lm3, lm4, lm10, lm11, test = "Chisq")

# fixed effects (if already removed subject means from RT and baselinenms)
lm0 = lm(RT ~ 1, df)
lm1 = lm(RT ~ contrast, df)
lm2 = lm(RT ~ eccentricity, df)
lm3 = lm(RT ~ contrast + eccentricity, df)
lm4 = lm(RT ~ contrast + eccentricity + contrast:eccentricity, df)
lm5 = lm(RT ~ contrast + eccentricity + contrast:eccentricity + baselinenms, df)
anova(lm0, lm1, lm2, lm3, lm4, lm5, test = "Chisq")


# model baseline dynamics from stimulus
lms0 = lmer(baselinenms ~ (1 | subject), df)
lms1 = lmer(baselinenms ~ (1 | subject) + contrast, df)
lms2 = lmer(baselinenms ~ (1 | subject) + eccentricity, df)
lms3 = lmer(baselinenms ~ (1 | subject) + contrast + eccentricity, df)
lms4 = lmer(baselinenms ~ (1 | subject) + contrast + eccentricity + contrast:eccentricity, df)

lms5 = lmer(baselinenms ~ (1 + contrast + eccentricity | subject) + contrast + eccentricity, df)
lms6 = lmer(baselinenms ~ (1 + contrast + eccentricity | subject) + contrast*eccentricity, df)

anova(lms0, lms1, lms2, lms3, lms4, lms5, lms6, test = "Chisq")

# calculate Bayes factor
exp((BIC(lm11) - BIC(lm20))/2)


## anova
df2 <- df %>%
  convert_as_factor(subject, contrast, eccentricity)

# boxplots
bxp <- ggboxplot(
  df2, x = "eccentricity", y = "immobduration",
  color = "contrast", palette = "jco"
)
bxp



# mixed effect models of retinal slip
lm0 = lmer(retinalslip ~ (1 | subject), df)
lm1 = lmer(retinalslip ~ (1 | subject) + contrast, df)
lm2 = lmer(retinalslip ~ (1 | subject) + eccentricity, df)
lm3 = lmer(retinalslip ~ (1 | subject) + contrast + eccentricity, df)
lm4 = lmer(retinalslip ~ (1 | subject) + contrast + eccentricity + contrast:eccentricity, df)


# view results of a model
summary(lm3)
tab_model(lm3)

## anova
df2 <- df %>%
  convert_as_factor(subject, contrast, eccentricity)

# summary stats
df2 %>%
  group_by(contrast, eccentricity) %>%
  get_summary_stats(retinalslip, type = "mean_sd")

# boxplots
bxp <- ggboxplot(
  df2, x = "eccentricity", y = "retinalslip",
  color = "contrast", palette = "jco"
)
bxp

# check assumptions
df2 %>%
  group_by(contrast, eccentricity) %>%
  identify_outliers(retinalslip)

df2 %>%
  group_by(contrast, eccentricity) %>%
  shapiro_test(retinalslip)

ggqqplot(df2, "RT", ggtheme = theme_bw()) +
  facet_grid(contrast ~ eccentricity, labeller = "label_both")

# run anova
res.aov <- anova_test(
  data = df2, dv = retinalslip, wid = subject,
  within = c(contrast, eccentricity)
)
get_anova_table(res.aov)

# pair wise post hoc
df2 %>%
  pairwise_t_test(
    RT ~ contrast, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )
df2 %>%
  pairwise_t_test(
    RT ~ eccentricity, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )


# model minimum nms by trial
lm0 = lmer(minnms ~ (1 | subject), df)
lm1 = lmer(minnms ~ (1 | subject) + contrast, df)
lm2 = lmer(minnms ~ (1 | subject) + eccentricity, df)
lm3 = lmer(minnms ~ (1 | subject) + contrast + eccentricity, df)
lm4 = lmer(minnms ~ (1 | subject) + contrast*eccentricity, df)
lm5 = lmer(minnms ~ (1 | subject) + contrast + eccentricity + baselinenms, df)
lm5a = lmer(minnms ~ (1 + baselinenms | subject) + contrast + eccentricity + baselinenms, df)
lm6 = lmer(minnms ~ (1 | subject) + baselinenms*contrast + eccentricity, df)
lm7 = lmer(minnms ~ (1 | subject) + baselinenms*contrast*eccentricity, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5, lm5a, lm6, lm7, test = "Chisq")

# same but fixed effects
lm0 = lm(minnms ~ 1, df)
lm1 = lm(minnms ~ contrast, df)
lm2 = lm(minnms ~ eccentricity, df)
lm3 = lm(minnms ~ contrast + eccentricity, df)
lm4 = lm(minnms ~ contrast*eccentricity, df)
lm5 = lm(minnms ~ contrast + eccentricity + baselinenms, df)
lm6 = lm(minnms ~ baselinenms*contrast + eccentricity, df)
lm7 = lm(minnms ~ baselinenms*contrast*eccentricity, df)

anova(lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, test = "Chisq")


## look into "immobilization duration"
lm0 = lmer(immobduration ~ (1 | subject), df)
lm1 = lmer(immobduration ~ (1 | subject) + contrast, df)
lm2 = lmer(immobduration ~ (1 | subject) + eccentricity, df)
lm3 = lmer(immobduration ~ (1 | subject) + contrast + eccentricity, df)
lm3a = lmer(immobduration ~ (1 + contrast | subject) + contrast + eccentricity, df) # BEST FIT.
lm4 = lmer(immobduration ~ (1 | subject) + contrast*eccentricity, df)
lm4b = lmer(immobduration ~ (1 + contrast | subject) + contrast*eccentricity, df)
lm5 = lmer(immobduration ~ (1 + contrast + eccentricity | subject), df)
lm5b = lmer(immobduration ~ (1 + contrast | subject) + eccentricity, df)

anova(lm0, lm1, lm2, lm3, lm3a, lm4, lm4b, lm5, lm5b, test = 'Chisq')

exp((BIC(lm5b) - BIC(lm5))/2)


# predict rt from immob
lm0 = lmer(RT ~ (1 + immobduration | subject) + immobduration, df)
lm1 = lmer(RT ~ (1 + contrast + eccentricity + immobduration | subject) + immobduration + contrast + eccentricity, df)
lm2 = lmer(RT ~ (1 + contrast + eccentricity + immobduration | subject) + immobduration + contrast*eccentricity, df)

# do mixed effects
# first remove subject effects
dfres = df
dfres$RT = residuals(lmer(RT ~ (1 | subject), df))
dfres$immobduration = residuals(lmer(immobduration ~ (1 | subject), df))

lm0 = lm(immobduration ~ 1, dfres)
lm1 = lm(immobduration ~ contrast, dfres)
lm2 = lm(immobduration ~ eccentricity, dfres)
lm3 = lm(immobduration ~ contrast + eccentricity, dfres)
lm4 = lm(immobduration ~ contrast*eccentricity, dfres)

anova(lm0, lm1, lm2, lm3, lm4, test = 'Chisq')

# mediation
lmmed = lm(immobduration ~ contrast + eccentricity, dfres) # mediation model
lmout = lm(RT ~ contrast + eccentricity + immobduration, dfres) # outcome model
med_out_ecc = mediate(lmmed, lmout, treat = "eccentricity", mediator = "immobduration", control.value = 2, treat.value = 6, boot = TRUE, sims = 1000)
med_out_con = mediate(lmmed, lmout, treat = "contrast", mediator = "immobduration", control.value = 1, treat.value = 3, boot = TRUE, sims = 1000)

lmmed = lm(RT ~ contrast + eccentricity, dfres) # mediation model
lmout = lm(immobduration ~ contrast + eccentricity + RT, dfres) # outcome model
med_out_ecc = mediate(lmmed, lmout, treat = "eccentricity", mediator = "RT", control.value = 2, treat.value = 6, boot = TRUE, sims = 1000)
med_out_con = mediate(lmmed, lmout, treat = "contrast", mediator = "RT", control.value = 1, treat.value = 3, boot = TRUE, sims = 1000)


# predict rt from immob
lm0 = lm(RT ~ immobduration + contrast + eccentricity, dfres)


## plot
df2 <- dfres %>%
  convert_as_factor(subject, contrast, eccentricity)

# summary stats
df2 %>%
  group_by(contrast, eccentricity) %>%
  get_summary_stats(immobduration, type = "mean_sd")

# boxplots
bxp <- ggboxplot(
  df2, x = "contrast", y = "immobduration",
  #color = "contrast", palette = "jco",
  #outlier.shape = NA # hide outliers
)
bxp
# don't look at outliers
ylim1 = boxplot.stats(df2$immobduration)$stats[c(1,5)]
bxp1 = bxp + coord_cartesian(ylim = ylim1*1.05)
bxp1



