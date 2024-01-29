# load packages
library(lme4)
library(lmeresampler)
library(rstatix)
library(foreach)
library(doParallel)
numCores <- 20
numReps <- 10000
numRepsPerCore = numReps/numCores
cl <- makeCluster(numCores, type = "PSOCK")
doParallel::registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths()) # I think register loaded packages in each cluster

## ---- load individual trials ---- ##
df <- read.csv('./trial_data_withblinks.csv',sep = ",",stringsAsFactors = FALSE)
df$FT = df$RT
df$baselinemsrate = df$baselinenms / df$FT
df <- df %>%
  convert_as_factor(subject, mspresent, earlymspresent, latemspresent, blinkpresent, session)
dforig = df

# rescale and zero-mean all continuous predictors (but not FT or immobduration)
numcols <- sapply(df, is.numeric)
dfs <- df
dfs[,numcols] <- scale(dfs[,numcols])
#dfs$FT = dforig$FT / 10
dfnocenter <- df
dfnocenter[,numcols] <- scale(dfnocenter[,numcols], center = FALSE, scale = apply(dfnocenter[,numcols], 2, sd, na.rm = TRUE))
dfs$FT = dfnocenter$FT
dfs$immobduration = dfnocenter$immobduration
dfs$baselinemsrate = dfnocenter$baselinemsrate
logdfs = dfs
logdfs$FT = log(dfs$FT)
logdfs$immobduration = log(dfs$immobduration)
logdfs$baselinemsrate = log(dfs$baselinemsrate)


# predict baseline RATE, only from trials with at least 1 microsaccade.
logdfswithms = logdfs
logdfswithms = logdfswithms[(dforig$baselinenms > 0),]
#lm6log = lmer(baselinemsrate ~ (contrast + eccentricity + trial | subject) + contrast + eccentricity + trial, logdfswithms, control=lmerControl("bobyqa"))
#boo6 <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(lm6log, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#}
#print(confint(boo6, level = 0.9875, type = "perc"))
#save(boo6, file = "baselinerate_bootstrap.Rdata")


# immobilization duration
#lm5log = lmer(immobduration ~ (contrast + eccentricity + immobmsamplitude + trial | subject) + contrast + eccentricity + immobmsamplitude + trial, logdfs, control=lmerControl("bobyqa"))
#boo5 <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(lm5log, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#}
#print(confint(boo5, level = 0.9875, type = "perc"))
#save(boo5, file = "immob_bootstrap.Rdata}")

## 3 step model creation.
# Step 1: model with no eye movements.
#loglmnoeyes = lmer(FT ~ (trial + contrast + eccentricity | subject) + contrast*eccentricity + contrast*trial + eccentricity*trial, logdfs, control=lmerControl(optimizer="bobyqa"))
#boonoeyes <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(loglmnoeyes, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#}
#print(confint(boonoeyes, level = 0.9917, type = "perc"))
#save(boonoeyes, file = "noeyes_bootstrap.Rdata")

# Step 2: use selected parameters from 1, and add eye movement main effects.
#loglmmain = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity  + baselineretinalslip | subject) + mspresent + nblinks + baselineretinalslip + trial*contrast + eccentricity, logdfs, control=lmerControl(optimizer="bobyqa"))
#boomain <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(loglmmain, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#print(confint(boomain, level = 0.9929, type = "perc"))
#save(boomain, file = "main_bootstrap.Rdata")

# Step 3: use selected parameters from 2, and add interactions bw C/E and eye movements. Bonferroni-corrected interaction selection.
#loglmselected = lmer(FT ~ (mspresent + nblinks + trial + contrast + eccentricity | subject) + mspresent*contrast + mspresent*eccentricity + nblinks*contrast + nblinks*eccentricity + trial*contrast, logdfs, control=lmerControl(optimizer="bobyqa"))
#booselected <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(loglmselected, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#}
#print(confint(booselected, level = 0.9955, type = "perc"))
#save(booselected, file = "selected_bootstrap.Rdata")


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
#loglmnomsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfsnomssubset, control=lmerControl("bobyqa"))
#boonomsnoblinks <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(loglmnomsnoblinks, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#}
#print(confint(boonomsnoblinks, level = 0.95, type = "perc"))
#save(boonomsnoblinks, file = "nomssubset_bootstraps.Rdata")

# Next: trial subset with at least 1 microsaccade, no blinks.
logdfswithmssubset = dfswithmssubset
logdfswithmssubset$FT = log(dfswithmssubset$FT)
#loglmwithmsnoblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfswithmssubset, control=lmerControl("bobyqa"))
#boowithmsnoblinks <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(loglmwithmsnoblinks, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#}
#print(confint(boowithmsnoblinks, level = 0.95, type = "perc"))
#save(boowithmsnoblinks, file = "withmssubset_bootstraps.Rdata")

# Next: trials with at least one blink, no microsaccades (only 118 trials)
# no interactions done because there are so, so few trials. And just C + E parameters.
dfsnomswithblinks = dfs[(dforig$baselinenms == 0) & (dforig$nblinks > 0),]
dforignomswithblinks = dforig[(dforig$baselinenms == 0) & (dforig$nblinks > 0),]
logdfsnomswithblinks = dfsnomswithblinks
logdfsnomswithblinks$FT = log(dfsnomswithblinks$FT)
#loglmnomswithblinks = lmer(FT ~ (1 | subject) + contrast + eccentricity + trial + baselineretinalslip, logdfsnomswithblinks, control=lmerControl("bobyqa"))
#boonomswithblinks <- foreach(B = rep(numRepsPerCore, numCores), .combine = combine_lmeresamp, .packages = c("lmeresampler", "lme4")) %dopar% {
#    bootstrap(loglmnomswithblinks, .f = fixef, type = "case", B = B, resample = c(TRUE, FALSE))
#}
#print(confint(boonomswithblinks, level = 0.95, type = "perc"))
#save(boonomswithblinks, file = "nomswithblinks_bootstraps.Rdata")



stopCluster(cl)

