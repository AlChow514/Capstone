# ---
# title: "filtering and summarizing sets"
# author: "Casey Cazer"
# Last updated: "January 6, 2020"
# ---

#filter sets based on eCSR P-value

best_setNames=c("SA2008_bestsets", "SA2009_bestsets", "SA2010_bestsets", "SA2011_bestsets", "SA2012_bestsets", "SA2013_bestsets", "SA2014_bestsets", "SA2015_bestsets", "SA2016_bestsets", "SA2017_bestsets", "SA2018_bestsets")
for (i in seq_along(setNames)){ #for each set db
  sets <- get(setNames[i]) #get sets
  bestsets <- subset(sets, subset=Pval.eCSR<=0.05) #select only those sets with eCSR P-value <= 0.05
  assign(best_setNames[i], bestsets) #save as 'bestsets'
  rm(sets, bestsets)
}

#present data as table of sets and quality measures
#summarize all sets found in one or more bestsets collections
all.sets <- all_sets(best_setNames, varNames)
#all.sets returns several summarizing dataframes as a list
##first dataframe 'all.sets' has a row for each set found in each db (e.g. year) and columns for quality measures. number of rows = sum of number of sets across all db (e.g. year)
##other dataframes (e.g. 'all.sets.sup') have one row for each unique set and columns for each db, giving the specific quality measure (e.g. sup) for the set in that db (NA if the set is not in the bestsets for that db)

#get AM classifications for each set
all.sets_classes <- lapply(all.sets, itemset_to_class_set, AM_class) #appends columns of AM classes, number of classes, and class codes to each dataframe in the list
names(all.sets_classes) <- names(all.sets)

################################
#MRSA v MSSA

mbest_setNames=c("MRSA_bestsets", "MSSA_bestsets")
for (i in seq_along(msetNames)){
  sets=get(msetNames[i])
  bestsets <- subset(sets, subset=Pval.eCSR<=0.05)
  assign(mbest_setNames[i], bestsets)
  rm(sets, bestsets)
}

m.all.sets <- all_sets(mbest_setNames, mvarNames)

m.all.sets_classes <- lapply(m.all.sets, itemset_to_class_set, AM_class)
names(m.all.sets_classes) <- names(m.all.sets)

###################
#Infection type

type.best_setNames=c("SA.blood_bestsets", "SA.inabd_bestsets", "SA.pneum_bestsets", "SA.SST_bestsets")
for (i in seq_along(type.setNames)){
  sets <- get(type.setNames[i])
  bestsets <- subset(sets, subset=Pval.eCSR<=0.05)
  assign(type.best_setNames[i], bestsets)
  rm(sets, bestsets)
}

type.all.sets <- all_sets(type.best_setNames, type.varNames)

type.all.sets_classes <- lapply(type.all.sets, itemset_to_class_set, AM_class)
names(type.all.sets_classes) <- names(type.all.sets)

##############################
#save all together, saving the first dataframe with all QM in columns and each set/category in a row
all.sets_combined <- rbind(all.sets_classes[[1]], m.all.sets_classes[[1]], type.all.sets_classes[[1]])
