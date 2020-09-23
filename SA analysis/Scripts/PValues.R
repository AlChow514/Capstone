# ---
# title: "P value calculations"
# author: "Casey Cazer"
# Last updated: "Jan 6, 2020"
# output: html_document
# ---
#use PValue.sets to calculate null distributions of QM for each db
#year db
rand_dbNames=c("rand_SA2008", "rand_SA2009", "rand_SA2010", "rand_SA2011", "rand_SA2012", "rand_SA2013", "rand_SA2014", "rand_SA2015", "rand_SA2016", "rand_SA2017", "rand_SA2018")

rand_setNames=c("rand_SA2008_sets", "rand_SA2009_sets", "rand_SA2010_sets", "rand_SA2011_sets", "rand_SA2012_sets", "rand_SA2013_sets", "rand_SA2014_sets", "rand_SA2015_sets", "rand_SA2016_sets", "rand_SA2017_sets", "rand_SA2018_sets")

Dist.out <- PValue.sets(dbNames, varNames, setNames, rand_dbNames, rand_setNames, 2, length(AM_col_db))
#returns the distribution of quality measures (CSR, eCSR, Lift, eLift, cLift) under null hypothesis of no association between resistances
#the distributions are given as the average quality measure at each percentile (0 to 0.95 by 0.01 and 0.95 to 1 by 0.001)
#Dist.out contains a list of dataframes, one data frame for each QM. In the dataframe, columns refer to the null distribution for each specific db (e.g. year)

#calculate P-value for actual sets based on null distribution. The P value is the percent of observations under null hypothesis that have a QM >= the acutal QM (this is a one-sided P value)
#First find the percentile of each actual QM value in the relevant null distribution. 
##detect_index searches from 100th percentile down, returning the row.number that is the first with a QM value less than the actual QM value. 
###The true percentile is slightly greater than the returned value (lies between the returned value and the next greatest percentile). Since the percentiles are unequally spaced (greater resolution at 0.95 to 1), this uses the row name to report the percentile

#must recode 0 as 1 (mapvvalues) to correctly reference row.names. position=0 indicates that the QM was smaller than the smallest QM in null distribution
#P-value is 1 - percentile. It is an upper bound/conservative estimate; the true P-value is slightly less than reported since the true percentile is slightly greater than reported
for (i in seq_along(setNames)){ #for each db of sets
  sets <- get(setNames[i]) #get sets
  
  #find percentile of each actual QM
  perc.ecsr <- as.numeric(row.names(Dist.out[["set.ecsr_pct_avg"]])
                          [adply(sets@quality$eCSR, .margins=1, #margins=1 splits by rows--so for each set eCSR value x
                                 function(x) detect_index(Dist.out[["set.ecsr_pct_avg"]][,i], #search within the null distribution for the QM and the relevant db i
                                      function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% #find where this function is true--the avg null QM value is <= set eCSR value
                                          mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.csr <- as.numeric(row.names(Dist.out[["set.csr_pct_avg"]])
                         [adply(sets@quality$CrossSupRatio, .margins=1, 
                                function(x) detect_index(Dist.out[["set.csr_pct_avg"]][,i], 
                                      function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                                          mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.eLift <- as.numeric(row.names(Dist.out[["set.elift_pct_avg"]])
                           [adply(sets@quality$eLift, .margins=1, 
                                  function(x) detect_index(Dist.out[["set.elift_pct_avg"]][,i], 
                                      function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                                          mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.lift <- as.numeric(row.names(Dist.out[["set.lift_pct_avg"]])
                          [adply(sets@quality$lift, .margins=1, 
                                 function(x) detect_index(Dist.out[["set.lift_pct_avg"]][,i], 
                                      function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                                          mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.cLift <- as.numeric(row.names(Dist.out[["set.clift_pct_avg"]])
                          [adply(sets@quality$cLift, .margins=1, 
                                 function(x) detect_index(Dist.out[["set.clift_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                              mapvalues(0,1, warn_missing=FALSE)]) 
  
  #save as P-values
  sets@quality$Pval.eCSR <- 1-perc.ecsr
  sets@quality$Pval.CSR <- 1-perc.csr
  sets@quality$Pval.eLift <- 1-perc.eLift
  sets@quality$Pval.lift <- 1-perc.lift
  sets@quality$Pval.cLift <- 1-perc.cLift
  
  assign(setNames[i], sets)
  rm(sets, perc.ecsr, perc.csr, perc.eLift, perc.lift, perc.cLift)
}

#################################
#MRSA v MSSA, same procedure as above
mrand_dbNames=c("rand_MRSA", "rand_MSSA")
mrand_setNames=c("rand_MRSA_sets", "rand_MSSA_sets")

mDist.out <- PValue.sets(mdbNames, mvarNames, msetNames, mrand_dbNames, mrand_setNames, 2, length(m.AM_col_db))

#calculate P-values
for (i in seq_along(msetNames)){ #for each db of sets
  sets <- get(msetNames[i])
  
  #find percentile of each actual QM
  perc.ecsr <- as.numeric(row.names(mDist.out[["set.ecsr_pct_avg"]])
                          [adply(sets@quality$eCSR, .margins=1, 
                                 function(x) detect_index(mDist.out[["set.ecsr_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>%
                              mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.csr <- as.numeric(row.names(mDist.out[["set.csr_pct_avg"]])
                         [adply(sets@quality$CrossSupRatio, .margins=1, 
                                function(x) detect_index(mDist.out[["set.csr_pct_avg"]][,i], 
                                                         function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                             mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.eLift <- as.numeric(row.names(mDist.out[["set.elift_pct_avg"]])
                           [adply(sets@quality$eLift, .margins=1, 
                                  function(x) detect_index(mDist.out[["set.elift_pct_avg"]][,i], 
                                                           function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                               mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.lift <- as.numeric(row.names(mDist.out[["set.lift_pct_avg"]])
                          [adply(sets@quality$lift, .margins=1, 
                                 function(x) detect_index(mDist.out[["set.lift_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                              mapvalues(0,1, warn_missing=FALSE)])
  
  perc.cLift <- as.numeric(row.names(mDist.out[["set.clift_pct_avg"]])
                          [adply(sets@quality$cLift, .margins=1, 
                                 function(x) detect_index(mDist.out[["set.clift_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                              mapvalues(0,1, warn_missing=FALSE)])
  
  #save as P-values
  sets@quality$Pval.eCSR <- 1-perc.ecsr
  sets@quality$Pval.CSR <- 1-perc.csr
  sets@quality$Pval.eLift <- 1-perc.eLift
  sets@quality$Pval.lift <- 1-perc.lift
  sets@quality$Pval.cLift <- 1-perc.cLift
  
  assign(msetNames[i], sets)
  rm(sets, perc.ecsr, perc.csr, perc.eLift, perc.lift, perc.cLift)
}

###########################
#Infection Type
type.rand_dbNames=c("rand_SA.blood", "rand_SA.inabd", "rand_SA.pneum", "rand_SA.SST")

type.rand_setNames=c("rand_SA.blood_sets", "rand_SA.inabd_sets", "rand_SA.pneum_sets", "rand_SA.SST_sets")

type.Dist.out <- PValue.sets(type.dbNames, type.varNames, type.setNames, type.rand_dbNames, type.rand_setNames, 2, length(AM_col_db))

#calculate P-values
for (i in seq_along(type.setNames)){ #for each db
  sets <- get(type.setNames[i])
  
  #find percentile of each actual QM
  perc.ecsr <- as.numeric(row.names(type.Dist.out[["set.ecsr_pct_avg"]])
                          [adply(sets@quality$eCSR, .margins=1, 
                                 function(x) detect_index(type.Dist.out[["set.ecsr_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>%
                              mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.csr <- as.numeric(row.names(type.Dist.out[["set.csr_pct_avg"]])
                         [adply(sets@quality$CrossSupRatio, .margins=1, 
                                function(x) detect_index(type.Dist.out[["set.csr_pct_avg"]][,i], 
                                                         function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                             mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.eLift <- as.numeric(row.names(type.Dist.out[["set.elift_pct_avg"]])
                           [adply(sets@quality$eLift, .margins=1, 
                                  function(x) detect_index(type.Dist.out[["set.elift_pct_avg"]][,i], 
                                                           function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                               mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.lift <- as.numeric(row.names(type.Dist.out[["set.lift_pct_avg"]])
                          [adply(sets@quality$lift, .margins=1, 
                                 function(x) detect_index(type.Dist.out[["set.lift_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                              mapvalues(0,1, warn_missing=FALSE)])
  
  perc.cLift <- as.numeric(row.names(type.Dist.out[["set.clift_pct_avg"]])
                          [adply(sets@quality$cLift, .margins=1, 
                                 function(x) detect_index(type.Dist.out[["set.clift_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                              mapvalues(0,1, warn_missing=FALSE)])
  
  #save as P-values
  sets@quality$Pval.eCSR <- 1-perc.ecsr
  sets@quality$Pval.CSR <- 1-perc.csr
  sets@quality$Pval.eLift <- 1-perc.eLift
  sets@quality$Pval.lift <- 1-perc.lift
  sets@quality$Pval.cLift <- 1-perc.cLift
  
  assign(type.setNames[i], sets)
  rm(sets, perc.ecsr, perc.csr, perc.eLift, perc.lift, perc.cLift)
}
