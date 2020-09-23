# ---
# title: "Mining Sets"
# author: "Casey Cazer"
# Last Updated: "Jan 6, 2020"
# ---

#mining sets

#create transaction databases
varNames <- c("SA2008", "SA2009", "SA2010", "SA2011", "SA2012", "SA2013", "SA2014", "SA2015", "SA2016", "SA2017", "SA2018")
for (i in 1:length(varNames)){ #for each year
  x <- as(get(dbNames[i]), "transactions") #make binary data into transactions
  label <- paste('SA',as.character(paste(2007+i)), sep="")
  assign(label, x) #save
  rm(x)
  rm(label)
}


#optimistic transaction databases for calculating expected support, CSR, lift
#optiministic transaction databases replace missing data (NA) with "true" (resistant)
oTransNames <- c("SA2008_o", "SA2009_o", "SA2010_o", "SA2011_o", "SA2012_o", "SA2013_o", "SA2014_o", "SA2015_o", "SA2016_o", "SA2017_o", "SA2018_o")
for (i in 1:length(oTransNames)){ #for each year
  x <- SA.db[which(SA.db$Study.Year==as.numeric(paste(2007+i))), itemLabels(get(varNames[i]))] #select only AM columns that are in the actual transaction database (e.g have at least one isolate resistant), otherwise columns with all NA will be included and then turned into all resistant
  x[is.na(x)] <- as.logical("TRUE") #replace NA with resistant 
  x <- as(x, "transactions") #make binary data into transactions
  label <- paste('SA',as.character(paste(2007+i)), "_o", sep="")
  assign(label, x) #save
  rm(x)
  rm(label)
}


#mine sets
setNames=c("SA2008_sets", "SA2009_sets", "SA2010_sets", "SA2011_sets", "SA2012_sets", "SA2013_sets", "SA2014_sets", "SA2015_sets", "SA2016_sets", "SA2017_sets", "SA2018_sets")

for (i in seq_along(varNames)){ #for each year
  #get relevant transaction, optimistic transaction, and binary databases
  data <- get(varNames[i]) 
  odata <- get(oTransNames[i])
  db <- get(dbNames[i])
  minsup <- 1/length(data) #minimum support = 1 isolate resistant
  #mine sets
  sets <- apriori(data, parameter=list(support=minsup,  maxlen=length(AM_col_db), minlen=2, target="frequent itemsets"))
  
  #get quality measures
  itemset_list <- LIST(items(sets), decode = FALSE) #required for expectedQM
  CrossSupRatio <- interestMeasure(sets, "crossSupportRatio", data, reuse=TRUE)
  lift <- interestMeasure(sets, "lift", data, reuse=TRUE)
  sets@quality$CrossSupRatio <- CrossSupRatio #save to set database
  sets@quality$lift <- lift
  eQM <- expectedQM(sets, data, odata, db, itemset_list) #calculate expected quality measures (QM)
  sets@quality$eSup <- eQM$eSup
  sets@quality$eCSR <- eQM$eCSR
  sets@quality$eLift <- eQM$eLift
  sets@quality$cLift <- eQM$cLift
  
  #save
  assign(setNames[i], sets)
  
  rm(data, odata, db, minsup, sets, itemset_list, CrossSupRatio, lift, eQM)
}



###################################
#MRSA v MSSA
#create transactions
mvarNames=c("MRSA", "MSSA")
MRSA=as(MRSA_db, "transactions")
MSSA=as(MSSA_db, "transactions")

msetNames=c("MRSA_sets", "MSSA_sets")

#optimistic transaction databases for calculating expected support, CSR, lift
omTransNames <- c("MRSA_o", "MSSA_o")
  MRSA_o <- SA.db[which(SA.db$Oxacillin==TRUE), itemLabels(MRSA)] #select only AM columns that are in the actual transaction database (otherwise columns with all NA will be included)
  MRSA_o[is.na(MRSA_o)] <- as.logical("TRUE") #replace NA with resistant
  MRSA_o <- as(MRSA_o, "transactions") #binary data to transactions
  
  MSSA_o <- SA.db[which(SA.db$Oxacillin==FALSE), itemLabels(MSSA)] #select only AM columns that are in the actual transaction database (otherwise columns with all NA will be included)
  MSSA_o[is.na(MSSA_o)] <- as.logical("TRUE") #replace NA with resistant
  MSSA_o <- as(MSSA_o, "transactions") #binary data to transactions



for (i in seq_along(mvarNames)){ #for each category
  #get relevant databases
  data=get(mvarNames[i])
  odata <- get(omTransNames[i])
  db <- get(mdbNames[i])
  
  minsup=1/length(data) #minimum support = 1 isolate resistant
  sets=apriori(data, parameter=list(support=minsup,  maxlen=length(m.AM_col_db), minlen=2, target="frequent itemsets"))
 
  #get quality measures (QM)
  itemset_list <- LIST(items(sets), decode = FALSE)
  CrossSupRatio <- interestMeasure(sets, "crossSupportRatio", data, reuse=TRUE)
  lift <- interestMeasure(sets, "lift", data, reuse=TRUE)
  sets@quality$CrossSupRatio <- CrossSupRatio
  sets@quality$lift <- lift
  eQM <- expectedQM(sets, data, odata, db, itemset_list)
  sets@quality$eSup <- eQM$eSup
  sets@quality$eCSR <- eQM$eCSR
  sets@quality$eLift <- eQM$eLift
  sets@quality$cLift <- eQM$cLift
  
  #save
  assign(msetNames[i], sets)
  
  rm(data, odata, db, minsup, sets, itemset_list, CrossSupRatio, lift, eQM)
}
  

 #Infection type
  type.varNames <- c("SA.blood", "SA.inabd", "SA.pneum", "SA.SST")
  for (i in 1:length(type.varNames)){ #build transaction databases
    x <- as(get(type.dbNames[i]), "transactions")
    label <- paste('SA.',type.abbrev[i], sep="")
    assign(label, x)
    rm(x)
    rm(label)
  }
  
  
  #optimistic transaction databases for calculating expected support, CSR, lift
  type.oTransNames <- c("SA.blood_o", "SA.inabd_o", "SA.pneum_o", "SA.SST_o")
  for (i in 1:length(type.oTransNames)){
    x <- SA.db[which(SA.db$Infection.Type==type[i]), itemLabels(get(type.varNames[i]))] #select only AM columns that are in the actual transaction database (otherwise columns with all NA will be included)
    x[is.na(x)] <- as.logical("TRUE")
    x <- as(x, "transactions")
    label <- paste('SA.', type.abbrev[i], "_o", sep="")
    assign(label, x)
    rm(x)
    rm(label)
  }
  
  
  #mine sets
  type.setNames=c("SA.blood_sets", "SA.inabd_sets", "SA.pneum_sets", "SA.SST_sets")
  
  for (i in seq_along(type.varNames)){ #for each category
    #get relevant databases
    data <- get(type.varNames[i])
    odata <- get(type.oTransNames[i])
    db <- get(type.dbNames[i])
    minsup <- 1/length(data) #minimum support = 1 isolate resistant
    
    #mine sets
    sets <- apriori(data, parameter=list(support=minsup,  maxlen=length(AM_col_db), minlen=2, target="frequent itemsets"))
    
    #get quality measures
    itemset_list <- LIST(items(sets), decode = FALSE)
    CrossSupRatio <- interestMeasure(sets, "crossSupportRatio", data, reuse=TRUE)
    lift <- interestMeasure(sets, "lift", data, reuse=TRUE)
    sets@quality$CrossSupRatio <- CrossSupRatio
    sets@quality$lift <- lift
    eQM <- expectedQM(sets, data, odata, db, itemset_list)
    sets@quality$eSup <- eQM$eSup
    sets@quality$eCSR <- eQM$eCSR
    sets@quality$eLift <- eQM$eLift
    sets@quality$cLift <- eQM$cLift
    
    #save
    assign(type.setNames[i], sets)
    
    rm(data, odata, db, minsup, sets, itemset_list, CrossSupRatio, lift, eQM)
  }
  