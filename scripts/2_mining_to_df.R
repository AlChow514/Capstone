# Load Libraries
library(tidyverse)
library(arules)

# find indexes of antimicrobial column to use in mining (exclude those missing due to no bp); ecv and clsi data is the same
am_col <- eu_mdr_df %>% 
  select(1:15) #collection number, infection type, drugs. exclude AM class columns and mdr, mdro
am_col <- match(names(am_col)[sapply(am_col, is.logical)], names(am_col)) #include only logical T/F columns (exclude collection number and inf type), get column number

# eu sets
eu_db_sets <- c("eu_all_db", "eu_bld_db", "eu_inab_db", "eu_pneu_db", "eu_skin_db")

# create eu transaction databases
for (i in seq_along(eu_db_sets)) {
  x <- as(get(eu_db_sets[i]), "transactions")
  label <- paste(as.character(paste(eu_db_sets[i], "_trans", sep = "")))
  assign(label, x)
  rm(label, x)
}

# Mine sets
eu_trans_names <- c("eu_all_db_trans", "eu_bld_db_trans", "eu_inab_db_trans", "eu_pneu_db_trans", "eu_skin_db_trans")
eu_set_names <- c("eu_all_set", "eu_bld_set", "eu_inab_set", "eu_pneu_set", "eu_skin_set")
for (i in seq_along(eu_trans_names)) {
  data <- get(eu_trans_names[i])
  # minimum support
  min_sup <- 1 / length(data)
  sets <- apriori(data,
                  parameter = list(support = min_sup,
                                   maxlen = length(am_col),
                                   minlen = 2,
                                   target = "frequent itemsets"))
  # Quality Measures
  CrossSupRatio <- interestMeasure(sets,
                                   "crossSupportRatio",
                                   data,
                                   reuse = TRUE)
  
  sets@quality$CrossSupRatio <- CrossSupRatio
  
  lift <- interestMeasure(sets,
                          "lift",
                          data,
                          reuse = TRUE)
  
  sets@quality$lift <- lift
  
  # save
  assign(eu_set_names[i], sets)
  rm(data, min_sup, sets, CrossSupRatio, lift)
}

# clsi sets
clsi_db_sets <- c("clsi_all_db", "clsi_bld_db", "clsi_inab_db", "clsi_pneu_db", "clsi_skin_db")

# create clsi transaction databases
for (i in seq_along(clsi_db_sets)) {
  x <- as(get(clsi_db_sets[i]), "transactions")
  label <- paste(as.character(paste(clsi_db_sets[i], "_trans", sep = "")))
  assign(label, x)
  rm(label, x)
}

# Mine sets
clsi_trans_names <- c("clsi_all_db_trans", "clsi_bld_db_trans", "clsi_inab_db_trans", "clsi_pneu_db_trans", "clsi_skin_db_trans")
clsi_set_names <- c("clsi_all_set", "clsi_bld_set", "clsi_inab_set", "clsi_pneu_set", "clsi_skin_set")
for (i in seq_along(clsi_trans_names)) {
  data <- get(clsi_trans_names[i])
  # minimum support
  min_sup <- 1 / length(data)
  sets <- apriori(data,
                  parameter = list(support = min_sup,
                                   maxlen = length(am_col),
                                   minlen = 2,
                                   target = "frequent itemsets"))
  # Quality Measures
  CrossSupRatio <- interestMeasure(sets,
                                   "crossSupportRatio",
                                   data,
                                   reuse = TRUE)
  
  sets@quality$CrossSupRatio <- CrossSupRatio
  
  lift <- interestMeasure(sets,
                          "lift",
                          data,
                          reuse = TRUE)
  
  sets@quality$lift <- lift
  
  # save
  assign(clsi_set_names[i], sets)
  rm(data, min_sup, sets, CrossSupRatio, lift)
}

# Casey's code to convert to DF
all_sets <- function (best_setNames, labelNames){
  #required packages
  require(tidyr)
  require(stringr)
  
  #first, sets must be transformed from class itemsets to data frame. uses setsAsDataFrame function
  setsAsDataFrame <- function(sets, cat) {
    itemsets <- labels(items(sets)) #itemset names
    
    #create dataframe with category and relevant quality measures
    data.frame(
      Category <- rep(cat, length(sets)),
      items <- itemsets,
      support <- quality(sets)$support,
      count <- quality(sets)$count,
      csr <- quality(sets)$CrossSupRatio,
      lift<-quality(sets)$lift
    )
  }
  
  
  #for storing itemsets
  all.sets <- data.frame()
  
  #turn each itemset database into dataframe and label with dataset name, then combine with other itemset dataframes
  for (i in seq_along(best_setNames)){ 
    dat <- get(best_setNames[i])
    cat <- labelNames[i] #labels for each dataset
    all.sets <- rbind(all.sets, setsAsDataFrame(dat, cat))
  }
  #appropriate column names
  colnames(all.sets) <- c("Category", "items", "support", "count", "csr", "lift")
  
  #calculate the itemset order (number of AM in the set)
  all.sets$order <- str_count(all.sets$items, ",")+1 #in the itemset string, AM are divided by ",": count the commas and add 1 for the order
  
  #tabulate quality measures for each itemset in each database; display with each unique itemset as a row and the QM value in each category as the columns
  all.sets.sup <- select(all.sets, "Category", "items", "support", "order") #support in each itemset. drop other QM columns
  all.sets.sup <- spread(all.sets.sup, Category, support, drop=TRUE) #create one column for each category, place support of itemset, within category, in those columns. reduces to one row per itemset
  
  #repeat for other QM
  
  all.sets.csr <- select(all.sets, "Category", "items", "csr", "order") #csr in each itemset
  all.sets.csr <- spread(all.sets.csr, Category, csr, drop=TRUE)
  
  all.sets.lift <- select(all.sets, "Category", "items", "lift", "order") #lift in each itemset
  all.sets.lift <- spread(all.sets.lift, Category, lift, drop=TRUE)
  
  #save all dataframes in a list
  all.sets.out <- list(all.sets, all.sets.sup, all.sets.csr, all.sets.lift) #put all results into a list
  
  names(all.sets.out) <- c("all.sets", "all.sets.sup", "all.sets.csr", "all.sets.lift")
  
  all.sets.out #return the list of dataframes
}


# Split list into individual DFs
all_mine_df <- c("eu_all_set", "eu_bld_set", "eu_inab_set", "eu_pneu_set", "eu_skin_set",
                 "clsi_all_set", "clsi_bld_set", "clsi_inab_set", "clsi_pneu_set", "clsi_skin_set")

df_labels <- c("eu_all_inf", "eu_bld_inf", "eu_inab_inf", "eu_pneumo_inf", "eu_skin_tissue_inf",
               "clsi_all_set", "clsi_bld_inf", "clsi_inab_inf", "clsi_pneumo_inf", "clsi_skin_tissue_inf")

all_df <- all_sets(all_mine_df, df_labels)

split_labels <- c("df_category", "df_support", "df_csr", "df_lift")
for (i in seq_along(all_df)) {
  x <- as.data.frame(all_df[i])
  label <- paste0(as.character(split_labels[i]))
  assign(label, x)
  rm(label, x)
}

save(df_category, file = "Data/df_category.RData")
save(df_support, file = "Data/df_support.RData")
save(df_csr, file = "Data/df_csr.RData")
save(df_lift, file = "Data/df_lift.RData")
