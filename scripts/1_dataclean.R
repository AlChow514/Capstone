install.packages("checkpoint")
library(checkpoint)
checkpoint("2021-4-30")

# Load Libraries
library(tidyverse)

# Load Functions
source("scripts/dataclean_fun.R", local = FALSE, echo = TRUE, spaced = TRUE)

# Load Data
## main DF
ast_data <- read.table("Data/Linelist SENTRY site_136.csv", header=TRUE, skipNul =TRUE, sep=",", na.strings = c(""))


## Breakpoint tables
#eucast ECOFF as of 5-11-21 and breakpoint version 11.0 (1-1-21)
bp_ec_eucast <- read.csv("Data/ec_bp_eucast.csv") #note that AMC ECOFF changed to "ID" on 5-27-21; meropenem changed to 0.06 on 5-15-21
#if missing ECOFF, use S/R bp
bp_ec_eucast$ecoff <- ifelse(is.na(bp_ec_eucast$ecoff)==TRUE, bp_ec_eucast$S..., bp_ec_eucast$ecoff)

bp_ec_clsi <- read.csv("Data/ec_bp_clsi.csv")
#note: parenteral breakpoints used over oral breakpoints when applicable. M100 breakpoints current as of 7-12-21

# Clean data
## Select E.coli data
ast_ecoli <- ast_data %>% 
  filter(.[, 11] == "EC") %>% #select E. coli from Org Code
  select(which(colSums(!is.na(.)) > 0)) %>%  #select columns that are not all missing
  select(c(1, 4, 11, 17:57)) #select collection number, year, org code, infection type, and AMs


## Removes UTI data
ast_ecoli <- ast_ecoli %>% 
  filter(.[, 4]  != "urinary tract infection") #col 4 = infection type


## Tabulate the percent of missing data for each anti-microbials by study year
ec_tab <- ast_ecoli %>% 
  select(c(2, 5:length(.))) %>% #year and AMs
  group_by(.[, 1]) %>%  #group by year
  summarise(across(2:(length(.) - 1), pct_missing)) #summarise AMs (last column is grouping column ",[,1]")
#note: gives reminder about summarise ungrouping but output ok
#gives 0 for no missing data in the year-AM combo, 1 for all missing


## Tabulate only data that was tested in every year from 2008-2018
drugs <- ec_tab %>% 
  select(-c(1)) %>%  #drop the group column
  summarise(across(.cols = everything(), sum)) %>% #sum down columns (across all years)
  select(which(. <= 0)) #keep only drug columns that = 0 (no missing data in any year)


## Index of antimicrobials that was tested consistently
am_included <- names(drugs)


## DF with only antimitcrobials that was tested
ast_ecoli_cleaned <- ast_ecoli %>% 
  select(c(1, 4), all_of(am_included)) #include collection number, infection type, and drugs


## Drop AM without s and r, NSbp is <= from breakpoints tables
bp_ec_eucast <- bp_ec_eucast %>% 
  filter(!is.na(.[, 7]) & !is.na(.[, 6])) %>% #keep only drugs with a value in R column and ECOFF column
  mutate('NSbp' = as.numeric(.[, 7])) #rename ECOFF column as the non-susceptible breakpoint (NSbp). MICs <= the NSbp will be susceptible
##CC: check that the ECOFF is wt <= ECOFF and non-wt > ECOFF

bp_ec_clsi <- bp_ec_clsi %>% 
  filter(!is.na(.[, 4]) & !is.na(.[, 6])) %>% #keep only drugs with a value in S and R columns
  mutate('NSbp' = as.numeric(.[, 4])) #S column becomes NSbp, such that MICs <= NSbp are S and MICs > NSbp are R/I


## Numerical index for mic to interpretation function. This will be the columns interpreted as S/R in the MIC_to_interpretation function
eucast_am_bp_ind <- match(bp_ec_eucast$Antimicrobial, names(ast_ecoli_cleaned)) #column index of each AM with breakpoint in the ast_ecoli_cleaned
eucast_am_bp_ind <- eucast_am_bp_ind[!is.na(eucast_am_bp_ind)] #remove NAs from AMs with breakpoints not in ast_ecoli_cleaned

clsi_am_bp_ind <- match(bp_ec_clsi$Antimicrobial, names(ast_ecoli_cleaned)) #repeat for CLSI
clsi_am_bp_ind <- clsi_am_bp_ind[!is.na(clsi_am_bp_ind)]


## Run MIC_to_interpretation to interpret breakpoints into TRUE (resistant), and FALSE (susceptible)
names(bp_ec_eucast)[1] <- "Antimicrobial" #bp dataframe needs col named "Antimicrobial" for the function
ec_interp_eucast <- MIC_to_interpretation(ast_ecoli_cleaned, eucast_am_bp_ind, bp_ec_eucast) ##many dilutions not interpretable with ECOFF; values replaced by NA
ec_interp_clsi <- MIC_to_interpretation(ast_ecoli_cleaned, clsi_am_bp_ind, bp_ec_clsi)


## Drop not logical values from drugs (i.e., those not interpreted because there isn't a ECOFF or CLSI bp)
ec_interp_eucast <- ec_interp_eucast %>% 
  select(1:2, names(ec_interp_eucast[sapply(ec_interp_eucast, is.logical)])) #keep collection number (1) and infection type (2)

ec_interp_clsi <- ec_interp_clsi %>% 
  select(1:2, names(ec_interp_clsi[sapply(ec_interp_clsi, is.logical)]))


## Group AM classes to determine MDR with greater than 3 resistances, main DF
  #if only one drug per class, don't group (Aztreonam, Gentamicin, Levofloxacin, Tigecycline, Trimeth-sulfa)
eu_mdr_df <- ec_interp_eucast %>% 
  mutate(
    BL = if_else(
      (Ampicillin.sulbactam | Piperacillin.tazobactam) == TRUE, TRUE, FALSE #beta lactam
    ),
    CS = if_else(
      (Cefepime | Ceftazidime | Ceftriaxone) == TRUE, TRUE, FALSE #cephalosporin
    ),
    CP = if_else(
      (Doripenem | Imipenem | Meropenem) == TRUE, TRUE, FALSE #carbapenem
    )
  ) %>%
  mutate(
    mdr = rowSums(across(c(4, 9, 11, 14:18)) == TRUE, na.rm = TRUE), #sum TRUEs across aztreonam, gentamicin, levo, tigecyclein, trimeth-sulfa, BL, CS, CP
    mdro = if_else(mdr >= 3, TRUE, FALSE) #MDR is >=3 resistances
  )

#repeat for CLSI
clsi_mdr_df <- ec_interp_clsi %>% 
  mutate(
    BL = if_else(
      (Ampicillin.sulbactam | Piperacillin.tazobactam) == TRUE, TRUE, FALSE
    ),
    CS = if_else(
      (Cefepime | Ceftazidime | Ceftriaxone) == TRUE, TRUE, FALSE
    ),
    CP = if_else(
      (Doripenem | Imipenem | Meropenem) == TRUE, TRUE, FALSE
    )
  ) %>%
  mutate(
    mdr = rowSums(across(c(4, 8:10, 12, 14:18)) == TRUE, na.rm = TRUE), #sum TRUEs across aztreonam, 
    ##CC: error? 8 is doripenem, which is also counted in CP (18). shouldn't affect results b/c only 3 isolates have Doripenem R and they are MDRO without considering the carbapenems
    mdro = if_else(mdr >= 3, TRUE, FALSE)
  )


## Vector for infection types
infection_types <- c("bloodstream infection", "intra-abdominal infection", "pneumonia in hospitalized patients", "skin/soft tissue infection")
#drop other infection type levels
eu_mdr_df$Infection.Type <- droplevels(eu_mdr_df$Infection.Type)
clsi_mdr_df$Infection.Type <- droplevels(clsi_mdr_df$Infection.Type)

## Generate prevalence tables
#fun_sum_ab calculates the percent susceptible (FALSE)
eu_prevalence_data <- eu_mdr_df %>% 
  select(2:15, 20) %>%  #columns: infection type, all drugs, mdro
  group_by(Infection.Type) %>% 
  summarize(
    N = n(),
    SAM = fun_sum_ab(Ampicillin.sulbactam),
    AZT = fun_sum_ab(Aztreonam),
    CEP = fun_sum_ab(Cefepime),
    CAZ = fun_sum_ab(Ceftazidime),
    CRO = fun_sum_ab(Ceftriaxone),
    DOR = fun_sum_ab(Doripenem),
    GM = fun_sum_ab(Gentamicin),
    IMI = fun_sum_ab(Imipenem),
    LVX = fun_sum_ab(Levofloxacin),
    MEM = fun_sum_ab(Meropenem),
    TZP = fun_sum_ab(Piperacillin.tazobactam),
    TGC = fun_sum_ab(Tigecycline),
    SXT = fun_sum_ab(Trimethoprim.sulfamethoxazole),
    MDR = 1 - fun_sum_ab(mdro) #MDR is TRUE, so need to do 1- fun_sum_ab
  )



clsi_prevalence_data <- clsi_mdr_df %>% 
  select(2:15, 20) %>% #columns: infection type, drugs, mdro
  group_by(Infection.Type) %>% 
  summarize(
    N = n(),
    SAM = fun_sum_ab(Ampicillin.sulbactam),
    AZT = fun_sum_ab(Aztreonam),
    CEP = fun_sum_ab(Cefepime),
    CAZ = fun_sum_ab(Ceftazidime),
    CRO = fun_sum_ab(Ceftriaxone),
    DOR = fun_sum_ab(Doripenem),
    GM = fun_sum_ab(Gentamicin),
    IMI = fun_sum_ab(Imipenem),
    LVX = fun_sum_ab(Levofloxacin),
    MEM = fun_sum_ab(Meropenem),
    TZP = fun_sum_ab(Piperacillin.tazobactam),
    DOX = fun_sum_ab(Doxycycline),
    SXT = fun_sum_ab(Trimethoprim.sulfamethoxazole),
    MDR = 1 - fun_sum_ab(mdro)
  )


# Data Splitting, divide into dataframes by infection type
inf_type_split <- split(eu_mdr_df, eu_mdr_df$Infection.Type)
inf_type_split_clsi <- split(clsi_mdr_df, clsi_mdr_df$Infection.Type)


## Rename with original antimicrobial names (replace abbreviations)
eu_abx_names <- eu_mdr_df %>% 
  select(3:15) %>% #drug columns
  names()

clsi_abx_names <- clsi_mdr_df %>% 
  select(3:15) %>% #drug columns
  names()


## eucast split, rename drugs in each dataframe and save separately. cols 3:15 are drugs
eu_bld_db <- as.data.frame(inf_type_split[1]) %>% 
  select(3:15) %>% 
  setNames(eu_abx_names) 

eu_inab_db <- as.data.frame(inf_type_split[2]) %>% 
  select(3:15) %>% 
  setNames(eu_abx_names) 

eu_pneu_db <- as.data.frame(inf_type_split[3]) %>% 
  select(3:15) %>% 
  setNames(eu_abx_names) 

eu_skin_db <- as.data.frame(inf_type_split[4]) %>% 
  select(3:15) %>% 
  setNames(eu_abx_names)

eu_all_db <- eu_mdr_df %>% 
  select(3:15) %>% 
  setNames(eu_abx_names) 

# removes drugs with more than 25% missing
eu_bld_db <- fun_remove(eu_bld_db)
eu_inab_db <- fun_remove(eu_inab_db)
eu_pneu_db <- fun_remove(eu_pneu_db)
eu_skin_db <- fun_remove(eu_skin_db)
eu_all_db <- fun_remove(eu_all_db)

## clsi split
clsi_bld_db <- as.data.frame(inf_type_split_clsi[1]) %>% 
  select(3:15) %>% 
  setNames(clsi_abx_names)

clsi_inab_db <- as.data.frame(inf_type_split_clsi[2]) %>% 
  select(3:15) %>% 
  setNames(clsi_abx_names)

clsi_pneu_db <- as.data.frame(inf_type_split_clsi[3]) %>% 
  select(3:15) %>% 
  setNames(clsi_abx_names)

clsi_skin_db <- as.data.frame(inf_type_split_clsi[4]) %>% 
  select(3:15) %>% 
  setNames(clsi_abx_names)

clsi_all_db <- clsi_mdr_df %>% 
  select(3:15) %>% 
  setNames(clsi_abx_names)