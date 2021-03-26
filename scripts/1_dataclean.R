# Load Libraries
library(tidyverse)

# Load Functions
source("scripts/dataclean_fun.R", local = FALSE, echo = TRUE, spaced = TRUE)

# Load Data
## main DF
ast_data <- read.table("Data/Linelist SENTRY site_136.csv", header=TRUE, skipNul =TRUE, sep=",", na.strings = c(""))


## Breakpoint tables
bp_ec_eucast <- read.csv("Data/ec_bp_eucast.csv")
bp_ec_clsi <- read.csv("Data/ec_bp_clsi.csv")


# Clean data
## Select E.coli data
ast_ecoli <- ast_data %>% 
  filter(.[, 11] == "EC") %>% 
  select(which(colSums(!is.na(.)) > 0)) %>% 
  select(c(1, 4, 11, 17:57))


## Removes UTI data
ast_ecoli <- ast_ecoli %>% 
  filter(.[, 4]  != "urinary tract infection")


## Tabulate the percent of missing data for each anti-microbials by study year
ec_tab <- ast_ecoli %>% 
  select(c(2, 5:length(.))) %>% 
  group_by(.[, 1]) %>% 
  summarise(across(2:(length(.) - 1), pct_missing))


## Tabulate only data that was tested in every year from 2008-2018
drugs <- ec_tab %>% 
  select(-c(1)) %>% 
  summarise(across(.cols = everything(), sum)) %>% 
  select(which(. <= 0))


## Index of antimicrobials that was tested
am_included <- names(drugs)


## DF with only antimitcrobials that was tested
ast_ecoli_cleaned <- ast_ecoli %>% 
  select(c(1, 4), all_of(am_included))


## Drop am without s and r, NSbp is <= from breakpoints tables
bp_ec_eucast <- bp_ec_eucast %>% 
  filter(!is.na(.[, 7]) & !is.na(.[, 6])) %>% 
  mutate('NSbp' = as.numeric(.[, 7]))

bp_ec_clsi <- bp_ec_clsi %>% 
  filter(!is.na(.[, 4]) & !is.na(.[, 6])) %>% 
  mutate('NSbp' = as.numeric(.[, 4]))


## Numerical index for mic to interpretation function
eucast_am_bp_ind <- match(bp_ec_eucast$Antimicrobial, names(ast_ecoli_cleaned))
eucast_am_bp_ind <- eucast_am_bp_ind[!is.na(eucast_am_bp_ind)]

clsi_am_bp_ind <- match(bp_ec_clsi$Antimicrobial, names(ast_ecoli_cleaned))
clsi_am_bp_ind <- clsi_am_bp_ind[!is.na(clsi_am_bp_ind)]


## Runs Casey's function to interpret breakpoints into TRUE (resistant), and FALSE (susceptible)
names(bp_ec_eucast)[1] <- "Antimicrobial"
ec_interp_eucast <- MIC_to_interpretation(ast_ecoli_cleaned, eucast_am_bp_ind, bp_ec_eucast)
ec_interp_clsi <- MIC_to_interpretation(ast_ecoli_cleaned, clsi_am_bp_ind, bp_ec_clsi)


## Drop not logical values
ec_interp_eucast <- ec_interp_eucast %>% 
  select(1:2, names(ec_interp_eucast[sapply(ec_interp_eucast, is.logical)]))

ec_interp_clsi <- ec_interp_clsi %>% 
  select(1:2, names(ec_interp_clsi[sapply(ec_interp_clsi, is.logical)]))


## Group AM classes to determine MDR with greater than 3 resistances, main DF
eu_mdr_df <- ec_interp_eucast %>% 
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
    mdr = rowSums(across(c(4, 9, 11, 14:18)) == TRUE, na.rm = TRUE),
    mdro = if_else(mdr >= 3, TRUE, FALSE)
  )

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
    mdr = rowSums(across(c(4, 8:10, 12, 14:18)) == TRUE, na.rm = TRUE),
    mdro = if_else(mdr >= 3, TRUE, FALSE)
  )


## Vector for infection types
infection_types <- c("bloodstream infection", "intra-abdominal infection", "pneumonia in hospitalized patients", "skin/soft tissue infection")


## Generate prevalence tables
eu_prevalence_data <- eu_mdr_df %>% 
  select(2:15, 20) %>% 
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
    MDR = 1 - fun_sum_ab(mdro)
  )



clsi_prevalence_data <- clsi_mdr_df %>% 
  select(2:15, 20) %>% 
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


# Data Splitting
inf_type_split <- split(eu_mdr_df, eu_mdr_df$Infection.Type)
inf_type_split_clsi <- split(clsi_mdr_df, clsi_mdr_df$Infection.Type)


## Rename with original antimicrobial names
eu_abx_names <- eu_mdr_df %>% 
  select(3:15) %>% 
  names()

clsi_abx_names <- clsi_mdr_df %>% 
  select(3:15) %>% 
  names()


## eucast split  
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

rm(ast_data, ast_ecoli, ast_ecoli_cleaned, drugs, inf_type_split, inf_type_split_clsi)
