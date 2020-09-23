#---
#  title: "AMR Data Cleaning"
#author: "Casey Cazer"
#date: Copied 11/18/2019 from Association Mining Functions Git. Last Updated 2/3/2020
#---
#   
# this script contains functions to clean MIC data and prepare it for rule or set mining
# 
# requires: .csv of breakpoints (interpretive criteria). AM names in csv must be same as AM names in MIC data. The breakpoint is the only other column required. It should be the non-susceptible breakpoint such that bacteria are NS if they are > breakpoint. Therefore, the breakpoint value is the highest MIC for susceptible bacteria. 
# 
# Abbreviations: 
#   MIC: Minimum inhibitory concentration
# bp: breakpoint, i.e. interpretive criteria for MIC value
# S: susceptible
# NS: non-susceptible (may be used interchangeably with resistant/R)
# AM: antimicrobial
  

################################
#convert MIC values to interpretation criteria (S vs NS). this function removes "<=", "=", ">", and whitespace. modifications required if there are other symbols in MIC values

#data: dataset to analyze. suggested to exclude AM that won't be analyzed
#index: indices of MIC columns to convert to numeric and then interpret as S or NS.
#bp: csv of breakpoints with column called Antimicrobial and column called NSbp. The Antimicrobial column needs to include The NSbp should be the highest susceptible MIC value. Then any MIC > NSbp will be NS.

MIC_to_interpretation <- function (data, index, bp){
  #required packages
  require (stringr)
  require (dplyr)
  require (tidyr)
  
  data[,index]<-sapply(data[,index], as.character) #first convert MIC to strings
  data[,index] <-sapply(data[,index], function(x) gsub(" ", "", x, fixed = TRUE)) #remove whitespace in the MIC columns
  
  for (i in index){ #for each MIC column
    bp.index<-match(names(data)[i],bp$Antimicrobial) #index of which AM bp should be applied to the i MIC column
    
    #Check that the dilutions tested cover the breakpoint appropriately (bp must be > smallest MIC value and < largest MIC value tested). Print warning if the bp does not fall within the tested dilutions
    dilution_bp_check1=as.numeric(str_match(data[,i], "<=(.*)")[,2])
    #this splits the MIC into two columns at "<=". The original string is in the first column returned. If there is a "<=", the MIC is in the second column. if there is not a "<=", it generates an NA in both columns
    
    dilution_bp_check2=as.numeric(str_match(data[,i], ">(.*)")[,2])
    #this splits the MIC into two columns at ">". The original string is in the first column. If there is a ">", the MIC is in the second column. if there is not a ">", it generates an NA in both columns
    
    #if any of the minimum MIC values are greater than the breakpoint, warn
    if (any(dilution_bp_check1>bp$NSbp[bp.index], na.rm=TRUE)){
      warning("minimum dilution greater than breakpoint for ", colnames(data)[i], "\n")
    }
    
    #if any of the maximum MIC values are greater than the breakpoint, warn
    if (any(dilution_bp_check2<bp$NSbp[bp.index], na.rm=TRUE)){
      warning("maximum dilution less than breakpoint for ", colnames(data)[i], "\n")
    }
    
    
    #next remove >, <=, =; then make MIC numerica
    data[,i]<-as.numeric(
      str_replace_all(data[,i], c(
        "<=" = "", 
        "=" = "",
        ">" = "1000"))) #treat any ">" MIC values as large number to clearly distinguish from "-" MIC values in case different dilution series are tested for the same AM.
    
    #categorize each MIC value as S or NS based on the breakpoint
    if (is.na(bp.index)==TRUE){ #if there bp.index for the i MIC column is NA, then there is no bp for that AM
      data[,i]<-NA #replace all MIC values with NA
      warning("missing breakpoint for ", colnames(data)[i], "\n") #warn that there is no bp for a column
    } else{ #if the bp.index is not NA, then there is a bp for that AM
      data[,i]<-discretize(data[,i], method="fixed", 
                           breaks=c(-Inf, bp$NSbp[bp.index], Inf), 
                           right=TRUE, labels=c("FALSE", "TRUE")) #discretize the MIC values at the bp. Intervals are right closed: (-Inf,bp]; (bp, Inf). Hence, MIC values >bp are given "True" (NS) and MIC values <=bp are given "False" (S)
      data[,i]=as.logical(data[,i]) #make logical
    } 
    
  }
  
  
  
  data #return the discretized data
}

#######################################
#generate antibiogram: prevalence of NS for each AM; data grouped by user-defined variable

#data: dataset to analyze. must be in logical form (F=S, T=NS), i.e. from MIC_to_interpretation
#index: indices of AM columns to analyze.
#group: index of the grouping variable column used to split the database (e.g. year)

antibiogram <- function (data, index, group){
  #required packages
  require(dplyr)
  require(tidyr)
  
  abgm <- group_by(data[,c(group, index)] , data[,group]) %>%
    summarise_all(list(Prevalence=~mean(.,na.rm=TRUE),Number_Isolates_Tested=~sum(!is.na(.)))) #drop all columns except grouping column and AM columns, group data by the grouping column. Summarise each column with prevalence (# NS / # tested = mean of the logical variable with na.rm=TRUE) and number of isolates tested. This returns the prevalence and number tested within each group value. Columns include: AM_Prevalence and AM_Number_Isolates_Tested, with "AM" being each AM name
  
  #overall prevalences (across all group values)
  overall_prev <- summarise_all(data[,c(group,index)], list(Prevalence=~mean(.,na.rm=TRUE),Number_Isolates_Tested=~sum(!is.na(.)))) #same as above but without grouping
  overall_prev <- cbind("data[, group]"="Overall", overall_prev) #add in a 'group' value of "Overall"
  abgm <- rbind(abgm, overall_prev) #add overall values to the grouped antibiogram
  
  #if an AM is tested but there is no breakpoint, num of isolates tested will show 0 since the data (logical form) will have all NA. doesn't differentiate from AM that were not tested and were in dataset but with no MIC values (logical data also NA and num of isolates tested will be 0).
  #warn if overall num of isolates tested columns (second half of antibiogram) are equal to 0. These AM could be excluded from 'index'
  num_tested_col=c(((ncol(overall_prev)-1)/2+2):ncol(overall_prev)) #column indices for the "Number_Isolates_Tested" columns in the overall antibiogram
  if (any(overall_prev[num_tested_col]==0)){ #if any of these columns equal 0, warn
    warning("All interpretations NA for:", "\n", paste(str_split(names(overall_prev[num_tested_col])[overall_prev[num_tested_col]==0], "_", simplify = TRUE)[,1], collapse=", "), "\n", "May be missing breakpoint or not tested.") #list out AM names of columns with missing bp or no isolates tested
  }
  
  #sort alphabetically to put prevalence next to number tested for each drug.
  abgm <- abgm[,order(names(abgm))]
  
  abgm #return antibiogram
  
}


#######################################
#describe AM susceptibility testing panels
#code originally from Will Love (Github EpidemiologyDVM)

#data: dataset of MICs or AM interpretations. Note that if there is no breakpoint, then a dataset of AM interpretations (from MIC_to_interpretation) can exclude AM that were tested but don't have a breakpoint by excluding the NA columns.
#index: indices of MIC or AM interpretation columns

abx_profile <- function (data, index){
  ABX_NAMES <- names(data)[index] #get the AM names that may be in each isolate's susceptibility test panel
  abx_panels <- rep("", dim(data)[1])  #Initializes a blank character vector as long as the number of isolates in the data
  
  for(abx in ABX_NAMES) {	#for each AM
    x <- with(data, get(abx)) #x is column of values for the AM
    x <- ifelse(is.na(x), '', paste(abx, '')) #replace MIC values with AM name, if NA then blank
    abx_panels <- paste(abx_panels, x, sep = '') #paste the x column to the blank character vector. for each isolate, abx_panels will accumulate the names of the AM that were tested for that isolate
  }
  abx_panels #returns the AM susceptibility test profile for each isolate as a character vector of names of the AM that were tested
}


########################################
#MDR profiles
#For each isolate, identifies the resistance profile (all the drugs to which the isolate is resistant), tabulates the number of drugs to which it is resistant, identifies the class resistance profile (all the AM classes to which the isolate is resistant), the number of drug classes to which it is resistant

#data: dataset of AM interpretations. Note that if there is no breakpoint, then a dataset of AM interpretations will exclude AM that were tested but don't have a breakpoint.
#index: indices of AM columns
#AM_class is a two column dataframe. first column (name = "AM") is the name of the antimicrobial (must match AM names from data). second column (name = "Class") is the assigned class (number or character that can be sorted alphanumerically)

mdr_profile <- function (data, index, AM_class){
  #required packages
  require(stringr)
  require(mgsub)
  
  z <- data[,index] #look at only AM columns
  
  mdr_profile <- data.frame(Profile=rep("", dim(z)[1]), NumRes=numeric(dim(z)[1]), NumClass=numeric(dim(z)[1])) #create dataframe to store mdr profiles
  profile <- rep("", dim(z)[1]) #Initializes a blank character as long as the number of isolates in the data
  
  
  mdr_profile[,2] <- rowSums(z, na.rm=TRUE) #number of drug resistances for each isolate
  
  for(i in 1:ncol(z)) { #for each AM 
    z[which(z[,i]==TRUE), i] <- as.character(colnames(z)[i]) #if resistant, paste AM name into cell instead of True/False (interpreted MIC)
    z[which(z[,i]==FALSE), i] <- as.character("") #if susceptible, paste blank into cell
    z[is.na(z[,i]), i] <- as.character("") #those that weren't tested or have no bp (have NA instead of T/F), paste blank into cell
    profile <- paste(profile, z[,i], sep = " ") #paste all AM columns together to get a character string of all the AM to which each isolate (row) is resistant
  }
  profile <- trimws(profile, which="both") #trim whitespace
  mdr_profile[,1] <- profile #store the resistance profile
  
  
  #warn if the data contains an AM that is not in the AM_class--cannot proceed with class resistance tabulation
  if (any(is.na(match(colnames(z), AM_class$AM)))){ #if there are any AM names in the data that do not have a corresponding entry in AM_class, warn:
    warning("Missing AM class for ", colnames(z)[is.na(match(colnames(z), AM_class$AM))], "\n Inappropriate data from this column is pasted in results.")
  }
  
  #tabulate the class resistance profile and number of classes to which each isolate is resistant
  mdr_profile$Classes <- mgsub(mdr_profile$Profile, as.character(AM_class$AM), as.character(AM_class$Class))#in a new column, replace the name of each AM from the resistance profile with the name of the corresponding class
  
  mdr_profile$Classes.sort <- sapply(sapply(str_extract_all(mdr_profile$Classes, "[a-z0-9]+"), unique), sort) #list unique classes for each isolate then sort alphanumerically
  
  mdr_profile$NumClass <- sapply(lapply(str_extract_all(mdr_profile$Classes, "[a-z0-9]+"), unique), length) #number of unique class resistances
  
  mdr_profile #output
}