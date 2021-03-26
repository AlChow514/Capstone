# Casey's code to convert mic to interpretation


MIC_to_interpretation <- function (data, index, bp){
  #required packages
  require (stringr)
  require (dplyr)
  require (tidyr)
  require(arules)
  
  data[,index]<-sapply(data[,index], as.character) #first convert MIC to strings
  data[,index] <-sapply(data[,index], function(x) gsub(" ", "", x, fixed = TRUE)) #remove whitespace in the MIC columns
  
  for (i in index){ #for each MIC column
    bp.index<-match(names(data)[i],bp$Antimicrobial) #index of which AM bp should be applied to the i MIC column
    
    #Check that the dilutions tested cover the breakpoint appropriately (bp must be > smallest MIC value and < largest MIC value tested). Print warning if the bp does not fall within the tested dilutions
    #dilution_bp_check1=as.numeric(sapply(data[,i], function(x) str_match(x, "<=(.*)")[,2]))
    min_dilution_too_big = as.numeric(str_match(as.character(data[,i]), "<=(.*)")[,2])>bp$NSbp[bp.index]
    if (any(min_dilution_too_big==TRUE, na.rm=TRUE)){
      warning("minimum dilution greater than breakpoint for ", colnames(data)[i], ". ", sum(min_dilution_too_big==TRUE, na.rm=TRUE), " values replaced by NA \n")
    }
    
    max_dilution_too_small=as.numeric(str_match(as.character(data[,i]), ">(.*)")[,2])<bp$NSbp[bp.index]
    if (any(max_dilution_too_small==TRUE, na.rm=TRUE)){
      warning("maximum dilution less than breakpoint for ", colnames(data)[i], ". ", sum(max_dilution_too_small==TRUE, na.rm=TRUE)," values replaced by NA \n")
    }
    
    #turn NA into FALSE for indexing
    min_dilution_too_big <- replace_na(min_dilution_too_big, FALSE)
    max_dilution_too_small <- replace_na(max_dilution_too_small, FALSE)
    #replace MIC values with NA when it is a min dilution that is too big or a max dilution that is too small
    data[min_dilution_too_big,i] <- as.character(NA)
    data[max_dilution_too_small,i] <- as.character(NA)
    
    #this splits the MIC into two columns at "<=". The original string is in the first column returned. If there is a "<=", the MIC is in the second column. if there is not a "<=", it generates an NA in both columns
    
    #dilution_bp_check2=as.numeric(sapply(data[,i], function(x) str_match(x, ">(.*)")[,2]))
    #this splits the MIC into two columns at ">". The original string is in the first column. If there is a ">", the MIC is in the second column. if there is not a ">", it generates an NA in both columns
    
    #if any of the minimum MIC values are greater than the breakpoint, warn
    # if (any(dilution_bp_check1>bp$NSbp[bp.index], na.rm=TRUE)){
    #   warning("minimum dilution greater than breakpoint for ", colnames(data)[i], "\n")
    # }
    
    #if any of the maximum MIC values are greater than the breakpoint, warn
    # if (any(dilution_bp_check2<bp$NSbp[bp.index], na.rm=TRUE)){
    #   warning("maximum dilution less than breakpoint for ", colnames(data)[i], "\n")
    # }
    
    
    #next remove >, <=, =; remove anything that isn't a digit; then make MIC numeric
    data[,i]<-as.numeric(sapply(data[[i]],
                                function(x)
                                  str_replace_all(x, c(
                                    "<=" = "", 
                                    "=" = "",
                                    ">" = "1000",
                                    "[^.,0-9]" = "")))) #treat any ">" MIC values as large number to clearly distinguish from "-" MIC values in case different dilution series are tested for the same AM.
    
    #categorize each MIC value as S or NS based on the breakpoint
    if ((is.na(bp.index)+is.na(bp$NSbp[bp.index]))>0){ #if there bp.index for the i MIC column is NA, then there is no bp for that AM. OR if the bp is NA
      data[,i]<-NA #replace all MIC values with NA
      warning("missing breakpoint for ", colnames(data)[i], "\n") #warn that there is no bp for a column
    } else{ #if the bp.index is not NA, then there is a bp for that AM
      data[,i]<-discretize((data[[i]]), method="fixed",
                           breaks=c(-Inf, bp$NSbp[bp.index], Inf),
                           right=TRUE, labels=c("FALSE", "TRUE")) #discretize the MIC values at the bp. Intervals are right closed: (-Inf,bp]; (bp, Inf). Hence, MIC values >bp are given "True" (NS) and MIC values <=bp are given "False" (S)
      data[,i]=as.logical(data[[i]]) #make logical
    }
    
  }
  
  
  
  data #return the discretized data
}





# function to summarize percent of isolates with missing MIC
pct_missing <- function (x){
  round(sum(is.na(x)) / length(x), 2)
}

# Functions to construct prevalence table
fun_sum_ab <- function(x) {
  prev <- round(sum(x == FALSE, na.rm = TRUE) / (sum(x == TRUE, na.rm = TRUE) + sum(x == FALSE, na.rm = TRUE)), digits = 2)
  prev
}

#function to remove data with more than 25% missing
fun_remove <- function(df){
  rm_names <- df %>% 
    summarise(across(everything(),
                     pct_missing)) %>%
    select(which(. <= 0.25))
  
  rm_names <- names(rm_names)
  
  df <- df %>% 
    select(all_of(rm_names))
  
  return(df)
}

