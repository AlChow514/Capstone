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
    dilution_bp_check1=as.numeric(sapply(data[,i], function(x) str_match(x, "<=(.*)")[,2]))
    #this splits the MIC into two columns at "<=". The original string is in the first column returned. If there is a "<=", the MIC is in the second column. if there is not a "<=", it generates an NA in both columns
    
    dilution_bp_check2=as.numeric(sapply(data[,i], function(x) str_match(x, ">(.*)")[,2]))
    #this splits the MIC into two columns at ">". The original string is in the first column. If there is a ">", the MIC is in the second column. if there is not a ">", it generates an NA in both columns
    
    #if any of the minimum MIC values are greater than the breakpoint, warn
    if (any(dilution_bp_check1>bp$NSbp[bp.index], na.rm=TRUE)){
      warning("minimum dilution greater than breakpoint for ", colnames(data)[i], "\n")
    }
    
    #if any of the maximum MIC values are greater than the breakpoint, warn
    if (any(dilution_bp_check2<bp$NSbp[bp.index], na.rm=TRUE)){
      warning("maximum dilution less than breakpoint for ", colnames(data)[i], "\n")
    }
    
    
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

