# ======================================================================================#
# This Script contains functions designed to return a 
# numerical value indicating the strength of polarity, or the polatiry value
# of a word or phrase.

# The input taken by the functions are in the form of 2x2 Contingency matrix, eg.
#                              pos          neg
#                   c   |    f(c,pos)    f(c,neg)
#                   ~c  |   f(~c,pos)   f(~c,neg)
# ======================================================================================#


# ================================================================================
# Function: calcChiSqrPV
# Task: Calculate the Chi Square based Polarity Value from a contingency matrix
# Input: Numeric vector/Matrix representing the contingeny matrix
# Output: Polarity Value
# ================================================================================

calcChiSqrPV<-function(conMat){
  
  # Chi Square based polarity value
  mat<- matrix(as.numeric(conMat), ncol=2, byrow=TRUE)
  expmat<- matrix(data=NA, nrow=2, ncol=2)
  
  # Sum of all elements of contingency matrix
  sum<-0
  
  #Loop to get the sum
  for(i in 1:2){ # row
    for(j in 1:2){ #coloumn
      sum<-sum+mat[i,j]
    }
  }
  
  # For-loop to genetare matrix of expected values
  for(i in 1:2){ # row
    for(j in 1:2){ #coloumn
      rowsum<-0
      colsum<-0
      
      for(k in 1:2){
        rowsum<- rowsum+mat[i,k]
      }
      
      for(k in 1:2){
        colsum<- colsum+mat[k,j]
      }
      
      expmat[i,j]<- (rowsum*colsum)/sum
    }
  }
  
  chisqr<-0
  
  # Loop to calculate the chi-square value
  for(i in 1:2){ # row
    for(j in 1:2){ #coloumn
      
      chisqr<- chisqr+ ((mat[i,j]-expmat[i,j])^2)/expmat[i,j]
      
    }
  }
  
  Pc_pos<-mat[1,1]/(mat[1,1]+mat[2,1])
  Pc_neg<-mat[1,2]/(mat[1,2]+mat[2,2])
  
  if(Pc_pos<Pc_neg){
    chisqr<- chisqr*-1
  }
  
  return(chisqr)
  
}


# ================================================================================
# Function: calcPmiPV
# Task: Calculate the PMI based Polarity Value from a contingency matrix
# Input: Numeric vector/Matrix representing the contingeny matrix
# Output: Polarity Value
#=================================================================================

calcPmiPV<- function(conMat){
  
  # Pointwise Mutual Information
  # Contingency Table
  mat<- matrix(as.numeric(conMat), ncol=2, byrow=TRUE)
  
  # cpos, cneg, nocpos, nocneg
  cpos<-mat[1,1]
  cneg<-mat[1,2]
  nocpos<-mat[2,1]
  nocneg<-mat[2,2]
  
  sum<-cpos+cneg+nocpos+nocneg
  
  # 
  Pcpos<- cpos/sum
  Pcneg<- cneg/sum
  Ppos<- cpos+nocpos/sum
  Pneg<- cneg+nocneg/sum
  Pc<- cpos+cneg/sum
  
  # PMI for cpos & cneg
  PMIcpos<- log2(Pcpos/(Pc*Ppos))
  PMIcneg<- log2(Pcneg/(Pc*Pneg))
  
  # Polarity Value (PV) is PMI(c,pos) - PMI(c,neg)
  PV<- (PMIcpos-PMIcneg)
  
  return(PV)
}
