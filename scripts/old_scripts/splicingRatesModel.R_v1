#!/usr/bin/env Rscript

# inputs:

args = commandArgs(trailingOnly = T)
JUNCTIONFILE = args[1] # file with junction counts, Dprime, and Rprime
TXNRATE = as.numeric(args[2]) #txn rate in nt/min
#TXNRATE = 1500
#OUTDIR=args[3]


##### function ######

# jointly model half-lives from all samples
sumsqequationsolve <- function(x = "file with junction counts, Dprimes, and Rprimes", txnrate = "transcription rate"){
  
  # separate out D' and R' values
  D_prime = x[c(grep("^Dprime_[0-9]+$", names(x)))]
  R_prime = x[c(grep("^Rprime_[0-9]+$", names(x)))]

  i = 1:length(R_prime)

  # build function as sum of functions for each sample
  hl <- function(h){sum(((h*(1 - 2^(-D_prime[i]/(h*txnrate)))) - R_prime[i])^2)}
  
  # initialize values for running optimization
  hold.row <- c(NA, NA)
  #fit.hold <- c(NA,NA)
  starth = 0
  if(sum(is.na(R_prime)) == length(R_prime)){ return(hold.row) } ### left it here but this conditional doesn't seem to make a difference 
  # try running function
  try(fit.hold <- optim(starth, hl), silent=T)
  try(hold.row <- c(fit.hold$par, fit.hold$value), silent=T)
  return(hold.row)
}


###### running script #####

# get file with junction counts, Dprime, and Rprime
juncratio.data <- read.table(JUNCTIONFILE, sep="\t", header=T)

# run model
sumsqfit.data <- t(apply(juncratio.data, 1, sumsqequationsolve, TXNRATE))

# fitvalue = half live, yvalue = error around fit from optimization function
juncratio.data$fitvalue <- sumsqfit.data[,1]
names(juncratio.data)[names(juncratio.data) == 'fitvalue'] <- 'half.life'
juncratio.data$yvalue <- sumsqfit.data[,2]
names(juncratio.data)[names(juncratio.data) == 'yvalue'] <- 'fit.error'

##### Parse modeled introns to a set with optimal power to confidently detect half-lives
juncratio.data$result <- "success"
ratioinds <- grep('_ratio', colnames(juncratio.data))

# Remove those without intron-exon junction coverage in the earliest timepoint
## determined by ratio_T1 > 0 (if ie = 0, then ratio = 0)
#juncratio.data.parsed <- juncratio.data[which(juncratio.data[,ratioinds[1]] > 0),]
juncratio.data$result[which(juncratio.data[,ratioinds[1]] > 0)] <- "noIE_veryFast"

# Remove those without exon-exon junction reads in the last timepoint
## determined by ratio_Tn = inf
#juncratio.data.parsed <- juncratio.data[which(!is.infinite(juncratio.data[,ratioinds[length(ratioinds)]])),]
juncratio.data$result[which(is.infinite(juncratio.data[,ratioinds[length(ratioinds)]]))] <- "noEE_verySlow"

# Tag those where fit failed
#juncratio.data.parsed <- subset(juncratio.data, !is.na(fitvalue))
juncratio.data$result[which(is.na(juncratio.data$half.life))] <- "no_fit"

# Keep only those with coverage in at least n timepoints (where we used n=2)
## to be added?

##### write out final file with splicing half lives #####
#write.table(juncratio.data.parsed, file=JUNCTIONFILE, quote=F, sep="\t", row.names=F, col.names=T)
write.table(juncratio.data, file=JUNCTIONFILE, quote=F, sep="\t", row.names=F, col.names=T)
