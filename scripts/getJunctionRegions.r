args = commandArgs(trailingOnly = T)
INTRONFILE = args[1]
READLENGTH = as.numeric(args[2])
OUTDIR=args[3]

allintrons = read.table(file=INTRONFILE, header=F, sep="\t")
colnames(allintrons) <- c("chr","start","end","name","type","strand")

exonsup <- exonsdown <-  introns <- c()
for(i in 1:nrow(allintrons)){
  #print(i)
  exonsup[i] <- paste(allintrons$chr[i], allintrons$start[i]-10, allintrons$start[i], 
                      allintrons$name[i], ".", allintrons$strand[i], sep="_")
  exonsdown[i] <- paste(allintrons$chr[i], allintrons$end[i], allintrons$end[i]+10,
                        allintrons$name[i], ".", allintrons$strand[i], sep="_")
  introns[i] <- paste(allintrons$chr[i], allintrons$end[i]-(READLENGTH-10), allintrons$end[i]-10,
                      allintrons$name[i], ".", allintrons$strand[i], sep="_")
  if(allintrons$strand[i] == "-"){
    exonsup[i] <- paste(allintrons$chr[i], allintrons$end[i], allintrons$end[i]+10,
                        allintrons$name[i], ".", allintrons$strand[i], sep="_")
    exonsdown[i] <- paste(allintrons$chr[i], allintrons$start[i]-10, allintrons$start[i],
                          allintrons$name[i], ".", allintrons$strand[i], sep="_")
    introns[i] <- paste(allintrons$chr[i], allintrons$start[i]+10, allintrons$start[i]+(READLENGTH-10),
                        allintrons$name[i], ".", allintrons$strand[i], sep="_")
  }
}

write.table(as.data.frame(matrix(unlist(strsplit(introns, split="_")), byrow=T,ncol=6)), 
            file=paste0(OUTDIR,"readlength", READLENGTH, "_intron_ie.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(as.data.frame(matrix(unlist(strsplit(exonsup, split="_")), byrow=T,ncol=6)), 
            file=paste0(OUTDIR,"exonup_ee.bed"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(as.data.frame(matrix(unlist(strsplit(exonsdown, split="_")), byrow=T,ncol=6)), 
            file=paste0(OUTDIR,"exondown_ee.bed"),sep="\t",quote=F,row.names=F,col.names=F)
