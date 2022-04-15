library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)

bams = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/noncod/07_fwd.bam")
path = "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/noncod"
j=2
for (i in 8:30){
  bams[j] = paste(path,paste(i,"_fwd.bam",sep=""),sep="/")
  j = j+1
}
bams[2] = "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/noncod/08_fwd.bam"
bams[3] = "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/noncod/09_fwd.bam"

# bams = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/fwd.bam")
# path = "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/"
# j=2
# for (i in 8:30){
#   bams[j] = paste(path,paste(i,"/fwd.bam",sep=""),sep="")
#   j = j+1
# }
# bams[2] = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/fwd.bam")
# bams[3] = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/09/fwd.bam")


samples = 7:30
Desc = c("ctl","ctl_AUX","siCPSF3","siCPSF3_AUX","siCPSF3L","siCPSF3L_AUX","siCPSF3/3L","siCPSF3/3L_AUX")
metadata = data.frame(samples, Desc, bams)
desmat = model.matrix(~factor(metadata$Desc))

desmat = matrix(rep(c(c(1,0,0,0),c(1,0,0,1),c(1,1,0,0),c(1,1,0,1),c(1,0,1,0),c(1,0,1,1),c(1,1,1,0),c(1,1,1,1)),3),
                nrow=4)
desmat = t(desmat)
colnames(desmat) = c("Intercept", "siCPSF3", "siCPSF3L", "AUX")

frag_len = 250
win_width = 100

#do for plus and minus spearately
param = readParam(minq=10, pe="both", max.frag = 500) # important - reads extracted from bams with these params
demo <- windowCounts(bams, ext=250, param=param, BPPARAM=MulticoreParam(), width = win_width, filter=10)

binned <- windowCounts(bams, bin=TRUE, width=10000, param=param, BPPARAM=MulticoreParam())
demo = normFactors(binned, se.out=demo, assay.id=1)

##average counts and stuff
abundances = aveLogCPM(asDGEList(demo,assay.id=1), lib.size=mean(demo$totals))

summary(abundances)

y = asDGEList(demo, assay.id=1)
keep = rowSums(cpm(y) > cpm(5, mean(y$samples$lib.size))[[1]]) >= 8 #cpm of 5 in at least 8 samples
y = y[keep,]
y = estimateDisp(y,desmat)
fit = glmQLFit(y, desmat, robust=T)
summary(fit$var.post)

pdf("dispersions_plus.pdf")
par(mfrow=c(1,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)
dev.off()

results = glmQLFTest(fit, coef="AUX")

#head(results$table)
demo = demo[keep,]

rowData(demo) <- cbind(rowData(demo), results$table)
saveRDS(demo, file = "allwins_plus.rds")

#check
#rowRanges(demo)[order(mcols(rowRanges(demo))$PValue),]

##### MDS
#par(mfrow=c(2,2), mar=c(5,4,2,2))
#adj.counts <- cpm(y, log=TRUE)
#for (top in c(100, 500, 1000, 5000)) {
#  plotMDS(adj.counts, main=top, col=rep(c("blue", "green", "red", "black"),3),
#          labels=c("SiCTL.1","siCPSF3.1","siCPSF3L.1","siCPSF3/3L.1",
#                   "SiCTL.2","siCPSF3.2","siCPSF3L.2","siCPSF3/3L.2",
#                   "SiCTL.3","siCPSF3.3","siCPSF3L.3","siCPSF3/3L.3"), top=top)
#}
###

#so now you need to cluster/join windows

merged = mergeWindows(demo, 2000, ignore.strand = T) #2nd is the max dist between adjacent windows
merged$regions

summary(width(merged$regions))

tabcom <- combineTests(merged$ids, rowData(demo))
is.sig.region <- tabcom$FDR <= 0.001
summary(is.sig.region)

mcols(merged$regions) = data.frame(tabcom$PValue, tabcom$FDR, tabcom$direction)

#check
fdrsig = merged$regions[mcols(merged$regions)$tabcom.FDR <= 0.001 & mcols(merged$regions)$tabcom.direction=="up"]
#merged$regions[mcols(merged$regions)$tabcom.PValue <= 0.01 & mcols(merged$regions)$tabcom.direction=="up"]

downsig = merged$regions[mcols(merged$regions)$tabcom.FDR <= 0.001 & mcols(merged$regions)$tabcom.direction=="down"]
##
fdrsig

#export.bed(downsig ,'regions_filt10_plus_down.bed')

df <- data.frame(seqnames=seqnames(fdrsig),
                 starts=start(fdrsig)-1,
                 ends=end(fdrsig),
                 scores=fdrsig$tabcom.FDR)

write.table(df, file="regions_filt10_fdr_plus.bed", quote=F, sep="\t", row.names=F, col.names=F)
