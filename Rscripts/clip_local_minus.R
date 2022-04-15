#adjust params to find peaks in clip seq data

library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)

bams = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/70_read2.sorted.bam", 
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/71_read2.sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/72_read2.sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/73_read2.sorted.bam")

bamsplus = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/70_read2_fwd.sorted.bam", 
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/71_read2_fwd.sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/72_read2_fwd.sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/73_read2_fwd.sorted.bam")

bamsminus = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/70_read2_rev.sorted.bam", 
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/71_read2_rev.sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/72_read2_rev.sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/73_read2_rev.sorted.bam")

samples = 1:4
Desc = c("INPUT","ZC3H4","ZC3H4","ZC3H4")
metadata = data.frame(samples, Desc, bams)
desmat = model.matrix(~factor(metadata$Desc))

colnames(desmat) = c("Intercept", "ZC3H4")
win_width = 50

#so, the data are kind of single-end now, but still stranded - so use
#strandedcounts but change the pe parameter
param4neighbour = readParam(minq=10, pe="none")
param = readParam(minq=10, pe="none", forward=logical(0)) # important - reads extracted from bams with these params
demo <- strandedCounts(bamsminus, ext=50, spacing=20, filter=5, param=param, BPPARAM=MulticoreParam(), width = win_width, bin=F)

#local background
surrounds <- 750
neighbor <- suppressWarnings(resize(rowRanges(demo), surrounds, fix="center"))
wider <- regionCounts(bamsminus, regions=neighbor, ext=50, param=param4neighbour)

#what threshold?
filter.stat <- filterWindowsLocal(demo, wider, assay.data = 1, assay.back = 1)
summary(filter.stat$filter)

keep <- filter.stat$filter > log2(2)

pdf("enrichment_minus.pdf")
hist(filter.stat$filter, xlab="Log-fold change from local background", 
     breaks=100, main="", col="grey80", xlim=c(0, 5))
abline(v=log2(2), col="red", lwd=2)
dev.off()


filt_demo = demo[keep,]

#demo <- regionCounts(bams, gr_obj, ext=frag_len, param=param)
binned <- strandedCounts(bamsminus, bin=TRUE, width=10000, param=param, BPPARAM=MulticoreParam())
filt_demo = normFactors(binned, se.out=filt_demo, assay.id=1)

##average counts and stuff
abundances = aveLogCPM(asDGEList(filt_demo,assay.id=1), lib.size=mean(demo$totals))

summary(abundances)

y = asDGEList(filt_demo, assay.id=1)
keep1 = rowSums(cpm(y) > cpm(10, mean(y$samples$lib.size))[[1]]) >= 2
y = y[keep1,]
y = estimateDisp(y,desmat)
fit = glmQLFit(y, desmat, robust=T)
summary(fit$var.post)

par(mfrow=c(1,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)

results = glmQLFTest(fit, coef="ZC3H4")

head(results$table)
demo2 = filt_demo[keep1,]

rowData(demo2) <- cbind(rowData(demo2), results$table)

#check
rowRanges(demo2)[order(mcols(rowRanges(demo2))$PValue),]
sigwin = rowRanges(demo2)[mcols(rowRanges(demo2))$logFC > 2 & mcols(rowRanges(demo2))$PValue < 0.001]
sigwin

##### MDS
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts <- cpm(y, log=TRUE)
for (top in c(100, 500, 1000, 5000)) {
  plotMDS(adj.counts, main=top, col=c("blue", "green", "green", "green"),
          labels=c("SiCTL.1","ZC3H4.1","ZC3H4.2","ZC3H4.3"), top=top)
}
###

#so now you need to cluster/join windows

merged = mergeWindows(demo2, 500,max.width = 1000, ignore.strand = F) #2nd is the max dist between adjacent windows
merged$regions

summary(width(merged$regions))

tabcom <- combineTests(merged$ids, rowData(demo2))
is.sig.region <- tabcom$FDR <= 0.01
summary(is.sig.region)

mcols(merged$regions) = data.frame(tabcom$PValue, tabcom$FDR, tabcom$direction)

#check
fdrsig = merged$regions[mcols(merged$regions)$tabcom.FDR <= 0.01 & mcols(merged$regions)$tabcom.direction=="up"]

fdrsig

#export.bed(sigwin ,'chr2_CLIP_local_wins_nomerge.bed')

df_sigwin <- data.frame(seqnames=seqnames(sigwin),
                 starts=start(sigwin)-1,
                 ends=end(sigwin),
                 name=sigwin$logFC,
                 scores=sigwin$PValue,
                 strand=strand(sigwin))


df_fdrsig <- data.frame(seqnames=seqnames(fdrsig),
                 starts=start(fdrsig)-1,
                 ends=end(fdrsig),
                 name=fdrsig$tabcom.PValue,
                 scores=fdrsig$tabcom.FDR,
                 strand=strand(fdrsig))

write.table(df_sigwin, file="sig_win_minus.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(df_fdrsig, file="merged_win_minus.bed", quote=F, sep="\t", row.names=F, col.names=F)

