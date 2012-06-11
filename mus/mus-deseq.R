library(DESeq)

c1 <- read.table(file="Control-3T3.fq.filtered.fq.trimmed.single.counts.txt", row.names=1)
c2 <- read.table(file="Stimulated_3T3.fq.filtered.fq.trimmed.single.counts.txt", row.names=1)

counts <- cbind(c1, c2)

colnames(counts) <- c("Ctrl","Stim")
conds <- factor( c("control", "stimulated") )

cds <- newCountDataSet( counts, conds )

cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind", sharingMode='fit-only', fitType="local" )

res <- nbinomTest(cds, "control", "stimulated")

write.table( res, file="mus-de-binom-table.txt", sep="\t" )

