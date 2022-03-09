### a script used to generate Figure 2
# Nadezhda Belonogova, 2022

#### prepare summary stats

library(sumFREGAT)
library(data.table)

df <- fread('../sumstats_neuro_sum_ctg_format.txt.gz', h = T, data.table = F) # neuroticism summary statistics from https://ctg.cncr.nl/software/summary_statistics
df <- df[, c('CHR', 'POS', 'SNP', 'A1', 'P', 'BETA', 'EAF')]
colnames(df)[3:4] <- c('ID', 'EA')
#fwrite(df, '4prep.txt', row = F, qu = F, sep = ' ', na = 'NA')

load('ref.ldstore265.all.RData') # 'ref', available at https://mga.bionet.nsc.ru/sumfregat/ukbb/
#df <- fread('4prep.txt', h = T, data.table = F)
df$rs <- df$ID
df <- df[df$rs %in% ref$rs, ]
df$ID <- ref$ID[match(df$rs, ref$rs)]
df$ref <- ref$REF[match(df$rs, ref$rs)]
df$alt <- ref$ALT[match(df$rs, ref$rs)]
fwrite(df, '4prep.ukbb.txt', row = F, qu = F, sep = ' ', na = 'NA') # only SNPs with correlations available
prep.score.files('4prep.ukbb.txt', reference.file = ref, output.file.prefix = 'scores.ukbb') # generates 'scores.ukbb.vcf.gz'


####### run the benchmarks

## no approximation

invisible(sapply(list.files(pattern = "[.]R$", path = "/home/lima/nadya/FREGAT/sumFREGAT_1.2.1_no_approx_benchmark/R/", full.names = TRUE), source))
# a version with uncommented lines in genewise.R

geneFile <- "/home/lima/nadya/FREGAT/sumFREGAT_1.0.1/inst/testfiles/NCBI37.3.geneFile.txt" 
scoreFile <- 'scores.ukbb.vcf.gz'
corpath <- '265K.all.RData/' # matrices available at https://mga.bionet.nsc.ru/sumfregat/ukbb/

gf <- read.table(geneFile)
gf <- gf[gf$V3 %in% paste0('chr', 1:22), ]
genes <- gf[, 1]

reg <- SKAT(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, genes = genes, beta.par = c(1, 1), method = "kuonen", write.file = paste0('SKAT11.time.out'), quiet = T)


## with approximation, a grid of threshold values

invisible(sapply(list.files(pattern = "[.]R$", path = "/home/lima/nadya/FREGAT/sumFREGAT_1.2.2_benchmark_threshold/R/", full.names = TRUE), source))
# a version with threshold to be set in options

geneFile <- "/home/lima/nadya/FREGAT/sumFREGAT_1.2.2/inst/testfiles/NCBI37.3.geneFile.txt" 
scoreFile <- 'scores.ukbb.vcf.gz'
corpath <- '265K.all.RData/' # matrices available at https://mga.bionet.nsc.ru/sumfregat/ukbb/

ini <- read.table('SKAT11.time.out', h = T)
genes <- ini$gene[which(ini$fil <= 100)]
genes <- ini$gene[ini$fil <= 10000]

thr0 <- c(0.05, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

for (thr in thr0) {

	reg <- SKAT(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, genes = genes, beta.par = c(1, 1), method = "kuonen", threshold = thr, write.file = paste0('SKAT11.approx.', thr, '.out'), quiet = T)

}


##### data for Figure 2

ini <- read.table('SKAT11.time.out', h = T)
genes <- ini$gene[ini$fil >= 500 & ini$fil <= 10000]
ini <- ini[ini$gene %in% genes, ]

tsum <- trel <- r2 <- up <- down <- c()
df0$ini <- data.frame(gene = ini$gene, m = ini$fil, p = ini$pval, t = ini$gene.elapsed)
logp.ini <- -log10(ini$pval)

thr0 = c(0.05, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
n <- length(thr0) - 1
thr0.names <- c('t05', 't30', 't40', 't50', 't60', 't65', 't70', 't75', 't80', 't85', 't90', 't95')


for (i in 1:n) {

	ap <- read.table(paste0('SKAT11.approx.', thr0[i], '.out'), h = T)}
	print(table(ap[, 1] == ini[, 1]))
	tsum[i] <- sum(ap$gene.elapsed)
	trel[i] <- mean(1 - ap$gene.elapsed/ini$gene.elapsed)
	logp.ap <- -log10(ap$pval)
	dif <- logp.ap - logp.ini
	up[i] <- sqrt(sum(dif[dif > 0] ^2))
	down[i] <- sqrt(sum(dif[dif < 0] ^2))
	a <- summary(lm(logp.ap ~ logp.ini))
	r2[i] <- a$adj.r
	df0[[as.character(thr0.names[i])]] <- data.frame(gene = ap$gene, m = ap$fil, p = ap$pval, t = ap$gene.elapsed, trel = 1 - ap$gene.elapsed/ini$gene.elapsed)

}

i <- length(thr0) # thr = 1, no approximation
tsum[i] <- sum(ini$gene.elapsed)
trel[i] <- up[i] <- down[i] <- 0
r2[i] <- 1


#### draw Figure 2

setEPS()

postscript('fig2.eps', width = 10.5, height = 5)

par(mar = c(5.1, 4.1, 4.1, 4.1))
par(mfrow = c(1, 2))
par(lwd = 2)

plot(thr, r2, typ = 'l', xlab = 'Threshold', ylab = '', col = 'red', main = 'A')
points(thr, r2, pch = 19, col = 'red')
mtext("R2 (red)", side = 2, line = 2, col = 'red')
par(new = TRUE)
plot(thr, tsum, typ = 'l', axes = FALSE, col = 'black', xlab = '', ylab = '')
points(thr, tsum, pch = 19, col = 'black')
axis(side = 4, at = pretty(range(tsum)))
mtext("Total time, sec (black)", side = 4, line = 2)

lim <- c(floor(min(c(up, down))), ceiling(max(c(up, down))))
plot(thr, down, typ = 'l', xlab = 'Threshold', ylim = lim, ylab = '', col = 'blue', main = 'B')
points(thr, down, pch = 19, col = 'blue')
mtext("Conservative deviance (blue)", side = 2, line = 2, col = 'blue')
points(thr, up, typ = 'l', col = 'red')
points(thr, up, pch = 19, col = 'red')
mtext("Inflative deviance (red)", side = 2, line = 3, col = 'red')
par(new = TRUE)
plot(thr, tsum, typ = 'l', axes = FALSE, col = 'black', xlab = '', ylab = '')
points(thr, tsum, pch = 19, col = 'black')
axis(side = 4, at = pretty(range(tsum)))
mtext("Total time, sec (black)", side = 4, line = 2)

dev.off()
