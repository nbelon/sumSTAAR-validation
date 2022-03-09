### a script used to generate Figure 3
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

geneFile <- "/home/lima/nadya/FREGAT/sumFREGAT_1.0.1/inst/testfiles/NCBI37.3.geneFile.txt" 
scoreFile <- 'scores.ukbb.vcf.gz'
corpath <- '265K.all.RData/' # matrices available at https://mga.bionet.nsc.ru/sumfregat/ukbb/

gf <- read.table(geneFile)
gf <- gf[gf$V3 %in% paste0('chr', 1:22), ]
genes <- gf[, 1]

reg <- SKAT(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, genes = genes, beta.par = c(1, 1), method = "kuonen", write.file = paste0('SKAT11.time.out'), quiet = T)
reg <- SKATO(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, genes = genes, beta.par = c(1, 1), method = "kuonen", write.file = paste0('SKATO11.time.out'), quiet = T)
reg <- PCA(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, n = 380506, genes = genes, write.file = paste0('PCA11.time.out'), quiet = T)
reg <- FLM(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, n = 380506, genes = genes, write.file = paste0('FLM11.time.out'), quiet = T)


## with approximation, p.threshold = 0.8

invisible(sapply(list.files(pattern = "[.]R$", path = "/home/lima/nadya/FREGAT/sumFREGAT_1.2.2_benchmark/R/", full.names = TRUE), source))

reg <- SKAT(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, genes = genes, beta.par = c(1, 1), method = "kuonen", write.file = paste0('SKAT11.approx.time.out'), quiet = T)
reg <- SKATO(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, genes = genes, beta.par = c(1, 1), method = "kuonen", write.file = paste0('SKATO11.approx.time.out'), quiet = T)
reg <- PCA(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, n = 380506, genes = genes, write.file = 'PCA11.approx.tome.out', quiet = T)
reg <- FLM(score.file = scoreFile, cor.path = corpath, gene.file = geneFile, n = 380506, genes = genes, write.file = paste0('FLM11.approx.time.out'), quiet = T)


######### output

df0 <- c()
ini <- read.table('SKAT11.time.out', h = T)
ap <- read.table('SKAT11.approx.time.out', h = T)
v <- which(ap$fil > 1); ap <- ap[v, ]; ini <- ini[v, ];table(ap$fil == ini$fil) # 17975 left
df0$skat <- data.frame(gene = ini$gene, m.ini = ini$fil, p.ini = ini$pval, t.ini = ini$gene.elapsed, m.ap = ap$fil, p.ap = ap$pval, t.ap = ap$gene.elapsed)

ini <- read.table('SKATO11.time.out', h = T)
ap <- read.table('SKATO11.approx.time.out', h = T)
v <- which(ap$fil > 1); ap <- ap[v, ]; ini <- ini[v, ];table(ap$fil == ini$fil) # 17975 left
df0$skato <- data.frame(gene = ini$gene, m.ini = ini$fil, p.ini = ini$pval, t.ini = ini$gene.elapsed, m.ap = ap$fil, p.ap = ap$pval, t.ap = ap$gene.elapsed)

ini <- read.table('PCA11.time.out', h = T)
ap <- read.table('PCA11.approx.time.out', h = T)
v <- which(ap$fil > 1); ap <- ap[v, ]; ini <- ini[v, ];table(ap$fil == ini$fil) # 17975 left
df0$pca <- data.frame(gene = ini$gene, m.ini = ini$fil, p.ini = ini$pval, t.ini = ini$gene.elapsed, m.ap = ap$fil, p.ap = ap$pval, t.ap = ap$gene.elapsed)

ini <- read.table('FLM11.time.out', h = T)
ap <- read.table('FLM11.approx.time.out', h = T)
v <- which(ap$fil >= 25); ap <- ap[v, ]; ini <- ini[v, ];table(ap$fil == ini$fil) # 7990 left
df0$flm <- data.frame(gene = ini$gene, m.ini = ini$fil, p.ini = ini$pval, t.ini = ini$gene.elapsed, m.ap = ap$fil, p.ap = ap$pval, t.ap = ap$gene.elapsed)


########### equations

### SKAT
df <- df0[['skat']]
m <- df$m.ini
m2 <- m^2
m3 <- m^3
m4 <- m^4

summary(lm(df$t.ini~m))
#Multiple R-squared:  0.4524,    Adjusted R-squared:  0.4523
#F-statistic: 1.485e+04 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m2))
#Multiple R-squared:  0.9246,    Adjusted R-squared:  0.9245
#F-statistic: 2.202e+05 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m3))
#Multiple R-squared:  0.9723,    Adjusted R-squared:  0.9723
#F-statistic: 6.304e+05 on 1 and 17973 DF,  p-value: < 2.2e-16

m <- df$m.ap
m2 <- m^2
m3 <- m^3
summary(lm(df$t.ap~m))
#Multiple R-squared:  0.554,     Adjusted R-squared:  0.554
#F-statistic: 2.232e+04 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m2))
#Multiple R-squared:  0.9369,    Adjusted R-squared:  0.9369
#F-statistic: 2.668e+05 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m3))
#Multiple R-squared:  0.8982,    Adjusted R-squared:  0.8982
#F-statistic: 1.586e+05 on 1 and 17973 DF,  p-value: < 2.2e-16

### SKATO
df <- df0[['skato']]

m <- df$m.ini
m2 <- m^2
m3 <- m^3
m4 <- m^4

summary(lm(df$t.ini~m))
#Multiple R-squared:  0.3684,    Adjusted R-squared:  0.3684
#F-statistic: 1.048e+04 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m2))
#Multiple R-squared:  0.8827,    Adjusted R-squared:  0.8827
#F-statistic: 1.353e+05 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m3))
#Multiple R-squared:  0.9999,    Adjusted R-squared:  0.9999
#F-statistic: 1.775e+08 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m4))
#Multiple R-squared:  0.9445,    Adjusted R-squared:  0.9445
#F-statistic: 3.059e+05 on 1 and 17973 DF,  p-value: < 2.2e-16

m <- df$m.ap
m2 <- m^2
m3 <- m^3
m4 <- m^4
summary(lm(df$t.ap~m))
#Multiple R-squared:  0.3881,    Adjusted R-squared:  0.3881
#F-statistic: 1.14e+04 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m2))
#Multiple R-squared:  0.8611,    Adjusted R-squared:  0.8611
#F-statistic: 1.114e+05 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m3))
#Multiple R-squared:  0.9662,    Adjusted R-squared:  0.9662
#F-statistic: 5.133e+05 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m4))
#Multiple R-squared:  0.9269,    Adjusted R-squared:  0.9269
#F-statistic: 2.28e+05 on 1 and 17973 DF,  p-value: < 2.2e-16

### PCA
df <- df0[['pca']]

m <- df$m.ini
m2 <- m^2
m3 <- m^3
m4 <- m^4

summary(lm(df$t.ini~m))
#Multiple R-squared:  0.3018,    Adjusted R-squared:  0.3018
#F-statistic:  7770 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m2))
#Multiple R-squared:  0.786,     Adjusted R-squared:  0.786
#F-statistic: 6.601e+04 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m3))
#Multiple R-squared:  0.9534,    Adjusted R-squared:  0.9534
#F-statistic: 3.676e+05 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m4))
#Multiple R-squared:  0.9518,    Adjusted R-squared:  0.9518
#F-statistic: 3.546e+05 on 1 and 17973 DF,  p-value: < 2.2e-16

m <- df$m.ap
m2 <- m^2
m3 <- m^3
m4 <- m^4
summary(lm(df$t.ap~m))
#Multiple R-squared:  0.3311,    Adjusted R-squared:  0.3311
#F-statistic:  8897 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m2))
#Multiple R-squared:  0.8011,    Adjusted R-squared:  0.8011
#F-statistic: 7.239e+04 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m3))
#Multiple R-squared:  0.9403,    Adjusted R-squared:  0.9403
#F-statistic: 2.829e+05 on 1 and 17973 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m4))
#Multiple R-squared:  0.9266,    Adjusted R-squared:  0.9266
#F-statistic: 2.269e+05 on 1 and 17973 DF,  p-value: < 2.2e-16

### FLM
df <- df0[['flm']]

m <- df$m.ini
m2 <- m^2
m3 <- m^3
m4 <- m^4

summary(lm(df$t.ini~m))
#Multiple R-squared:  0.4398,    Adjusted R-squared:  0.4398
#F-statistic:  6272 on 1 and 7988 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m2))
#Multiple R-squared:  0.9039,    Adjusted R-squared:  0.9039
#F-statistic: 7.512e+04 on 1 and 7988 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m3))
#Multiple R-squared:  0.9971,    Adjusted R-squared:  0.9971
#F-statistic: 2.729e+06 on 1 and 7988 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m4))
#Multiple R-squared:  0.9296,    Adjusted R-squared:  0.9296
#F-statistic: 1.054e+05 on 1 and 7988 DF,  p-value: < 2.2e-16
summary(lm(df$t.ini~m+m2+m3))
#Multiple R-squared:  0.9985,    Adjusted R-squared:  0.9985
#F-statistic: 1.798e+06 on 3 and 7986 DF,  p-value: < 2.2e-16

m <- df$m.ap
m2 <- m^2
m3 <- m^3
m4 <- m^4
summary(lm(df$t.ap~m))
#Multiple R-squared:  0.4868,    Adjusted R-squared:  0.4867
#F-statistic:  7577 on 1 and 7988 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m2))
#Multiple R-squared:  0.9252,    Adjusted R-squared:  0.9252
#F-statistic: 9.886e+04 on 1 and 7988 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m3))
#Multiple R-squared:  0.9813,    Adjusted R-squared:  0.9813
#F-statistic: 4.198e+05 on 1 and 7988 DF,  p-value: < 2.2e-16
summary(lm(df$t.ap~m4))
#Multiple R-squared:  0.8968,    Adjusted R-squared:  0.8967
#F-statistic: 6.938e+04 on 1 and 7988 DF,  p-value: < 2.2e-16


#### figure

meth0 <- c('skat', 'skato', 'pca', 'flm')

# pow determined by max Fstat above
pow.ini0 <- c(3, 3, 3, 3)
pow.ap0 <- c(2, 3, 3, 3)

setEPS()

postscript('fig3.eps', width = 6.3, height = 13)

par(mfrow = c(4, 2))

for (j in 1:4) {

cat(j, '')
meth <- meth0[j]
df <- df0[[as.character(meth)]]
df <- df[which(df$p.ini != 0 & df$p.ap != 0), ]

x <- -log10(df$p.ini)
y <- -log10(df$p.ap)

#lim <- ceiling(max(c(x, y)))
lim <- 31

plot(x, y, main = toupper(meth), xlab = 'Original', ylab = 'Approximation', xlim = c(0, lim), ylim = c(0, lim))
points(x[which(df$m.ini >= 500)], y[which(df$m.ini >= 500)], pch = 19)
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')

legend(x = lim / 100, y = lim - lim / 100, leg = c('< 500 variants', '>= 500 variants'), pch = c(1, 19), bty = 'n')

plot(df$m.ini, df$t.ini, main = toupper(meth), xlab = 'Number of SNPs in gene', ylab = 'Time elapsed (sec)')
points(df$m.ap, df$t.ap, xlab = 'Number of SNPs in gene', ylab = 'Time elapsed (sec)', col = 'red')

pow.ini <- pow.ini0[j]
mpow <- df$m.ini ^ pow.ini

a <- summary(lm(df$t.ini ~ mpow))
i.ini <- format(a$coef[1, 1], digits = 2)
bpow.ini <- -ceiling(abs(log(a$coef[2, 1], 10)))
b.ini <- round(a$coef[2, 1] / (10 ^ bpow.ini), 2)
r2.ini <- format(a$adj.r, digits = 3)

lim <- max(c(df$m.ini, df$m.ap))
x <- 0:(lim * 5) / 5
y.ini <- a$coef[1, 1] + a$coef[2, 1] * (x ^ pow.ini)

pow.ap <- pow.ap0[j]
mpow <- df$m.ini ^ pow.ap

a <- summary(lm(df$t.ap ~ mpow))
i.ap <- format(a$coef[1, 1], digits = 2)
bpow.ap <- -ceiling(abs(log(a$coef[2, 1], 10)))
b.ap <- round(a$coef[2, 1] / (10 ^ bpow.ap), 2)
r2.ap <- format(a$adj.r, digits = 3)

y.ap <- a$coef[1, 1] + a$coef[2, 1] * (x ^ pow.ap)

#draw the lines
points(x, y.ini, typ = 'l')
points(x, y.ap, typ = 'l', col = 'red')

lim <- ceiling(max(c(df$t.ini, df$t.ap)))
legend(x = 200, y = lim - lim / 100, leg = sapply(c('Original', bquote(y == .(i.ini) + .(b.ini)%.%10^.(bpow.ini)%.%x^.(pow.ini)), bquote(R^2 == .(r2.ini)), '',
	'Approximation', bquote(y == .(i.ap) + .(b.ap)%.%10^.(bpow.ap)%.%x^.(pow.ap)), bquote(R^2 == .(r2.ap))), as.expression),
	pch = c(1, NA, NA, NA, 1, NA, NA), lty = c(NA, 1, NA, NA, NA, 1, NA), col = c('black', 'black', 'white', 'white', 'red', 'red', 'white'), bty = 'n')

}

dev.off()
