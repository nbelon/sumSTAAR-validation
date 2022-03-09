### a script used to generate Figures S2, S3
# Nadezhda Belonogova, Irina Zorkoltseva 2022

##### STAAR null model

library(Matrix)
library(STAAR)
library(data.table)

phen <- fread('/home/lima/Zor/STAAR/CAD/I25cov.nonrel.ord.csv', h = T, data.table = F)
# UK Biobank phenotype of the Chronic ischaemic heart disease (IHD, ICD-10 code I25)
# containing 153,379 unrelated individuals (12,931 cases / 140,448 controls) with European ancestry (project #59345)

obj_nullmodel <- fit_null_glm(I25 ~ Sex + age + batch, data = phen,family = "binomial")
save(obj_nullmodel,file='obj_nullmodel.Rda') 

####### STAAR individual scores 

library(Matrix)
library(STAAR)
library(snpStats)
library(data.table)

load('obj_nullmodel.Rda')

mat = matrix(ncol = 0, nrow = 0)
ind.pval <- data.frame(mat)

M=756200
nvar = 10000
ndose =M%/%nvar

for (i in 1:ndose) {
start <-(i-1)*nvar + 1
end <- start + nvar - 1
Geno <- read.plink('/home/lima/Zor/STAAR/CAD/white_nonrel_chr01',select.snps = c(start:end))
gt <- as(Geno$genotype,'numeric')
pval <- Indiv_Score_Test_Region(genotype=gt,obj_nullmodel,rv_num_cutoff = 5,rare_maf_cutoff = 1)
ind.pval <- rbind (ind.pval,pval)
write.table(pval, file = 'pval.all.staar.chr01.tmp', row.names = F, qu = F, sep=' ', append = T, col.names = F)
print(paste('i=',i,'from ',start,' to ',end))
}

##Last dose
start = nvar*ndose+1
Geno <- read.plink('/home/lima/Zor/STAAR/CAD/white_nonrel_chr01',select.snps = c(start:M))
gt <- as(Geno$genotype,'numeric')
pval <- Indiv_Score_Test_Region(genotype=gt,obj_nullmodel,rv_num_cutoff = 5,rare_maf_cutoff = 1)
ind.pval <- rbind (ind.pval,pval)
print(paste('end of last dose', start,M))
write.table(pval, file = 'pval.all.staar.chr01.tmp', row.names = F, qu = F, sep=' ', append = T, col.names = F)

print('read map file')
map <- fread('/home/lima/Zor/STAAR/CAD/white_nonrel_chr01.bim', data = F, h = F)
map <- map[,-c(3)]
colnames (map)[1:5] <- c('CHR','SNP','POS','A1','A2')
total <- cbind(map,ind.pval)
print('write total result to file ind.staar.chr01.txt')
write.table(total, file = 'ind.pval.all.staar.chr01.txt', row.names = F, qu = F,sep=' ')

####### Figure S3: fastGWA vs STAAR individual scores

staar <- read.table('/home/lima/Zor/STAAR/staar/chr01/ind.pval.all.staar.chr01.txt', h = T)
fast <- read.table('/home/lima/Zor/STAAR/staar/chr01/fast4gb_mac5.maf_0.01.txt', h = T)
# CAD summary statistics from fastGWA, calculated from the UK Biobank whole exome sequencing data and phenotype of the Chronic ischaemic heart disease (IHD, ICD-10 code I25)
# containing 153,379 unrelated individuals (12,931 cases / 140,448 controls) with European ancestry (project #59345)

fast <- fast[fast$N == 153375, ] # genotypes all measured
fast$zt <- fast$T/fast$SE_T
fast$p.fastT <- 2*pnorm(abs(fast$zt),low=F)
staar$zt <- staar$Score/staar$SE
staar$p.staar <- 2*pnorm(abs(staar$zt),low=F)

df <- merge(staar, fast, by = 'SNP')
#105310 SNPs

#pdf('Fig.S3.pdf')
tiff('Fig.S3.tif', units="in", width=5, height=5.5, res=600, compression = 'lzw')

x <- -log10(df[, 'p.staar'])
y <- -log10(df[, 'p.fastT'])
cor <- cor(x, y)
lim <- ceiling(max(c(x, y)))
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'fastGWA-GLMM')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')
legend('topleft', legend = paste0('r = ', format(round(cor, 6), nsmall = 3)), pch = '', bty = "n")

dev.off()



#### prepare summary stats for sumSTAAR

fast <- read.table('fast4gb_mac5.maf_0.01.txt', h = T)

fast <- fast[fast$N == 153375, ] # genotypes all measured

fast$zt <- fast$T/fast$SE_T
fast$pz <- 2*pnorm(abs(fast$zt),low=F)

gwas <- fast[, c('CHR', 'SNP', 'POS', 'A1', 'AF1', 'BETA', 'SE', 'pz')]###fast
colnames(gwas)[1:8] <- c('CHR','ID','POS','EA','EAF','BETA','SE','P')###fast

write.table(gwas, file = 'fastT.all.measured.txt', row = F, qu = F)

ref <- '/home/common/DataStorage/UKBB/Project_18219/Matrix/RData/reference.RData' # available at https://mga.bionet.nsc.ru/ukbb_exome_matrix/
library(sumFREGAT)
prep.score.files('fastT.all.measured.txt', reference = ref, output.file.prefix = 'fastT.all.measured') # generates 'fastT.all.measured.vcf.gz'


############# sumSTAAR

cor.path <- '/home/common/DataStorage/UKBB/Project_18219/Matrix/RData/' # matrices available at https://mga.bionet.nsc.ru/ukbb_exome_matrix/
scoreFile <- 'fastT.all.measured.vcf.gz'

library(sumFREGAT)

res <- sumSTAAR(score.file = scoreFile, gene.file = 'hg38', genes = 'chr1', cor.path = cor.path, tests = c('ACAT', 'BT', 'SKAT'), beta = rbind(c(1, 1), c(1, 25)), mac = 10, n = 153379, prob.causal = NA, quiet = TRUE, write.file = F)
save(res, file = 'all.meas.RDa')


###### STAAR

library(Matrix)
library(STAAR)
library(snpStats)
library(data.table)

load('obj_nullmodel.Rda')

geneFile <- "/home/lima/nadya/FREGAT/sumFREGAT_1.2.3/inst/testfiles/ensembl.hg38.txt" 
gens <- read.table(geneFile)
gens <- gens[gens$V3 =='chr1',c(1,5,6)]
colnames(gens)[1:3] <- c('gene','start','end')

# take SNPs with call rate = 1
list <- fread('~/Zor/STAAR/staar/chr01/snplist4gb_mac5.maf_0.01.txt',h = T, data.table = F)
fast <- read.table('fastT.all.measured.txt', h = T)
list <- list[list$SNP %in% fast$ID, ]
dim(list)

mat = matrix(ncol = 0, nrow = 0)
total <- data.frame(mat)

for (i in 1:dim(gens)[1]) {

	start <- gens$start[i]
	end <- gens$end[i]
	g <- gens$gene[i]
	snps <- list$SNP[which(list$POS >= start & list$POS <= end)]

	if (length(snps)> 1){
		print(paste(i, g, length(snps)))
		Geno <- read.plink('/home/lima/Zor/STAAR/CAD/white_nonrel_chr01',select.snps = snps)
		gt <- as(Geno$genotype,'numeric')
		pval <- STAAR(genotype=gt,obj_nullmodel,annotation_phred=NULL,rare_maf_cutoff = 1)
		pval_gen <- c()
		pval_gen$gene <- g
		pval_gen$nvar <- pval$num_variant
		pval_gen$staar_o.pv <- pval$results_STAAR_O
		pval_gen$acat_o.pv <- pval$results_ACAT_O
		pval_gen$skat_1_1.pv <- as.numeric(pval$results_STAAR_S_1_1[1])
		pval_gen$skat_1_25.pv <- as.numeric(pval$results_STAAR_S_1_25[1])
		pval_gen$bt_1_1.pv <- as.numeric(pval$results_STAAR_B_1_1[1])
		pval_gen$bt_1_25.pv <- as.numeric(pval$results_STAAR_B_1_25[1])
		pval_gen$acat_1_1.pv <- as.numeric(pval$results_STAAR_A_1_1[1])
		pval_gen$acat_1_25.pv <- as.numeric(pval$results_STAAR_A_1_25[1])
		write.table(pval_gen, file = 'gb.staar.allmeas.out', row.names = F,qu = F,sep=' ',append=T,col.names = F) 
	}
}


###### draw Figure S4: STAAR vs sumSTAAR on real data

df1 <- read.table('gb.staar.allmeas.out')
colnames(df1) <- c('gene', 'nvar', 'STAAR.STAAR_O', 'STAAR.ACAT_O', 'STAAR.S11', 'STAAR.S125', 'STAAR.B11', 'STAAR.B125', 'STAAR.A11', 'STAAR.A125')
load('all.meas.RDa')
colnames(res) <- c('gene', 'sumSTAAR.A11', 'sumSTAAR.A125', 'sumSTAAR.B11', 'sumSTAAR.B125', 'sumSTAAR.S11', 'sumSTAAR.S125', 'sumSTAAR.ACAT_O')
res <- res[which(!is.na(res$sumSTAAR.B11)), ]

res0 <- merge (df1,res,by='gene')
dim(res0)
#[1] 1927   17
sum(res0$nvar)
# 110538
res0$gene <- NULL
res0 <- as.matrix(res0)

#pdf('Fig.S4.pdf', height = 7.5)
tiff('Fig.S4.tif', units="in", width=7, height=7.5, res=600, compression = 'lzw')

par(mfrow = c(2, 2))

x <- -log10(as.numeric(res0[, c('STAAR.B11', 'STAAR.B125')]))
y <- -log10(as.numeric(res0[, c('sumSTAAR.B11', 'sumSTAAR.B125')]))
cor <- cor(x, y)
lim <- ceiling(max(c(x, y)))
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'Burden tests')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')
legend('topleft', legend = paste0('r = ', format(round(cor, 6), nsmall = 3)), pch = '', bty = "n")

x <- -log10(as.numeric(res0[, c('STAAR.S11', 'STAAR.S125')]))
y <- -log10(as.numeric(res0[, c('sumSTAAR.S11', 'sumSTAAR.S125')]))
cor <- cor(x, y)
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'SKAT')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')
legend('topleft', legend = paste0('r = ', format(round(cor, 6), nsmall = 3)), pch = '', bty = "n")

x <- -log10(as.numeric(res0[, c('STAAR.A11', 'STAAR.A125')]))
y <- -log10(as.numeric(res0[, c('sumSTAAR.A11', 'sumSTAAR.A125')]))
cor <- cor(x, y)
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'ACAT')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')
legend('topleft', legend = paste0('r = ', format(round(cor, 6), nsmall = 3)), pch = '', bty = "n")

x <- -log10(as.numeric(res0[, 'STAAR.ACAT_O']))
y <- -log10(as.numeric(res0[, 'sumSTAAR.ACAT_O']))
cor <- cor(x, y)
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'Combined tests')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')
legend('topleft', legend = paste0('r = ', format(round(cor, 6), nsmall = 3)), pch = '', bty = "n")

dev.off()


######## ACAT on STAAR pvals

################ score file for sumFREGAT

fast <- read.table('fast4gb_mac5.maf_0.01.txt', h = T)
fast <- fast[fast$N == 153375, ] # genotypes all measured
staar <- read.table('ind.pval.all.staar.chr01.txt', h = T) # individual scores from STAAR

v <- match(fast$SNP, staar$SNP)
fast$pz <- staar$pval[v]

fast <- fast[fast$SNP %in% staar$SNP, ]

#gwas <- fast[,c(1:4,7,11,12,17)]###fast
gwas <- fast[, c('CHR', 'SNP', 'POS', 'A1', 'AF1', 'BETA', 'SE', 'pz')]###fast
colnames(gwas)[1:8] <- c('CHR','ID','POS','EA','EAF','BETA','SE','P')###fast

write.table(gwas, file = 'staar.p.txt', row = F, qu = F)

ref <- '/home/common/DataStorage/UKBB/Project_18219/Matrix/RData/reference.RData'
library(sumFREGAT)
prep.score.files('staar.p.txt', reference = ref, output.file.prefix='staar.p') # generates 'staar.p.vcf.gz'

######### ACAT

library(sumFREGAT)

scoreFile <- 'staar.p.vcf.gz'

flout <-'ACAT.staarp.out'
res <- sumSTAAR(score.file = scoreFile, gene.file = 'hg38', genes = 'chr1', cor.path = cor.path, tests = c('ACAT'), beta = rbind(c(1, 1), c(1, 25)), mac = 10, n = 153379, prob.causal = NA, quiet = TRUE, write.file = flout)

########## compare

df1 <- read.table('gb.staar.allmeas.out') # STAAR gene-based tests output
colnames(df1) <- c('gene', 'nvar', 'STAAR.STAAR_O', 'STAAR.ACAT_O', 'STAAR.S11', 'STAAR.S125', 'STAAR.B11', 'STAAR.B125', 'STAAR.A11', 'STAAR.A125')
df2 <- read.table('ACAT.staarp.out', h = T)
df <- merge (df1,df2,by='gene')

p1<- -log10(c(df$STAAR.A11, df$STAAR.A125))
p2 <- -log10(c(df$ACAT.1.1.PROB0, df$ACAT.1.25.PROB0))
cor <- cor(p1[!is.na(p2)],p2[!is.na(p2)])
cor
#[1] 0.9985428
