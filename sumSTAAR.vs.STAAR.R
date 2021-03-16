# based on https://github.com/xihaoli/STAAR/blob/master/docs/STAAR_vignette.html by Xihao Li and Zilin Li

############ 
load('example.rda') # https://github.com/xihaoli/STAAR/blob/master/data/example.rda

library(Matrix)
library(STAAR)

## Loading required package: quadprog
library(MASS)
library(matrixStats) 
library(rje)      # for expit() function

library(sumFREGAT)

genotype <- example$genotype
maf <- example$maf
snploc <- example$snploc

grid <- 1
Npercell <- 10000
ndiv <- 1
N <- round(grid*grid*Npercell/ndiv)

gene.file <- system.file("testfiles/NCBI37.3.geneFile.txt.gz", package = "sumFREGAT")
gf <- read.table(gene.file)[1, ]

res0 <- matrix(NA, 0, 149)
colnames(res0) <- c('iter', paste0('STAAR.B11.A', 0:10), 'STAAR.B11', paste0('STAAR.B125.A', 0:10), 'STAAR.B125', paste0('STAAR.S11.A', 0:10), 'STAAR.S11', paste0('STAAR.S125.A', 0:10), 'STAAR.S125', paste0('STAAR.A11.A', 0:10), 'STAAR.A11', paste0('STAAR.A125.A', 0:10), 'STAAR.A125', 'STAAR.ACAT_O', 'STAAR.STAAR_O', paste0('sumSTAAR.B11.A', 0:10), 'sumSTAAR.B11', paste0('sumSTAAR.B125.A', 0:10), 'sumSTAAR.B125', paste0('sumSTAAR.S11.A', 0:10), 'sumSTAAR.S11', paste0('sumSTAAR.S125.A', 0:10), 'sumSTAAR.S125', paste0('sumSTAAR.A11.A', 0:10), 'sumSTAAR.A11', paste0('sumSTAAR.A125.A', 0:10), 'sumSTAAR.A125', 'sumSTAAR.ACAT_O', 'sumSTAAR.STAAR_O')

#write.table(t(colnames(res0)), file = 'staar.out', col = F, qu = F, row = F)

for (kkk in 1:10) { # 10 simulations

cat(kkk, '')
X1 <- rnorm(N)
X2 <- rbinom(N,1,0.5)
eps <- rnorm(N)

numVar <- dim(snploc)[1]
Z1 <- rnorm(numVar); Z2 <- rnorm(numVar)
Z3 <- rnorm(numVar); Z4 <- rnorm(numVar)
Z5 <- rnorm(numVar); Z6 <- rnorm(numVar)
Z7 <- rnorm(numVar); Z8 <- rnorm(numVar)
Z9 <- rnorm(numVar); Z10 <- rnorm(numVar)
Z <- cbind(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10)

annotation_rank <- colRanks(Z,preserveShape = TRUE) 
M <- dim(annotation_rank)[1]                        ###  M is total number of variants sequenced across the whole genome
PHRED <- -10 * log10(1 - annotation_rank/M)         ### annotation_phred

maxLength <- 1000000                                   # Length of the sequence
sigLength <- 5000                                      # Length of the signal region
startloc <- sample(1 : (maxLength - sigLength + 1),1)  # Location of the signal region
endloc   <- startloc + sigLength - 1

# Extracting the signal region from snploc
snploc <- example$snploc
snplist <- which(snploc$CHROM_POS >= startloc & snploc$CHROM_POS <= endloc)
snpRegion <- snploc[snplist,]

numSNPs <- dim(snpRegion)[1]     # number of variants in the signal region
Geno <- genotype[,snplist]       # genotype of the signal region
numSNPs

# Simulate causal variants from this region
b0 <- rje::logit(0.015)
b  <- rep(log(5),10)

causalprob <- apply(Z[snplist,],1,function(z){
  ind <- sample(1:10,5)
  rje::expit(b0 + b[ind] %*% z[ind])
})
isCausal <- as.logical(rbinom(numSNPs,1,causalprob))

# effects size of covariates 
alpha0 <- 0; alpha1 <- 0.5; alpha2 <- 0.5

# effect size of causal variants
c0 <- 0.13
beta <- -c0 * log10(maf[snplist][which(isCausal==1)])

# generating the phenotype data
Y <- alpha0 + alpha1 * X1 + alpha2 * X2 + as.matrix(genotype)[,snplist][,which(isCausal==1)] %*% beta + eps

### Fitting null model for unrelated samples using fit_null_glm()  ###
pheno         <- data.frame(Y = Y, X1 = X1, X2 = X2)
obj_nullmodel <- fit_null_glm(Y ~ X1 + X2,data = pheno,family = "gaussian")

rare_maf_cutoff  <- 1                  ###0.05
maf0             <- maf[snplist]
annotation_phred <- as.matrix(PHRED[snplist,])

RV_label         <- as.vector((maf0 < rare_maf_cutoff) & (maf0 > 0))
annotation_phred <- annotation_phred[RV_label,,drop=FALSE]
colnames(annotation_phred) <- paste0('PROB', 1:10)
#annotation_phred [3, 5] <- Inf

Geno <- as.matrix(Geno[,RV_label])
maf0 <- maf0[RV_label]

# STAAR p-values
pvalues <- STAAR(genotype=Geno,obj_nullmodel=obj_nullmodel,annotation_phred=annotation_phred,rare_maf_cutoff=1)

# calculating the summary statistics
NNN  <- length(Y); nnn <- sqrt(NNN)
X    <- cbind(1,X1,X2)
Pr   <- diag(rep(1,NNN)) - X %*% (solve(t(X) %*% X)) %*% t(X)        ### projection matrix
CenY <- as.vector(Pr %*% Y)
CenG <- as.matrix(Pr %*% Geno)
sdY  <- sqrt(mean(CenY^2))  ### sd(CenY)
sdG  <- as.vector(sqrt(diag(t(CenG) %*% CenG)/NNN))  ### #sdG  <- SD(CenG)

sebetas <- sdY/(sdG * nnn)
betas   <- as.vector((t(CenG) %*% CenY)) / as.vector(diag(t(CenG) %*% CenG))
zscores  <- betas/sebetas   ### zscores2 <- as.vector((t(CenG) %*% CenY)/(nnn * sdG * sdY))

U <- diag(1/sdG) %*% (t(CenG) %*% CenG) %*% diag(1/sdG)/NNN  ### corr matrix among SNPs  ### U1 <- cov2cor(cov(CenG))
colnames(U) <- rownames(U) <- snpRegion$SNP
save(U, file = 'gene.RDa')

ref <- data.frame(ID = colnames(U), REF = 0, ALT = 1)
input.data <- cbind(data.frame(CHROM = snpRegion$CHROM, POS = snpRegion$CHROM_POS, ID = snpRegion$SNP, EA = 1, P = (pnorm(abs(zscores), lower.tail = FALSE) * 2), BETA = betas, EAF = sapply(1:dim(Geno)[2], function(x)sum(Geno[,x])/20000)), annotation_phred)

prep.score.files(input.data, ref = ref) # write score file 'scores.vcf.gz'

# write dummy genefile
gf[1, 1] <- 'gene'
gf[1, 5] <- min(snpRegion$CHROM_POS) - 1
gf[1, 6] <- max(snpRegion$CHROM_POS) + 1
write.table(gf, row = F, qu = F, file = 'genefile.txt')

## sumSTAAR p-values
res <- sumSTAAR(score.file = 'scores.vcf.gz', gene.file = 'genefile.txt', genes = 'gene', cor.path = '', tests = c('BT', 'SKAT', 'ACAT'), prob = paste0('PROB', 1:10), beta.par.matrix <- rbind(c(1, 1), c(1, 25)), quiet = TRUE)

res$gene <- NULL
res <- c(kkk, c(as.numeric(pvalues$results_STAAR_B_1_1), as.numeric(pvalues$results_STAAR_B_1_25), as.numeric(pvalues$results_STAAR_S_1_1), as.numeric(pvalues$results_STAAR_S_1_25), as.numeric(pvalues$results_STAAR_A_1_1), as.numeric(pvalues$results_STAAR_A_1_25)), pvalues$results_ACAT_O, pvalues$results_STAAR_O, res)

#write.table(t(res), file = 'staar.out', col = F, qu = F, row = F, app = T)

res0 <- rbind(res0, res)

}

##### plot

pdf('STAAR_vs_sumFREGAT.pdf', height = 8)

par(mfrow = c(2, 2))

x <- -log10(as.numeric(res0[, c(paste0('STAAR.B11.A', 0:10), paste0('STAAR.B125.A', 0:10))]))
y <- -log10(as.numeric(res0[, c(paste0('sumSTAAR.B11.A', 0:10), paste0('sumSTAAR.B125.A', 0:10))]))
lim <- ceiling(max(c(x, y)))
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'Burden tests')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')

x <- -log10(as.numeric(res0[, c(paste0('STAAR.S11.A', 0:10), paste0('STAAR.S125.A', 0:10))]))
y <- -log10(as.numeric(res0[, c(paste0('sumSTAAR.S11.A', 0:10), paste0('sumSTAAR.S125.A', 0:10))]))
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'SKAT')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')

x <- -log10(as.numeric(res0[, c(paste0('STAAR.A11.A', 0:10), paste0('STAAR.A125.A', 0:10))]))
y <- -log10(as.numeric(res0[, c(paste0('sumSTAAR.A11.A', 0:10), paste0('sumSTAAR.A125.A', 0:10))]))
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'ACAT')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')

x <- -log10(as.numeric(res0[, c('STAAR.B11', 'STAAR.B125', 'STAAR.S11', 'STAAR.S125', 'STAAR.A11', 'STAAR.A125', 'STAAR.ACAT_O', 'STAAR.STAAR_O')]))
y <- -log10(as.numeric(res0[, c('sumSTAAR.B11', 'sumSTAAR.B125', 'sumSTAAR.S11', 'sumSTAAR.S125', 'sumSTAAR.A11', 'sumSTAAR.A125', 'sumSTAAR.ACAT_O', 'sumSTAAR.STAAR_O')]))
plot(x, y, ylim = c(0, lim), xlim = c(0, lim), xlab = 'STAAR', ylab = 'sumFREGAT', main = 'Combined tests')
abline(0, 1)
abline(reg = lm(y ~ x), col = 'red')

dev.off()

