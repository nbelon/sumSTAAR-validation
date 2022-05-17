### a script used to generate the data for Table S1 and Figure S5
# Nadezhda Belonogova, 2022

library(sumFREGAT)
library(data.table)

df <- fread('../../sumstats_neuro_sum_ctg_format.txt.gz', h = T, data.table = F) # neuroticism summary statistics from https://ctg.cncr.nl/software/summary_statistics
df <- df[, c('CHR', 'POS', 'SNP', 'A1', 'P', 'BETA', 'EAF')]
colnames(df)[3:4] <- c('ID', 'EA')
fwrite(df, '4prep.txt', row = F, qu = F, sep = ' ', na = 'NA')

####### link annotations from 1KG 
library(data.table)
df <- fread('4prep.txt', h = T, data.table = F)
df$ind <- paste(df$CHR, df$POS, sep = ':')

for (ch in 1:22) {

cat(ch, '')
anno <- read.table(paste0('ALL.chr', ch, '.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz'), as.is = T)
anno$ind <- paste(anno$V1, anno$V2, sep = ':')
v <- match(anno$ind, df$ind)
anno$my <- as.character(df$EA[v])

tmp <- strsplit(anno$V8, ';')
tmp <- unlist(tmp)
csq0 <- tmp[grepl('CSQ', tmp)]
csq <- strsplit(csq0, '|', fixed = T)

l <- status <- rep(NA, length(csq))

for (i in 1:length(csq)) {
	l[i] <- length(csq[[i]])
	if (l[i] > 19 & !is.na(anno$my[i])) {
		tmp <- strsplit(csq0[i], '=', fixed = T)[[1]][2]
		tmp <- unlist(strsplit(tmp, ',', fixed = T))
		my.allele <- anno$my[i]
		for (j in 1:length(tmp)) {
			tmp1 <- unlist(strsplit(tmp[j], '|', fixed = T))
			if (tmp1[1] == my.allele) {
				status[i] <- tmp1[5]
			}
		}
	} else {
		status[i] <- csq[[i]][5]
	}
}

anno$status <- status
df$ANNO[df$CHR == ch] <- as.character(anno$status[match(df$ind[df$CHR == ch], anno$ind)])
}

fwrite(df, '4prep.anno.txt', row = F, qu = F, sep = ' ', na = 'NA')


##### select exonic SNPs

library(data.table)
df <- fread('4prep.anno.txt', h = T, data.table = F)
df <- df[!is.na(df$ANNO), ]
df$ANNO0 <- df$ANNO
noncoding <- c('3_prime_UTR_variant', '5_prime_UTR_variant', 'downstream_gene_variant', 'upstream_gene_variant', 'intron_variant', 'regulatory_region_variant', 'TF_binding_site_variant') # remove these
df$ANNO <- !grepl(paste(noncoding, collapse = "|"), df$ANNO)
df <- df[df$ANNO, ]
df$ANNO <- df$ANNO0 <- NULL
write.table(df, file = '4prep.ACAT.exon-only.txt', row = F, qu = F)


###### polygene pruning for exonic SNPs

library(data.table)

# read genefile
geneFile <- "/home/lima/nadya/FREGAT/sumFREGAT_1.0.1/inst/testfiles/NCBI37.3.geneFile.txt" 
df <- read.table(geneFile)
df$V2 <- df$V4 <- NULL
colnames(df)[1:4] <- c('gene', 'chr', 'start', 'stop')
df$chr <- gsub('chr', '', df$chr)
df$gene <- as.character(df$gene)

tp <- 'r2.0.2'
load(paste0('../res/clist.', tp, '.RDa')) # file from clumping.R
ran <- read.table(paste0('../res/r2/', tp, '.clumped.ranges'), h = T) # file from plink
ran$SNP <- as.character(ran$SNP)

pref1 <- '.exon'
pref2 <- paste0(pref1, '-only')
pref3 <- paste0('.ACAT', pref2)
gwas <- fread(paste0('../4prep', pref3, '.txt'), h = T, data.table = F)# variable file name, specific for a scenario 
# keep only SNPs present in 1KG
map1kg <- fread('/home/common/g1000p3_EUR/origin/g1000p3_EUR.bim', data.table = F) # 1KG snp list
v <- gwas$ID %in% map1kg$V2
table(v)
gwas <- gwas[v, ]
# remove SNPs monomorphic in 1KG
freq <- fread('../res/r2/g1000p3_EUR.frq', data.table = F, h = T)
freq <- freq[freq$MAF > 0, ]
v <- gwas$ID %in% freq$SNP
table(v)
gwas <- gwas[v, ]
#ran has monomorphic SNPs with no clumps, that is ok
#save(gwas, file = '../gwas.1kg.RDa')
save(gwas, file = '../gwas.1kg.exon.RDa')

# 'acat' is a full stats data table
acat <- fread(paste0('../4prep', pref3, '.txt'), h = T, data.table = F)
prep0 <- fread(paste0('../4prep.ukbb', pref2, '.txt'), h = T, data.table = F) # SNP ids as in UKBB matrices
# rs is a common ind, all ukbb rs are in acat ID
v <- match(prep0$rs, acat$ID)
table(acat$P[v] == prep0$P)
#    TRUE
# 4395683
acat$ini.ID <- acat$ID
acat[v, 1:7] <- prep0[, 1:7]
#save(acat, file = '../acat.combinedID.RDa')
save(acat, file = '../acat.combinedID.exon.RDa')

#load('../gwas.1kg.exon.RDa') # 'gwas'
#load('../acat.combinedID.exon.RDa') # 'acat'

# pruning
write.table(t(c('gene', 'snp')), file = 'all.genes.snps2keep.txt', row = F, qu = F, col = F)

for (i in 1:dim(df)[1]) {
	g <- df[i, 'gene']
	cat(g, '')
	topsnps <- ran[grep(paste0('\\b', g, '\\b'), ran[, 'RANGES']), 'SNP']
	vvv <- gwas$CHR == df$chr[i] & gwas$POS >= df$start[i] & gwas$POS <= df$stop[i]
	genesnps <- gwas[vvv, 'ID']
	if (length(genesnps) == 0) next
	topsnps <- topsnps[!topsnps %in% genesnps] # outside topsnps only
	snps2waste <- c()
	for (topsnp in topsnps) {
		snps2waste <- c(snps2waste, clist[[topsnp]]) # correlated snps from clumps
	}
	if (sum(genesnps %in% snps2waste) == 0) { # take all snps regardless 1kg
		vvv <- acat$CHR == df$chr[i] & acat$POS >= df$start[i] & acat$POS <= df$stop[i]
		genesnps <- acat[vvv, 'ini.ID']
	} else { # if there are snps to exclude
		cat('!!! ')
		write.table(g, file = 'all.genes.pruned.exon.txt', row = F, qu = F, col = F, app = T)
		v <- !genesnps %in% snps2waste
		genesnps <- genesnps[v]
		if (sum(v) == 0) { # all snps excluded
			cat('no snps ')
			next
		}
		write.table(cbind(g, genesnps), file = 'all.genes.snps2keep.exon.txt', row = F, qu = F, col = F, app = T)
	}

	######## write new files with pruned snps only
	prep.acat <- acat[acat$ini.ID %in% genesnps, ]
	write.table(prep.acat, file = paste0('all.genes.pruned/exon/', g, '.pruned.txt'), row = F, qu = F)
	write.table(g, file = 'all.genes.exon.txt', row = F, qu = F, col = F, app = T)

}



############# write rsIDs to get annotations

geneFile <- "/home/lima/nadya/FREGAT/sumFREGAT_1.0.1/inst/testfiles/NCBI37.3.geneFile.txt" 
df <- read.table(geneFile)
df$V2 <- df$V4 <- NULL
colnames(df)[1:4] <- c('gene', 'chr', 'start', 'stop')
df$chr <- gsub('chr', '', df$chr)
df$gene <- as.character(df$gene)

m <- 0
n <- 1

for (i in 1:dim(df)[1]) {

	g <- df[i, 'gene']
	fn <- paste0('all.genes.pruned/exon/', g, '.pruned.txt')
	if (!file.exists(fn)) next
	rs <- read.table(fn, h = T)$ini.ID
	m <- m + length(rs)
	if (m > 10000) {
		m <- length(rs)
		n <- n + 1
	}
	write.table(rs, file = paste0('snps.', n, '.txt'), col = F, row = F, qu = F, app = T)
}


######### liftover for SNPs not found by rsIDs

files <- list.files(pattern = "\\.csv.gz$")
write.table(t(colnames(an)), file = 'favor.exon.txt', col = F, row = F, qu = F)

for (n in 1:11) {
	cat(n, '')
	rs <- unique(read.table(paste0('snps.', n, '.txt'))$V1)
	an <- read.table(files[n], sep = ',', h = T)
	rs.not <- rs[!rs %in% an$rsID]
	write.table(rs.not, file = 'rs.not.found.txt', col = F, row = F, qu = F, app = T)
	write.table(an, file = 'favor.exon.txt', col = F, row = F, qu = F, app = T)
}

tmp <- unique(read.table('rs.not.found.txt')$V1)
tmp <- paste(tmp, collapse = '", "')
tmp <- paste0('"', tmp, '"')
system(paste0("mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg38 -e 'select chrom,chromStart,chromEnd,name from snp147 where name in (", tmp, ")' > snps.not.found.hg38.txt"))

tmp <- read.table('snps.not.found.hg38.txt', h = T)
hg38 <- paste(tmp$chrom, tmp$chromEnd, sep = '-')
write.table(hg38, file = 'hg38.txt', col = F, row = F, qu = F)


##### link annotations downloaded from favor.genohub.org/

library(data.table)

anno <- c("Variant", "Chromosome", "Position", "TOPMed.QC.Status", "rsID", "LINSIGHT", "FATHMM.XF.coding", "FATHMM.XF.noncoding", "CADD.RawScore", "CADD.PHRED", "aPC.Epigenetics", "aPC.Conservation", "aPC.Protein.Function", "aPC.Local.Nucleotide.Diversity", "aPC.Mutation.Density", "aPC.Transcription.Factor", "aPC.Mappability", "aPC.Proximity.To.TSS.TES", "Aloft.Value", "Aloft.Description", "Funseq.Value", "Funseq.Description")
files <- list.files(pattern = "\\.csv.gz$") # downloaded from favor.genohub.org/

write.table(t(colnames(an)[colnames(an) %in% anno]), file = 'favor.exon.txt', col = F, row = F, qu = F, sep = ',')
for (n in 1:11) {
	cat(n, '')
	fn <- files[n]
	if (fn == '2022-01-21 12 12 00 _annotations.csv.gz') next # liftover SNPs
	an <- read.table(fn, sep = ',', h = T)
	an <- an[, colnames(an) %in% anno]
	write.table(an, file = 'favor.exon.txt', col = F, row = F, qu = F, app = T, sep = ',')
}

an <- read.table('2022-01-21 12 12 00 _annotations.csv.gz', sep = ',', h = T) # liftover SNPs
an <- an[, colnames(an) %in% anno]
hg38 <- read.table('snps.not.found.hg38.txt', h = T)

hg38$ind <- paste(hg38$chrom, hg38$chromEnd, sep = '-')
an$ind <- paste0('chr', an$Chromosome, '-', an$Position)
an$rsID <- hg38$name[match(an$ind, hg38$ind)]
an$ind <- NULL

an0 <- read.table('favor.exon.txt', h = T, sep = ',')
an0 <- rbind(an0, an)
write.table(an0, file = 'favor.exon.txt', row = F, qu = F, sep = ',')

# check alt ref
df0 <- fread('../../sumstats_neuro_sum_ctg_format.txt.gz', h = T, data.table = F) # original summary statistics for neuroticism
df0 <- df0[df0$SNP %in% an0$rsID, ]

tmp <- unlist(strsplit(an0$Variant, '-'))
an0$A1 <- tmp[(1:dim(an0)[1])*4-1]
an0$A2 <- tmp[(1:dim(an0)[1])*4]
an0$ind1 <- paste(an0$rsID, an0$A1, an0$A2, sep = ':')
an0$ind2 <- paste(an0$rsID, an0$A2, an0$A1, sep = ':')
df0$ind <- paste(df0$SNP, df0$A1, df0$A2, sep = ':')
v1 <- an0$ind1 %in% df0$ind
v2 <- an0$ind2 %in% df0$ind
v <- v1 | v2
an1 <- an0[v, ]

# alternative chain
snps <- df0[!(df0$SNP %in% an1$rsID), 'SNP']
an2 <- an0[an0$rsID %in% snps, ]
an2$A1[an2$A1 == 'G'] <- 'O'
an2$A1[an2$A1 == 'C'] <- 'G'
an2$A1[an2$A1 == 'O'] <- 'C'
an2$A1[an2$A1 == 'A'] <- 'O'
an2$A1[an2$A1 == 'T'] <- 'A'
an2$A1[an2$A1 == 'O'] <- 'T'
an2$A2[an2$A2 == 'G'] <- 'O'
an2$A2[an2$A2 == 'C'] <- 'G'
an2$A2[an2$A2 == 'O'] <- 'C'
an2$A2[an2$A2 == 'A'] <- 'O'
an2$A2[an2$A2 == 'T'] <- 'A'
an2$A2[an2$A2 == 'O'] <- 'T'
an2$ind1 <- paste(an2$rsID, an2$A1, an2$A2, sep = ':')
an2$ind2 <- paste(an2$rsID, an2$A2, an2$A1, sep = ':')
v1 <- an2$ind1 %in% df0$ind
v2 <- an2$ind2 %in% df0$ind
v <- v1 | v2
an1 <- rbind(an1, an2[v, ])
an1 <- an1[-which(duplicated(an1$rsID)), ]
write.table(an1, file = 'favor.exon.alleles.txt', row = F, qu = F, sep = ',')



######### write score files with annotations for sumSTAAR

library(sumFREGAT)

genes <- read.table('all.genes.exon.txt')$V1 # genes list
an <- read.table('favor.exon.alleles.txt', h = T, sep = ',') # annotations
anno <- c("FATHMM.XF.coding", "CADD.PHRED", "aPC.Epigenetics", "aPC.Conservation", "aPC.Protein.Function", "aPC.Local.Nucleotide.Diversity", "aPC.Mutation.Density", "aPC.Transcription.Factor", "aPC.Mappability", "aPC.Proximity.To.TSS.TES")
v <- colnames(an) %in% anno
load('/home/belon/ya/ref.ldstore265.all.RData') # 'ref', available at https://mga.bionet.nsc.ru/sumfregat/ukbb/

for (g in genes) {

	df <- read.table(paste0('all.genes.pruned/exon/', g, '.pruned.txt'), h = T)
	df <- cbind(df, an[match(df$ini.ID, an$rsID), v])
	colnames(df)[10:dim(df)[2]] <- paste0('PROB', 1:10, '.', colnames(df)[10:dim(df)[2]])
	write.table(df, file = 'tmp.txt', row = F, qu = F)

	prep.score.files('tmp.txt', reference = ref, output.file.prefix = paste0('all.genes.pruned/exon/', g)) # generates paste0(g, '.vcf.gz') score files

}


######## run sumSTAAR

library(sumFREGAT)
genes <- read.table('all.genes.exon.txt')$V1 # genes list
anno <- c("FATHMM.XF.coding", "CADD.PHRED", "aPC.Epigenetics", "aPC.Conservation", "aPC.Protein.Function", "aPC.Local.Nucleotide.Diversity", "aPC.Mutation.Density", "aPC.Transcription.Factor", "aPC.Mappability", "aPC.Proximity.To.TSS.TES")
corpath <- '265K.all.RData/' # matrices available at https://mga.bionet.nsc.ru/sumfregat/ukbb/

for (g in genes[1:length(genes)]) {

	cat(g, '')
	tmp <- read.table(paste0('all.genes.pruned/exon/', g, '.vcf.gz'))
	if (length(grep('SE', tmp[1, ])) == 0) next
	reg <- sumSTAAR(score.file=paste0(g, '.vcf.gz'), genes = g, tests = c('ACAT', 'SKAT', 'BT', 'PCA'), cor.path = corpath, gene.file = 'hg19', phred = TRUE, prob.causal= paste0('PROB', 3:10, '.', anno[3:10]), n = 380506, write.file = paste0('all.genes.pruned/exon/', g, '.out'), quiet = T)

}


######### read the results

genes <- read.table('all.genes.exon.txt')$V1

df <- c()
i = 0

for (g in genes) {

	i = i + 1
	if (!i %% 1000) cat (i, '')
	fn <- paste0(g, '.out')
	if (file.exists(fn)){
		df <- rbind(df, read.table(fn, h = T))
	}

}
#save(df, file = 'favor.exon.RDa')


##### combine

library(sumFREGAT)

df$anno11 <- sapply(1:dim(df)[1], function(x) ACATO(c(df$PCA.1.1.STAAR[x], df$SKAT.1.1.STAAR[x], df$BT.1.1.STAAR[x], df$ACAT.1.1.STAAR[x])))
df$anno125 <- sapply(1:dim(df)[1], function(x) ACATO(c(df$PCA.1.25.STAAR[x], df$SKAT.1.25.STAAR[x], df$BT.1.25.STAAR[x], df$ACAT.1.25.STAAR[x])))

df$PCAanno <- sapply(1:dim(df)[1], function(x) ACATO(c(df$PCA.1.1.STAAR[x], df$PCA.1.25.STAAR[x])))
df$BTanno <- sapply(1:dim(df)[1], function(x) ACATO(c(df$BT.1.1.STAAR[x], df$BT.1.25.STAAR[x])))
df$SKATanno <- sapply(1:dim(df)[1], function(x) ACATO(c(df$SKAT.1.1.STAAR[x], df$SKAT.1.25.STAAR[x])))
df$ACATanno <- sapply(1:dim(df)[1], function(x) ACATO(c(df$ACAT.1.1.STAAR[x], df$ACAT.1.25.STAAR[x])))

out <- df[, c('gene', 'sumSTAAR.STAAR_O', 'sumSTAAR.ACAT_O', 'PCAanno', 'SKATanno', 'BTanno', 'ACATanno', 'anno11', 'anno125')]

v <- sapply(1:dim(out)[1], function(x) any(out[x, 2:dim(out)[2]] <= 2.5e-5, na.rm = T))
write.table(out[v, ], file = 'favor2.exon.2.5e-5.csv', sep = ';', row = F, qu = F) # data for Table S1


################ Figure S5: qqplots

ylb <- expression(paste(-log[10], "(observed P)"))
xlb <- expression(paste(-log[10], "(expected P)"))

draw.qq <- function(pval, main) {

	chi <- qchisq((1 - pval), 1)
	lam <- median(chi)/qchisq(0.5, 1)

	pval <- -log10(pval)

	n <- length(pval)
	a <- 1:n
	b <- a/n
	x <- -log10(b)

	upper <- qbeta(0.025, a, rev(a))
	lower <- qbeta(0.975, a, rev(a))

	plot(x, pval, type = "n", ylab = ylb, xlab = xlb, ylim = c(0, ceiling(max(pval))), main = main)
	# upper and lower have already been subset
	polygon(-log10(c(b, rev(b))), -log10(c(upper, rev(lower))), density=NA, col="gray80")
	points(x, pval, pch = 19, col = 'black')
	abline(0, 1)
	legend('topleft', legend = paste0('Lambda = ', round(lam, 3)), pch = 19, bty = "n")#, cex = .9)

}

tiff('Fig.S5.tif', units = "in", width = 5.5, height= 11, res = 600, compression = 'lzw')

par(mfrow = c(4, 2))

draw.qq(sort(out$sumSTAAR.STAAR_O), 'sumSTAAR, annotation included') ## sort() removes NAs
draw.qq(sort(out$sumSTAAR.ACAT_O), 'sumSTAAR, without annotation')
draw.qq(sort(out$PCAanno), 'PCA')
draw.qq(sort(out$BTanno), 'BT')
draw.qq(sort(out$SKATanno), 'SKAT')
draw.qq(sort(out$ACATanno), 'ACAT-V')
draw.qq(sort(out$anno11), 'Weighting function parameters (1, 1)')
draw.qq(sort(out$anno125), 'Weighting function parameters (1, 25)')

dev.off()
