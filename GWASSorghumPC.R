
###Adapted from https://github.com/Alice-MacQueen/switchgrassGWAS/blob/master/R/pvdiv_gwas.R
library(bigsnpr)
library(reshape2)
library(dplyr)
ncores <- 1
nPCs <- 4


args = commandArgs(trailingOnly = TRUE)
args = unlist(strsplit(args, ','))


tablePath = args[1]

index1 = as.numeric(args[2])
index2 = as.numeric(args[3])


#phenotype <- readRDS("/shares/ibaxter_share/private/cluebbert/Sorghum_Bellwethers/results/P.PhenotypeSTAND_subset50_forbenchmarking.rds")
phenotype <- read.csv(tablePath, stringsAsFactors = F)
phenotype <- phenotype[ , c(1,seq(index1,index2,1))]




########Do some quick cleanup############
colnames(phenotype)[1] <- "Genotype"
traits <- colnames(phenotype)[2:ncol(phenotype)]



#remove specific genotypes/traits 


phenotype[mapply(is.infinite, phenotype)] <- NA


snpFilePath <- "/shares/ibaxter_share/private/cluebbert/Sorghum_Bellwethers/data/genotype/3.1.IB009_maf.1.hetsIMP.het.25.Sorghum_LT_lines_chrom_only.recode.1maf.9missing.recode.rds"

#snpFilePath = '/home/lconnelly/sorgSNPfile/2.1.IB008het.25_Sorghum_LT_lines_chrom_only.recode.1maf.9missing.rds'


#newest 
#snpFilePath = '/shares/ibaxter_share/private/cluebbert/Sorghum_Bellwethers/data/genotype/3.1.IB009_maf.1.hetsIMP.het.25.Sorghum_LT_lines_chrom_only.recode.1maf.9missing.recode.rds'
#snpFilePath = B008_maf.1.maxmaf.9.hetsIMP.filteredSNPs.2kbDistThresh.0.5neighborLD.hetFilter0.25.recode.rds'


snp <- snp_attach(snpFilePath)

G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos

#Change chromsomes to numeric
CHRN <- as.numeric(gsub("Chr_","",CHR))

#covar <- readRDS("/shares/ibaxter_share/private/cluebbert/Sorghum_Bellwethers/data/genotype/4.IB009_hetIMP.10PCs.SVD_object.rds")
#covar = readRDS('/home/lconnelly/sorgSNPfile/covarFile/pcs_sorg_all_genotypes.rds')

#newest
covar = readRDS('/shares/ibaxter_share/private/cluebbert/Sorghum_Bellwethers/results/4.IB009_10PCs.SVD_object.rds')

# iteration counter and first row of storage matrix
#c = 2
#time_out <- system.time({ print("Starting Loop") })
#meta_out <- data.frame(num_genos = 0, num_snps = 0)
output = data.frame(matrix(ncol = 5, nrow = 0))
colnames(output) = c('metabolite', 'pval', 'chr', 'pos', 'pcs')


correlation_info = data.frame()


for(i in traits){
  # # create dummy files which indicate that a trait has already been analyzed, a way to do a type of parallel processing 
  # if(file.exists(paste0("../temp_progress/workingon", i))){
  #   next;
  # }else{
  #   write.table(i, paste0("../temp_progress/workingon", i)) 
  # }
  
  #this.time <- system.time({
  message("Peforming GWAS on ",i,"\n")
  
  
  
  # Make a numeric vector of the phenotype and the individuals, those that have genotype data are a subset of those that were phenotyped
  phenotyped_lines <- phenotype$Genotype[which(!is.na(phenotype[,i]))] # get genotypes with data for current trait
  #print(phenotyped_lines)
  pheno_lines.in_geno <- phenotyped_lines[which(phenotyped_lines %in% snp$fam$sample.ID)] # get the subset of these genos that have genotype data
  #print(pheno_lines.in_geno)

  thisPheno.vec <- as.numeric(as.matrix(phenotype[,i]))[phenotype$Genotype %in% pheno_lines.in_geno] # get vector of phenotypes for this subset of genos
   if(length(thisPheno.vec) < 75){
    message("Trait ",i,"only had ",length(thisPheno.vec)," observations. Skipping.\n")
    next;
  }

  print(thisPheno.vec)  
  test = shapiro.test(thisPheno.vec)
  if(test[2] > .05) {
	normal = TRUE
  } else {
	normal = FALSE
  }
  # get index in snp file of included genotypes for this phenotype
  genoLineIndx <- sapply(pheno_lines.in_geno, function(x){which(x == snp$fam$sample.ID)})
  all_pc = matrix(covar$u[genoLineIndx,1:10], ncol = 10)
  print(identical(snp$fam$sample.ID[genoLineIndx], pheno_lines.in_geno))
  # thisPheno = as.numeric(as.matrix(phenotype[,i])) 
  # phenoIndx <- which(!(is.na(thisPheno)))
  # thisPheno.vec <- thisPheno[phenoIndx]
  #print(ggdensity(thisPheno.vec))
  # check for appropriate amount of data
  #Sys.sleep(5)
  pcs_to_include = c(1,2,3,4)
  for(j in 5:8) {

  if(cor.test(thisPheno.vec, all_pc[,j])[3] < .05) {

	pcs_to_include = c(pcs_to_include, j)
	
  }
  }
  pc_correlations = cor(thisPheno.vec, all_pc)
  pc_correlations = c(i, pc_correlations, normal)
  correlation_info = rbind(correlation_info, pc_correlations)  
  colnames(correlation_info) = c('trait', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8', 'pc9', 'pc10', 'shapiro_test_normality')  
  # read in pc info for these genos
  ind_u <- matrix(covar$u[genoLineIndx,1:10], ncol = 10)
  ind_u = ind_u[,pcs_to_include]
  # run gwas
  gwaspc <- big_univLinReg(G, y.train = thisPheno.vec, covar.train = ind_u,
                           ind.train = genoLineIndx, ncores = ncores)
 
if(TRUE) { 
  pvals = predict(gwaspc)
  chromo = CHRN[which(pvals == min(pvals, na.rm = TRUE))][1]
  posi = POS[which(pvals == min(pvals, na.rm=TRUE))][1]
  pvals = min(pvals, na.rm = TRUE)
  pvals = -pvals
  i = substr(i, 2,nchar(i))
  print(paste(pvals, 'is max log10p pval for metabolite', i))

  print(i)
  print(" is i")
  print(pvals)
  print("is pvals")

  de<-data.frame(i,pvals,chromo,posi, paste(pcs_to_include, collapse = ', '))
  names(de)<-names(output)

  output <- rbind(output, de)
}

  # save gwas output as RDS
  png(paste0('/home/lconnelly/GWASparallel/most_recent/result_sorg/qq_plots/', i, '.png'), width = 720, height = 480, units = 'px') 
  print(snp_qq(gwaspc))
  dev.off()

#saveRDS(gwaspc, file = paste0("/home/lconnelly/GWASparallel/most_recent/result_sorg/gwas_obj/4.GWAS_object_X", i, ".rds"))
 
  #})
  
  
  #this.meta <- data.frame(num_genos = length(genoLineIndx),
   #                       num_snps = nrow(snp$map))
  #meta_out <- rbind(meta_out, this.meta)

  
  #time_out <- rbind(time_out, this.time)
  #rownames(time_out)[c] <- i
  #c <- c + 1
}



colnames(output) = c('metabolite', '-log10p', 'chromo', 'pos', 'pcs')
print(head(output))
write.csv(output, file = paste0('/home/lconnelly/GWASparallel/most_recent/result_sorg/GWASResultIndex',index1,'.csv'), row.names=F)


#time_out.meta <- cbind(time_out, meta_out)
#time_out.meta <- time_out.meta[-1, ]

#write.csv(time_out.meta, "../results/4.GWAS_benchmark_times_sorghum_subset_of_50_phenos.withmeta.csv", row.names = T)
#write.csv(time_out, "../results/4.GWAS_benchmark_times_sorghum_subset_of_50_phenos.csv", row.names = T)
 



