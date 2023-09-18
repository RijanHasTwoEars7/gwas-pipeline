###Adapted from https://github.com/Alice-MacQueen/switchgrassGWAS/blob/master/R/pvdiv_gwas.R
#library(bigsnpr)
#library(reshape2)
#library(dplyr)

library(bigsnpr)
library(bigstatsr)
library(reshape2)
library(dplyr)

ncores <- 1
nPCs <- 3



args = commandArgs(trailingOnly = TRUE)
args = unlist(strsplit(args, ','))

tablePath = args[1]
print(tablePath)
print("is tablePath")


index1 = as.numeric(args[2])
index2 = as.numeric(args[3])

#phenotype <- read.table("../data/phenotype/3.daily_setaria_exp.ALL_days.TRANS_STAND.csv", sep = ",", header = T, stringsAsFactors = F) 
#phenotype <- filter(phenotype, !is.na(exp.code))

#phenotype <- read.table("../data/phenotype/ALL_setaria_exp.ALL_days.TRANS_STAND.wet_dry_RATIO.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)


####Do just DAP18 and height in wet, from IB005,IB007,IB008
#phenotype <- phenotype[phenotype$exp.code %in% c("IB005","IB007","IB008"),]
#phenotype <- phenotype[phenotype$dap==18,]
#phenotype <- phenotype[,c("genotype","exp.code","dap",grep("height|area",colnames(phenotype),value=T))]
##########################


#phenotype$exp.day <- paste(phenotype$exp.code,phenotype$dap,sep=".")
#phenotype$exp.code <- NULL
#phenotype$dap <- NULL
#phenotype$exp.dap <- NULL


# single.measure.traits <- c("area.slope",
#                            "CR_num",
#                            "CR.Angle",
#                            "d13C",
#                            "CN_ratio",
#                            "gN.m2",
#                            "gC.m2",
#                            "leaf.area",
#                            "leaf.weight",
#                            "spec.leaf.area",
#                            "stom.density")
# allSingleTraits <- unique(grep(paste(single.measure.traits,collapse="|"), colnames(phenotype), value=TRUE))
#meltPhenotype <- reshape2::melt(phenotype,id.vars=c("genotype","exp.day"))

# add in single value traits
#pheno.single <- read.csv("../data/phenotype/3.singlevalue_setaria_exp.ALL_days.TRANS_STAND.csv", stringsAsFactors = F) %>% 
#  filter(!is.na(exp.code))
#pheno.single$dap <- "single"
#pheno.single$exp.day <- paste(pheno.single$exp.code,pheno.single$dap,sep=".")
#pheno.single$exp.code <- NULL
#pheno.single$dap <- NULL
#meltPhenotype_single <- reshape2::melt(pheno.single,id.vars=c("genotype","exp.day"))

# combine
#names(meltPhenotype_single); names(meltPhenotype)
#meltBoth <- rbind(meltPhenotype, meltPhenotype_single)

# spread to wide format 
#phenotype <- reshape2::dcast(meltBoth,...~variable+exp.day)
#phenotype <- dplyr::select_if(phenotype, function(x) !all(is.na(x)))
#rm(meltPhenotype, meltPhenotype_single, meltBoth)




#add Louis-
print(tablePath)
print(" is tablePath")
phenotype <- read.csv(tablePath, stringsAsFactors = F)
phenotype <- phenotype[ ,c(1,seq(index1,index2,1))]

#print(phenotype)
#print(" is the phenotype")

#phenotype <- reshape2::dcast(meltBoth,...~variable+exp.day)
#phenotype <- dplyr::select_if(phenotype, function(x) !all(is.na(x)))
#rm(meltPhenotype, meltPhenotype_single, meltBoth)


########Do some quick cleanup############
colnames(phenotype)[1] <- "Genotype"
traits <- colnames(phenotype)[2:ncol(phenotype)]


########Do some quick cleanup############
#colnames(phenotype)[1] <- "Genotype"
#traits <- colnames(phenotype)[2:ncol(phenotype)]

####Subset to just traits to run right nwo### - 4 CN traits and dap18 biomass and water use...
#traits <- grep("d13C|CN_ratio|gN.m2|gC.m2|WUE.+\\.18|.+\\.area.+\\.18",traits,value=TRUE)
#traits <- grep("transpire.+\\.18",traits,value=TRUE)
#traits <- grep("over",traits,value=TRUE,invert=TRUE)

# Choosing to cut off all experiments between dap 10 and dap 22.
# The later days get cutoff due to quality of photos at end of experiment relating 
# to plants growing out of the frame of the camera. 
# early days are truncated due to experiments first day differing. Want to 
# run only days with most number of pictures overlapping experiments
#traits <- grep("\\.1[0-9]|\\.2[0-2]", traits, value = T)
#traits <- grep("\\.1[0-9]|\\.2[0-2]|single", traits, value = T)

####These can be used to read a table of specific phenotypes to run#####
# missingPhenos <- read.table("../10b.skippedPhenos.spec.leaf.area.csv",sep=",",header=FALSE,stringsAsFactors = FALSE,comment.char = "")
# phenotype <- phenotype[,c("Genotype",missingPhenos$V1)]
# traits <- missingPhenos$V1
phenotype[mapply(is.infinite, phenotype)] <- NA



#snpFilePath <- "../data/genotype/7B.SUBSETIB008_MAF.1.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.recode_QC.rds"

#MAF .1:
#snpFilePath <- "/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/8.1.IB008_maf.1.hetsIMP.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.recode.rds"

#MAF .05:
#snpFilePath = "/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/lower_maf/8.1.IB008_maf.05.hetsIMP.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.recode.rds"


#Updated for paper:
#snpFilePath = '/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/shattering/8.1.shattering.data_maf.1.hetsIMP.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.recode.rds'

#New SNP File 4/19/2022
#snpFilePath = '/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/shattering/N8.1.shattering.data_maf.1.FilteredGenotypeFile.hetFilter0.25.recode.rds'


#snpFilePath = '/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/new_geno/N.8.1.IB008_maf.1.FilteredGenotypeFile.hetFilter0.25.recode.rds'

snpFilePath = '/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/genotype_filtering/option_III_1.5mil/final_results/8.2.IB008_maf.1.maxmaf.9.hetsIMP.filteredSNPs.2kbDistThresh.0.5neighborLD.hetFilter0.25.recode.rds'

snp <- snp_attach(snpFilePath)

G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos

#Change chromsomes to numeric
CHRN <- as.numeric(gsub("Chr_","",CHR))

##########Create a table with a row for each ID in the genotype file and the name it is translated to from the phenotype file####
#This can then be used to send a set of indices to the bigsnpr GWAS function#########
#######Open csv containing info about lines in Genotype file#####
genoInfo <- read.table("/home/cluebbert/from_github/SetariaGWASmash/data/genotype/Setaria_597_diversity_samples.csv",sep=",",header=TRUE,stringsAsFactors = FALSE,comment.char = "")

#genoInfo$Genotype <- gsub("_setaria_12","",genoInfo$New_name)

genoInfo$Genotype <- genoInfo$New_name
names(genoInfo)[1] = 'LIB'

print('crash before sapply')

genoInfo$repLib <- sapply(genoInfo$LIB,function(x){paste(rep(x,4),collapse = "_")})
genoInfo <- genoInfo[which(genoInfo$repLib %in% intersect(genoInfo$repLib,snp$fam$sample.ID)),]


#genoInfo = sapply(intersect(genoInfo$repLib,snp$fam$sample.ID), function(x){which(x==genoInfo$repLib)})

keepLines <- data.frame(from=sapply(genoInfo$LIB[genoInfo$Genotype %in% c("A10.1",intersect(genoInfo$Genotype,phenotype$Genotype))],function(x){paste(rep(x,4),collapse = "_")}),
                        to=genoInfo$Genotype[genoInfo$Genotype %in% c("A10.1",intersect(genoInfo$Genotype,phenotype$Genotype))],stringsAsFactors = FALSE)

# creat folder for "parallel" files
#dir.create("../temp_progress", showWarnings = FALSE)

#covar <- readRDS("../results/9A.10PCs.SVD_object_subsetIB008.rds")
#covar <- readRDS("/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/results/9.IB008_10PCs.SVD_object.rds")

#covar = readRDS('/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/shattering/9.shattering.data_10PCs.SVD_object.rds')

#covar = readRDS('/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/shattering/N9.shattering.data_10PCs.SVD_object.rds')

#covar = readRDS('/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/data/genotype/new_geno/N9.IB008_10PCs.SVD_object.rds')

covar = readRDS('/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash/genotype_filtering/option_III_1.5mil/9.III.IB008_10PCs.SVD_object.rds')

phenotype = phenotype[phenotype$Genotype %in% genoInfo$New_name,]

libsInOrder = c()
#add corresponding LIB to each genotype
for(genotype in phenotype$Genotype) {
	lib = genoInfo$LIB[genoInfo$New_name == genotype]
	libsInOrder = c(libsInOrder, lib)
}


print(length(phenotype$Genotype))
print(length(libsInOrder))

phenotype$LIB = libsInOrder

phenotype$Genotype <- sapply(phenotype$LIB,function(x){paste(rep(x,4),collapse = "_")})

#/shares/ibaxter_share/private/cluebbert/from_github/SetariaGWASmash


output = data.frame(matrix(ncol = 4, nrow = 0))
colnames(output) = c('metabolite', 'pval', 'chromo', 'pos')


#print(traits)
#print(" is the traits")
#traits = traits[traits %in% c('X252.0446_1642', 'X351.10453_238', 'X407.26827_233', 'X181.04926_614')]


for(i in sample(traits)){


 
  # create dummy files which indicate that a trait has already been analyzed, a way to do a type of parallel processing 
#  if(file.exists(paste0("/home/lconnelly/GWASparallel/temp_progress/workingon", i))){
#    next;
#  }else{
#    write.table(i, paste0("/home/lconnelly/GWASparallel/temp_progress/workingon", i)) 
#  }
  
  message("Peforming GWAS on ",i,"\n")
  ####Make a numeric vector of the phenotype and the individuals####
  thisPheno = as.numeric(as.matrix(phenotype[,i])) 
  print(thisPheno)
  print(" is thisPheno")
  thisPheno <- thisPheno[which(!(is.na(thisPheno)))]
  if(length(thisPheno) < 75){
    message("Trait ",i,"only had ",length(thisPheno)," observations. Skipping.\n")
    next;
  }

  phenotypedLines <- phenotype$Genotype[which(!(is.na(phenotype[,i])))]
#  phenotypedLines <- keepLines$from[which(keepLines$to %in% phenotype$Genotype[which(!(is.na(phenotype[,i])))])]
#  if(length(phenotypedLines) != length(thisPheno)){
#    thisPheno <- phenotype[which(phenotype$Genotype %in% keepLines$to),i]
#    thisPheno = as.numeric(as.matrix(thisPheno)) 
#    thisPheno <- thisPheno[which(!(is.na(thisPheno)))]
#  }
  genoLineIndx <- sapply(phenotypedLines, function(x){which(x == snp$fam$sample.ID)})
  print(genoLineIndx)
  print(identical(snp$fam$sample.ID[genoLineIndx], phenotypedLines)) 
  print('if ^^^^ is false, something went wrong')

  # if statement to allow for running without PC's if desired
  if(!is.na(covar[1])){
     ind_u <- matrix(covar$u[genoLineIndx,1:nPCs], ncol = nPCs)
     gwaspc <- big_univLinReg(G, y.train = thisPheno, covar.train = ind_u,
                              ind.train = genoLineIndx, ncores = ncores)
  }
#else {
#    gwaspc <- bigstatsr::big_univLinReg(G, y.train = thisPheno, ind.train = genoLineIndx,
#                             ncores = ncores)
#  }


  #png(paste0('/home/lconnelly/GWASparallel/most_recent/result_set/qq_plots/', i, '.png'), width = 720, height = 480, units = 'px')
  #print(snp_qq(gwaspc))
  #dev.off()



  
  #if(!file.exists(paste0("/home/lconnelly/GWAS/gwas/10.GWAS_object_", i, ".rds"))) {  
  saveRDS(gwaspc, file = paste0("/home/lconnelly/GWASparallel/most_recent/result_set/gwas_obj/10.GWAS_object_", i, ".rds"))
  #}


if(FALSE) {
  pvals = predict(gwaspc)
  chromo = CHRN[which(pvals == min(pvals, na.rm = TRUE))][1]
  posi = POS[which(pvals == min(pvals, na.rm = TRUE))][1]
  pvals = min(pvals, na.rm = TRUE)
  pvals = -pvals
  i = substr(i, 2,nchar(i))
  print(paste(pvals, 'is max log10p pval for metabolite', i))

  print(i)
  print(" is i")
  print(pvals)
  print("is pvals")

  de<-data.frame(i,pvals, chromo, posi)
  names(de)<-names(output)

  output <- rbind(output, de)
  #output = rbind(output, as.data.frame(c(i, pvals)))
}  

}
colnames(output) = c('metabolite', '-log10p', 'chromo', 'posi')
#print(head(output))
#print(paste0('/home/ahubbard/GWASparallel/GWASresults/GWASResultIndex',index1,'.csv'))
#write.csv(output, file = paste0('/home/lconnelly/GWASparallel/most_recent/result_set/GWASResultIndex',index1,'.csv'), row.names=F)













#for(i in sample(traits)){
  # create dummy files which indicate that a trait has already been analyzed, a way to do a type of parallel processing 
#  if(file.exists(paste0("../temp_progress/workingon", i))){
#    next;
#  }else{
#    write.table(i, paste0("../temp_progress/workingon", i)) 
#  }
#  
#  message("Peforming GWAS on ",i,"\n")
#  ####Make a numeric vector of the phenotype and the individuals####
#  thisPheno = as.numeric(as.matrix(phenotype[,i])) 
#  thisPheno <- thisPheno[which(!(is.na(thisPheno)))]
#  if(length(thisPheno) < 75){
#    message("Trait ",i,"only had ",length(thisPheno)," observations. Skipping.\n")
#    next;
#  }
#  phenotypedLines <- keepLines$from[which(keepLines$to %in% phenotype$Genotype[which(!(is.na(phenotype[,i])))])]
#  if(length(phenotypedLines) != length(thisPheno)){
#    thisPheno <- phenotype[which(phenotype$Genotype %in% keepLines$to),i]
#    thisPheno = as.numeric(as.matrix(thisPheno)) 
#    thisPheno <- thisPheno[which(!(is.na(thisPheno)))]
#  }
#  genoLineIndx <- which(snp$fam$sample.ID %in% phenotypedLines)
#  
#  # if statement to allow for running without PC's if desired
#  if(!is.na(covar[1])){
#     ind_u <- matrix(covar$u[genoLineIndx,1:nPCs], ncol = nPCs)
#     gwaspc <- big_univLinReg(G, y.train = thisPheno, covar.train = ind_u,
#                              ind.train = genoLineIndx, ncores = ncores)
#  } #else {
##    gwaspc <- bigstatsr::big_univLinReg(G, y.train = thisPheno, ind.train = genoLineIndx,
##                             ncores = ncores)
##  }
#  
#  saveRDS(gwaspc, file = paste0("../GWASresults/10.GWAS_object_", i, ".rds"))
#  
#}





