library(lme4)
library(tidyverse)
library(dplyr)
library(data.table)
library(stringr)
library(sva)
library(doParallel)

nCore = 4
registerDoParallel(cores = nCore)

#args = commandArgs(trailingOnly = TRUE)


args = "/home/lconnelly/MetabolomicsDOE/outputs/posPheno/full_rplc_doe_recentroid/"
args[2] = '/home/lconnelly/MetabolomicsDOE/outputs/posPheno/full_hilic_doe_recentroid/df_set.csv'
args[3] = '/home/lconnelly/MetabolomicsDOE/outputs/posPheno/full_hilic_doe_recentroid/df_sorg.csv'

phenoSet = data.frame(fread(args[2]))
phenoSorg = data.frame(fread(args[3]))
#summary_table = data.frame(fread(args[4]))

phenoSet[,5:ncol(phenoSet)] = log(phenoSet[,5:ncol(phenoSet)])
phenoSorg[,5:ncol(phenoSorg)] = log(phenoSorg[,5:ncol(phenoSorg)])

phenoSet[phenoSet == '-Inf'] = NA
phenoSorg[phenoSorg == '-Inf'] = NA

#phenoSet = phenoSet[,1:80]
#phenoSorg = phenoSorg[,1:80]

#set
#for each column, group by batch, take bottom 30% of non-zero. Sample these to replace NA's. Imputed values get +-5% of their value added/subtracted to induce variance.
imputed_columns = foreach(column = phenoSet[,5:ncol(phenoSet)], .combine = cbind, .errorhandling = "pass") %dopar% {

	batches = phenoSet$Batch
	for(b in unique(phenoSet$Batch)) {
		col_sub = column[!is.na(column) & (batches == b)]
		need_imputing = length(column[is.na(column) & (batches == b)])

                if(need_imputing > 10) {
                        #column[(batches == b)] = NA
                        next
                }

		col_sub = col_sub[col_sub < quantile(col_sub, .33)]
		imputation_values = sample(col_sub, need_imputing, replace = TRUE)
		noised_values = imputation_values + (imputation_values * runif(need_imputing, min = -.05, max = .05))
		column[is.na(column) & (batches == b)] = noised_values 
	}

	return(column)
}

#formatting
num_errors = 0
which_errors = c()
formatted_columns = data.frame(matrix(nrow = nrow(imputed_columns), ncol = 0))
for(i in 1:ncol(imputed_columns)) {
	
	error = 1
	vec = rep(NA, nrow(imputed_columns))
	tryCatch (
	{
		vec = as.numeric(imputed_columns[,i])
		error = 0
	}, 
	error = function(e) {
		which_errors <<- c(which_errors, i)	
	})
	num_errors = num_errors + error
	formatted_columns = cbind(formatted_columns, vec)
}

colnames(formatted_columns) = colnames(phenoSet)[5:ncol(phenoSet)]
pheno_set_blues_ready = formatted_columns
if(length(which_errors) > 0) {
issues = imputed_columns[,which_errors]
pheno_set_blues_ready = pheno_set_blues_ready[,-which_errors]
}
pheno_set_blues_ready = cbind(phenoSet[,1:4], pheno_set_blues_ready)









#sorg

#for each column, group by batch, take bottom 30% of non-zero. Sample these to replace NA's. Imputed values get +-5% of their value added/subtracted to induce variance.
imputed_columns = foreach(column = phenoSorg[,5:ncol(phenoSorg)], .combine = cbind, .errorhandling = "pass") %dopar% {

        batches = phenoSorg$Batch
        for(b in unique(phenoSorg$Batch)) {
                col_sub = column[!is.na(column) & (batches == b)]
                need_imputing = length(column[is.na(column) & (batches == b)])


		if(need_imputing > 10) {
                	column[(batches == b)] = NA
                        next
                }

                col_sub = col_sub[col_sub < quantile(col_sub, .33)]
                imputation_values = sample(col_sub, need_imputing, replace = TRUE)
                noised_values = imputation_values + (imputation_values * runif(need_imputing, min = -.05, max = .05))
                column[is.na(column) & (batches == b)] = noised_values
        }

        return(column)
}

#formatting
num_errors = 0
which_errors = c()
formatted_columns = data.frame(matrix(nrow = nrow(imputed_columns), ncol = 0))
for(i in 1:ncol(imputed_columns)) {

        error = 1
        vec = rep(NA, nrow(imputed_columns))
        tryCatch (
        {
                vec = as.numeric(imputed_columns[,i])
                error = 0
        },
        error = function(e) {
                which_errors <<- c(which_errors, i)
        })
        num_errors = num_errors + error
        formatted_columns = cbind(formatted_columns, vec)
}

colnames(formatted_columns) = colnames(phenoSorg)[5:ncol(phenoSorg)]
pheno_sorg_blues_ready = formatted_columns
if(length(which_errors) > 0) {

issues = imputed_columns[,which_errors]
pheno_sorg_blues_ready = pheno_sorg_blues_ready[,-which_errors]

}
pheno_sorg_blues_ready = cbind(phenoSorg[,1:4], pheno_sorg_blues_ready)




phenoSet = pheno_set_blues_ready
phenoSorg = pheno_sorg_blues_ready
species = rep('Setaria', nrow(phenoSet))
phenoSet = cbind(species, phenoSet)
species = rep('Sorghum', nrow(phenoSorg))
phenoSorg = cbind(species, phenoSorg)


#set_var_before = phenoSet %>% group_by(Batch) %>% summarise(across(c(6:ncol(phenoSet)-1),  ~ var(.x), .names = "variance-{.col}"))
#sorg_var_before = phenoSorg %>% group_by(Batch) %>% summarise(across(c(6:ncol(phenoSorg)-1),  ~ var(.x), .names = "variance-{.col}"))


#x = rbind(phenoSet, phenoSorg)
out = vector('list', 2)
iterator = 0
for(x in list(phenoSet, phenoSorg)) {
#x = phenoSet
iterator = iterator + 1
p1 = x$File
p2 = x$treatment
p3 = x$Batch
p4 = x$genotype
p5 = x$species

pheno_set = data.frame(Treatment = p2, Batch = p3, Genotype = p4, Species = p5)
rownames(pheno_set) = p1



e = x[,6:ncol(x)]
#e[e == '-Inf'] = 0
e = t(e)
colnames(e) = rownames(pheno_set)
expression_set = e
batch = pheno_set$Batch



#removes genotypes for features that are not present in both treatment conditions
tab = table(pheno_set$Genotype, pheno_set$Treatment)
one = tab[,1]
two = tab[,2]
tab2 = data.frame(drought = one, ww = two)
rownames(tab2) = rownames(tab)

to_remove = rownames(tab2)[c(which(tab2$drought == 0), which(tab2$ww == 0), which(tab2$drought == 1), which(tab2$ww == 1))]
to_remove_names = rownames(pheno_set)[pheno_set$Genotype %in% to_remove]
pheno_set = pheno_set[!rownames(pheno_set) %in% to_remove_names,]
expression_set = expression_set[,!colnames(expression_set) %in% to_remove_names]
batch = pheno_set$Batch



#need to make a list of features, and a 2nd list containing the sample names of 0's

#n = rowSums(is.na(expression_set))
#na = table(n)



feature_list = vector('list', nrow(expression_set))
sample_list = vector('list', nrow(expression_set))
batch_list = vector('list', nrow(expression_set))
mod_list = vector('list', nrow(expression_set))


#feature_list = vector('list', 10)
#sample_list = vector('list', 10)
#batch_list = vector('list', 10)
#mod_list = vector('list', 10)



for(i in 1:length(feature_list)) {

print(i / length(feature_list))

mat = matrix(nrow = 0, ncol = ncol(expression_set))

feat = expression_set[i,]

#On a feature by feature basis, check that for each genotype, treatment group, we have at least 1 sample. Or remove the genotype from the feature.
check = cbind(pheno_set, feat)
check2 = check %>% group_by(Genotype, Treatment) %>% summarise(across(feat, ~ sum(!is.na(.x)), .names = 'num_nonzero'))

geno_remove = c()
geno_remove = c(geno_remove, check2$Genotype[check2$num_nonzero == 0 | check2$num_nonzero == 1])
#geno_remove = geno_remove[duplicated(geno_remove)]



mat = data.frame(rbind(mat, feat))

rownames(mat) = rownames(expression_set)[i]

remove = which(is.na(mat))
if(length(remove) == 0) {

remove = 1000000000000

}
remove = c(remove, which(check$Genotype %in% geno_remove))



feature_list[[i]] = mat[,-remove]
sample_list[[i]] = unique(colnames(mat)[remove])
batch_list[[i]] = pheno_set$Batch[-remove]
mod_list[[i]] = pheno_set[-remove,]

}




corrected_expression = foreach(

p = mod_list,
b = batch_list,
f = feature_list,
s = sample_list,
.errorhandling = 'pass'

) %dopar% {



mod = model.matrix(~ Genotype * Treatment, data = p)
f = rbind(f,f)
print(rownames(f)[1])

corrected = ComBat(dat = f, batch = b, mod=mod, par.prior=FALSE, BPPARAM=SerialParam())

corrected = corrected[1,]
return(corrected)

}





cef = foreach(row = corrected_expression, to_add = sample_list, i = seq(1,length(corrected_expression),1), .combine='rbind') %dopar% {

#for(i in 1:length(corrected_expression)) {

	#print(i / length(corrected_expression))
	#row = corrected_expression[[i]]
	if(i %% 10 == 0) {

		print(i / length(corrected_expression))	

	}	
	if(class(row) == "numeric") {
	#to_add = sample_list[[i]]
	if(!is.na(to_add)) {
		
		start = length(row)
		row = c(row, rep(NA, length(to_add)))
		names(row)[(start+1):length(row)] = to_add


	}
	row = data.frame(as.list(row))
	row$corrected = TRUE
	#corrected = c(corrected, i)

	} else {

	row = data.frame(as.list(expression_set[i,]))	
	row$corrected = FALSE
	#not_corrected = c(not_corrected, i)
	}
	
	
	#cef = rbind(cef, row)
	return(row)
}

backup = cef



corr = cef$corrected
cef = cef[-length(colnames(cef))]
#cef = cbind(corr, cef)
combat_edata3 = cef

combat_formatted = data.frame(t(combat_edata3))
colnames(combat_formatted) = rownames(expression_set)
colnames(combat_formatted)[!corr] = paste0(colnames(combat_formatted)[!corr], '_not_corrected')


combat_back = combat_formatted

File = rownames(combat_formatted)
genotype = pheno_set$Genotype
species = pheno_set$Species
treatment = pheno_set$Treatment
batch = pheno_set$Batch


combat_back = cbind(treatment, combat_back)
combat_back = cbind(species, combat_back)
combat_back = cbind(genotype, combat_back)
combat_back = cbind(File, combat_back)
combat_back = cbind(batch, combat_back)









#pheno_out = args[1]
if(iterator == 1) {
write.csv(combat_back, paste0(args[1], 'g_bc_phenotype_object_for_set_blues_impute_by_batch_19_jul_23.csv'), row.names=TRUE)
} else {
write.csv(combat_back, paste0(args[1], 'g_bc_phenotype_object_for_sorg_blues_impute_by_batch_19_jul_23.csv'), row.names=TRUE)
}

combat_formatted = cbind(genotype, combat_formatted)
combat_formatted = cbind(batch, combat_formatted)
combat_formatted = cbind(treatment, combat_formatted)


out[[iterator]] = combat_formatted

}



#out = readRDS('/home/lconnelly/MetabolomicsDOE/outputs/posPheno/full_hilic_doe_recentroid/bc_pheno_objects_both_species.rds')

#phenoSet = combat_formatted[combat_formatted$species == 'Setaria',]
#phenoSorg = combat_formatted[combat_formatted$species == 'Sorghum',]

setInfo = data.frame(fread(paste0(args[1], 'g_bc_phenotype_object_for_set_blues_impute_by_batch_21_june_23.csv')))$treatment
sorgInfo = data.frame(fread(paste0(args[1], 'g_bc_phenotype_object_for_sorg_blues_impute_by_batch_21_june_23.csv')))$treatment

print('made it')

phenoSetAfter = out[[1]]
phenoSorgAfter = out[[2]]

keep_set = colnames(phenoSetAfter)[which(!str_detect(colnames(phenoSetAfter), 'not_corrected'))]
keep_sorg = colnames(phenoSorgAfter)[which(!str_detect(colnames(phenoSorgAfter), 'not_corrected'))]

phenoSetAfter = phenoSetAfter[,keep_set]
phenoSorgAfter = phenoSorgAfter[,keep_sorg]

#iterate thru each column, imputing 0's. Sample bottom 33% of data. Add small amount of noise (+-5% value).


#need to use linear model to extract variance explained by G term.
create_blues = function(pheno_object, info_object, species) {
phenoAfter = pheno_object
Info = info_object
heritability = c()
which_errors = c()
blues = data.frame(matrix(nrow = 0, ncol = length(unique(phenoAfter$genotype))))
colnames(blues) = unique(phenoAfter$genotype)[order(unique(phenoAfter$genotype))]
rows = colnames(phenoAfter)[4:length(colnames(phenoAfter))]
for(i in 4:ncol(phenoAfter)) {
#for(i in 4:12) {
  print(i / ncol(phenoAfter))

  #pheno = data.frame(Compound = as.numeric(phenoAfter[,i]), genotype = phenoAfter$genotype, treatment = Info, batch = phenoAfter$batch)  
  pheno = data.frame(Compound = as.numeric(phenoAfter[,i]), genotype = phenoAfter$genotype, treatment = phenoAfter$treatment, batch = phenoAfter$batch)

  pheno = pheno[pheno$Compound != 0,]
  pheno = pheno[!is.na(pheno$Compound),]

  genotypeCoef = 0

   tryCatch (
        {
		f1.1 = lm(I(Compound) ~ 0 + genotype * treatment, data = pheno)
		#f1.2 = lm(I(Compound) ~ genotype * treatment, data = pheno)

		f1 = lm(I(Compound - f1.1$coefficients[1]) ~ 0 + genotype * treatment, data = pheno)


		#png(paste0('/home/lconnelly/lm/', colnames(phenoAfter)[i], '_', species, '.png'), width = 1920, height = 1080)
		#print(hist(pheno$Compound, breaks = 30))
		#dev.off()
  		genotypeCoef = coef(f1)[1:length(unique(pheno$genotype))]
        },
        error = function(e) {
                which_errors <<- c(which_errors, i)
        })



  if(length(genotypeCoef) == 1) {
	
	blues = rbind(blues, c(rep(NA, 224)))
	next

  }


  names(genotypeCoef) = str_replace(names(genotypeCoef), "genotype", "")
  gorder = names(genotypeCoef)
  
  
  anovaTest = anova(f1)
  heritability =  c(heritability,sum(anovaTest$`Sum Sq`[1:(length(anovaTest$`Sum Sq`) - 1)]) / sum(anovaTest$`Sum Sq`))
  genotypeCoef = data.frame(t(genotypeCoef))

  blues = bind_rows(blues,genotypeCoef)
  
}


rownames(blues) = rows
colnames(blues) = gorder
blues2 = data.frame(t(blues))
blues2 = cbind(rownames(blues2), blues2)
colnames(blues2)[1] = 'genotype'
#write.csv(blues2, '/home/lconnelly/GWASparallel/most_recent/new_sorg_blues_impute_by_batch_may25.csv', row.names=FALSE)

return(list(blues2, heritability))

}

set_var_after = phenoSetAfter %>% group_by(batch) %>% summarise(across(c(4:ncol(phenoSetAfter)-1),  ~ var(.x), .names = "variance-{.col}"))
sorg_var_after = phenoSorgAfter %>% group_by(batch) %>% summarise(across(c(4:ncol(phenoSorgAfter)-1),  ~ var(.x), .names = "variance-{.col}"))




set_out = create_blues(phenoSetAfter, setInfo, 'setaria')
sorg_out = create_blues(phenoSorgAfter, sorgInfo, 'sorghum')

set_blues = set_out[[1]]
sorg_blues = sorg_out[[1]]

set_heritability = set_out[[2]]
sorg_heritability = sorg_out[[2]]

write.csv(set_blues, '/home/lconnelly/GWASparallel/most_recent_rplc_neg/set_blues_july_14_23.csv', row.names=FALSE)
write.csv(sorg_blues, '/home/lconnelly/GWASparallel/most_recent_rplc_neg/sorg_blues_july_14_23.csv', row.names=FALSE)

write.csv(set_heritability, '/home/lconnelly/GWASparallel/most_recent_rplc_neg/set_heritability_june_22_23.csv', row.names=FALSE)
write.csv(sorg_heritability, '/home/lconnelly/GWASparallel/most_recent_rplc_neg/sorg_heritability_june_22_23.csv', row.names=FALSE)


#test_out_set = phenoSet %>% group_by(genotype) %>% summarize_at(vars(colnames(phenoSet)[3:ncol(phenoSet)]), funs(mean(., na.rm=TRUE)) )
#test_out_sorg = phenoSorg %>% group_by(genotype) %>% summarize_at(vars(colnames(phenoSorg)[3:ncol(phenoSorg)]), funs(mean(.,na.rm=TRUE)))

#test_out_set = data.frame(test_out_set)
#test_out_sorg = data.frame(test_out_sorg)

#out_path_set = args[1]
#out_path_sorg = args[1]

#out_path_set = paste0(out_path_set, 'G/')
#out_path_sorg = paste0(out_path_sorg, 'G/')

#if(!dir.exists(out_path_set)) {
#	dir.create(out_path_set,recursive = TRUE)
#}
#if(!dir.exists(out_path_sorg)) {
#	dir.create(out_path_sorg, recursive = TRUE)
#}


#out_path_set = paste0(out_path_set, 'setaria_g_blues.csv')
#out_path_sorg = paste0(out_path_sorg, 'sorghum_g_blues.csv')


#write.csv(test_out_set, paste0(out_path_set), row.names = FALSE)
#write.csv(test_out_sorg, paste0(out_path_sorg), row.names = FALSE)

