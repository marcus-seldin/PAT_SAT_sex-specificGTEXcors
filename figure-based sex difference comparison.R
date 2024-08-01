
################GTEx Integration

load('G:/My Drive/lab files/sex-difference myokine study/github/raw files/GTEx NA included env.RData')

library(WGCNA)
library(RColorBrewer)
library(ggVennDiagram)
library(pheatmap)
library(MetBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(forcats)
library(enrichR)
working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
test1 = working_dataset[,grepl('ITIH5', colnames(working_dataset))]
colnames(test1)
sex_table = read.delim('G:/My Drive/lab files/sex-difference myokine study/github/raw files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]

working_datasetm = working_dataset[row.names(working_dataset) %in% sex_table$GTEx_ID[sex_table$sexMF=='M'],]
sec_prots = read.delim('G:/My Drive/Datasets/Human/genome files/human secreted proteins.tab')
sec_IDs = as.vector(paste0(sec_prots$Gene.names...primary.., '_Pancreas'))
target = working_dataset[,grepl('Pancreas', colnames(working_dataset) )]
targ_ids = paste0(res1$ID[res1$logFC>0.5 & res1$P.Value<0.01], '_Pancreas')
target = target[,colnames(target) %in% targ_ids]
origin = working_dataset[,colnames(working_dataset) %in% sec_IDs]
######################################################
#run cross-tissue correlations
#initially only run secreted proteins:
tissue.tissue.p = bicorAndPvalue(origin, target, use='pairwise.complete.obs')

tt11 = tissue.tissue.p$p
tt11[is.na(tt11)] = 0.5
tt11[tt11==0] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt11)))

colnames(cc3) = 'Ssec_score'
cc3$gene_symbol = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
summary(cc3$Ssec_score)

scores=cc3
scores$gene_tissue = scores$gene_symbol
scores$gene_symbol = gsub("\\_.*","", scores$gene_tissue)
scores$tissue = gsub(".*_","", scores$gene_tissue)
head(scores)

########################################################
#filter for tissue-specificity
tissue1 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset)) | grepl('Brain - Hippocampus', colnames(working_dataset))  | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]

melted_expr = reshape2::melt(as.matrix(tissue1))
colnames(melted_expr) = c('GTExID', 'gene_tissue', 'value')
melted_expr$gene_symbol = gsub("\\_.*","", melted_expr$gene_tissue)
melted_expr$tissue = gsub(".*_","", melted_expr$gene_tissue)
sum_stats = melted_expr %>% dplyr::group_by(gene_symbol, tissue) %>% dplyr::summarise(mean=mean(value, na.rm=T), sd=sd(value, na.rm=T), n=n())

cross_tissue_stats = sum_stats %>%  dplyr::group_by(gene_symbol) %>% summarise(mean=mean(mean), sdmean = sd(sd))
cross_tissue_stats$mean_plussd = cross_tissue_stats$mean + cross_tissue_stats$sdmean
sum_stats$sd_cut_off = cross_tissue_stats$mean_plussd[match(sum_stats$gene_symbol, cross_tissue_stats$gene_symbol)]
head(sum_stats)
expre_set = sum_stats[sum_stats$mean>sum_stats$sd_cut_off,]
expre_set$gene_tissue = paste0(expre_set$gene_symbol, '_', expre_set$tissue)
expre_set$gene_tissue[1:30]
scores=cc3
scores$gene_tissue = scores$gene_symbol
scores$gene_symbol = gsub("\\_.*","", scores$gene_tissue)
scores$tissue = gsub(".*_","", scores$gene_tissue)
scores = scores[scores$gene_tissue %in% expre_set$gene_tissue,]
head(scores)
scores$gene_symbol[1:20]
tt22 = scores
tt22$category = paste0('M')

males_ssec = tt22
################
#Females
working_datasetm = working_dataset[row.names(working_dataset) %in% sex_table$GTEx_ID[sex_table$sexMF=='F'],]
sec_prots = read.delim('G:/My Drive/Datasets/Human/genome files/human secreted proteins.tab')
sec_IDs = as.vector(paste0(sec_prots$Gene.names...primary.., '_Pancreas'))
target = working_dataset[,grepl('Pancreas', colnames(working_dataset) )]
targ_ids = paste0(res1$ID[res1$logFC<0.5 & res1$P.Value<0.01], '_Pancreas')
target = target[,colnames(target) %in% targ_ids]
origin = working_dataset[,colnames(working_dataset) %in% sec_IDs]
######################################################
#run cross-tissue correlations
#initially only run secreted proteins:
tissue.tissue.p = bicorAndPvalue(origin, target, use='pairwise.complete.obs')

tt11 = tissue.tissue.p$p
tt11[is.na(tt11)] = 0.5
tt11[tt11==0] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt11)))

colnames(cc3) = 'Ssec_score'
cc3$gene_symbol = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
summary(cc3$Ssec_score)

scores=cc3
scores$gene_tissue = scores$gene_symbol
scores$gene_symbol = gsub("\\_.*","", scores$gene_tissue)
scores$tissue = gsub(".*_","", scores$gene_tissue)
head(scores)

########################################################
#filter for tissue-specificity
tissue1 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset)) | grepl('Brain - Hippocampus', colnames(working_dataset))  | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]

melted_expr = reshape2::melt(as.matrix(tissue1))
colnames(melted_expr) = c('GTExID', 'gene_tissue', 'value')
melted_expr$gene_symbol = gsub("\\_.*","", melted_expr$gene_tissue)
melted_expr$tissue = gsub(".*_","", melted_expr$gene_tissue)
sum_stats = melted_expr %>% dplyr::group_by(gene_symbol, tissue) %>% dplyr::summarise(mean=mean(value, na.rm=T), sd=sd(value, na.rm=T), n=n())

cross_tissue_stats = sum_stats %>%  dplyr::group_by(gene_symbol) %>% summarise(mean=mean(mean), sdmean = sd(sd))
cross_tissue_stats$mean_plussd = cross_tissue_stats$mean + cross_tissue_stats$sdmean
sum_stats$sd_cut_off = cross_tissue_stats$mean_plussd[match(sum_stats$gene_symbol, cross_tissue_stats$gene_symbol)]
head(sum_stats)
expre_set = sum_stats[sum_stats$mean>sum_stats$sd_cut_off,]
expre_set$gene_tissue = paste0(expre_set$gene_symbol, '_', expre_set$tissue)
expre_set$gene_tissue[1:30]
scores=cc3
scores$gene_tissue = scores$gene_symbol
scores$gene_symbol = gsub("\\_.*","", scores$gene_tissue)
scores$tissue = gsub(".*_","", scores$gene_tissue)
scores = scores[scores$gene_tissue %in% expre_set$gene_tissue,]
head(scores)
scores$gene_symbol[1:20]
tt22 = scores
tt22$category = paste0('F')
females_ssec = tt22

gene_list1 = females_ssec$gene_symbol[1:10]
tt1 = as.data.frame(rbind(females_ssec[females_ssec$gene_symbol %in% gene_list1,], males_ssec[males_ssec$gene_symbol %in%gene_list1,]))
table(tt1$category)
pdf(file = paste0('top-ranked female proteins.pdf'))
ggplot(tt1, aes(x=fct_reorder2(gene_symbol, Ssec_score, Ssec_score, .desc = T), y=Ssec_score, fill=category)) + geom_col(position = 'dodge') +  ggtitle(paste0('top-ranked female secreted proteins')) + theme_minimal()  + xlab('gene') + ylab('Ssec score') + geom_hline(yintercept = mean(scores$Ssec_score), linetype='dashed', color='gray40', cex=3) + theme(axis.text.x = element_text(angle=90, size=8,vjust =0.5, hjust = 0.5), plot.title = element_text(hjust=0.5))  + xlab('')

dev.off()


gene_list1 = males_ssec$gene_symbol[1:10]
tt1 = as.data.frame(rbind(females_ssec[females_ssec$gene_symbol %in% gene_list1,], males_ssec[males_ssec$gene_symbol %in%gene_list1,]))
table(tt1$category)
pdf(file = paste0('top-ranked Male proteins.pdf'))
ggplot(tt1, aes(x=fct_reorder2(gene_symbol, Ssec_score, Ssec_score, .desc = T), y=Ssec_score, fill=category)) + geom_col(position = 'dodge') +  ggtitle(paste0('top-ranked Mmale secreted proteins')) + theme_minimal()  + xlab('gene') + ylab('Ssec score') + geom_hline(yintercept = mean(scores$Ssec_score), linetype='dashed', color='gray40', cex=3) + theme(axis.text.x = element_text(angle=90, size=8,vjust =0.5, hjust = 0.5), plot.title = element_text(hjust=0.5))  + xlab('')

dev.off()
intersect(males_ssec$gene_symbol[1:10], females_ssec$gene_symbol[1:10])
