getwd()
#load packages
library(tidyverse)
library(multtest)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(ggplot2)
library(plotly)
library(tibble)
library(DESeq2)
#library(biomaRt)
#library(pheatmap)
#library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
#Creating variable Sample_type to get coloums of interest.

# sample_table = read_csv("SraRunTable.txt") %>% 
#   # selecting interested coloumns from the data table
#   ###do no load Biomart otherwise encounter error
#   select(`Sample Name`, source_name,
#          Treatment, Cell_Line, 
#          Cell_type, 
#          Time_point) %>% 
#   #selecting one among the four replicates (DAta iD's)
#   unique()
# another way of doing it, just re-writting the variable again !!
sample_table = read_csv("SraRunTable.txt")
sample_table = select(sample_table, `Sample Name`, source_name,
                      Treatment, Cell_Line, 
                      Cell_type, 
                      Time_point )
sample_table = unique(sample_table)
#install tximport and load library !!!!!
#loading the data file quant.sf files into all for all samples
#pulling the data from directories 
# MAKE SURE YOU HAVE ALL FILES avilable, if not (eeror= all(file.exists(files)) is not TRUE)
######filtering the data we need ########
sample_table_a549=sample_table[sample_table$Cell_Line=='A549',]
samplefiles_a549 = paste0(
  pull(sample_table_a549, `Sample Name`),
  '/quant.sf')
names(samplefiles_a549)=pull(sample_table_a549, `Sample Name`)
#Naming the SF files, thats where the data will be stored otherwise instead of names will be no.'s
names(samplefiles_a549) = pull(sample_table_a549, `Sample Name`)
#extracting the transcript and gene ids from terminal>>>>>> 
#>>>>>>out of reference and stick then in 1 file.
genemap = read_csv("gene_map.csv", col_names = c('enstid','ensgid')) #show_col_types = FALSE)
#Running tximport for abundance, amke sure all files have same names.
countdata_a549 = tximport(files = samplefiles_a549, type = "salmon",
                     tx2gene = genemap,
                     ignoreTxVersion = TRUE)
sample_table_a549$conditions=factor(rep(c('mock', 'infected'), 
                                    each=3), 
                                levels=c('mock', 'infected'))
dds_a549=DESeqDataSetFromTximport(txi = countdata_a549,
                                  colData = sample_table_a549,
                                  design =~ conditions )
#DEseq fuction does all three functions in one go 
#1. sizefactors, dispersion, binomwaldtest
#lets make life easy with one function only

dds_a549=DESeq(dds_a549)
#PCA Plot
vst_a549=varianceStabilizingTransformation(dds_a549)
plotPCA(vst_a549, intgroup='conditions')

######making own PCA PLOT (Prcomp)##### instead from DEseq######

vst_matrix_a549=assay(vst_a549)
pca=prcomp(t(vst_matrix_a549))
df= as.data.frame(pca$x)
df$conditions=sample_table_a549$conditions
#####Percent variance(pvr)##########
pvr=round(pca$sdev^2/sum(pca$sdev^2)*100, 2)

ggplot(df,aes(x=PC1, y=PC2)) +
  geom_point(aes(color=conditions)) +
  xlab(label = paste0('PC1(',pvr[1],'%)'))+
  ylab(label = paste0('PC2(',pvr[2],'%)'))

a549_results=results(dds_a549, contrast = c('conditions', 
                                            'infected', 
                                            'mock'))
summary(a549_results)
View(summary(a549_results))

a549_df=data.frame(a549_results)
#a549_filter1= rownames_to_column(a549_filter1, var = 'ensmblID')
#a549_filter1= a549_results[complete.cases(a549_df),]
a549_filter1=a549_df[complete.cases(a549_df), ]
a549_filter2=a549_filter1[a549_filter1$padj<0.05, ]
a549_filter3=a549_filter2[abs(a549_filter2$log2FoldChange) > 1,]

#plotting results
plotMA(a549_results)
#volcano plot
a549_filter1$test=a549_filter1$padj < 0.05 & abs(a549_filter1$log2FoldChange) > 1
a549_filter1= rownames_to_column(a549_filter1, var = 'ensmblID')
#for colored just create the coloum which will be basis of color
View(a549_filter1)
vol_plot_a549=ggplot(a549_filter1,aes(x=log2FoldChange, y= -log10(padj), name=ensmblID))+
  geom_point(aes(color=test), size=1, alpha=0.3)+
  scale_colour_manual(values = c('darkblue','red'))+
  geom_vline(xintercept = 1, color='green',linetype=2)+
  geom_vline(xintercept = -1, color='green', linetype=2)+
  geom_hline(yintercept = -log10(0.05), color='darkgreen', linetype=3)+
  xlim(-3,3)+
  ylim(0,10)+
  theme_classic() +
  theme(legend.position = 'none')
vol_plot
ggplotly(vol_plot)

#DEseq2 shortcut
#dds_nhbe=DESeq(dds_nhbe)
result_a549 = results(dds_a549)
summary(result_a549)

# DATAFrame
res_a549_df=as.data.frame(result_a549)
View(res_a549_df)


plotCounts(dds_a549, gene = 'ENSG00000265794', intgroup = 'conditions')
sum(complete.cases(res_a549_df))
#filter data to remove na's 
filtered_1a549= res_a549_df[complete.cases(res_a549_df), ]
View(filtered_1a549)

filtered_1a549$padj < 0.05
filtered_2a549=filtered_1a549[filtered_1a549$padj < 0.05,]
View(filtered_2a549)
abs(filtered_2a549$log2FoldChange) > 1
filtered_3a549 = filtered_2a549[abs(filtered_2a549$log2FoldChange) > 1, ]
View(filtered_3a549)
plotMA(result_a549)
filtered_1a549$test=filtered_1a549$padj < 0.05 & abs(filtered_1a549$log2FoldChange) > 1
#volcano plot, -log10 just to not create - pvalue coz log10 0f 0-1 = negative.
filtered_1a549$test=filtered_1a549$padj < 0.05 & abs(filtered_1a549$log2FoldChange) > 1
filtered_1a549= rownames_to_column(filtered_1a549, var = 'ensmblID')
#for colored just create the coloum which will be basis of color
View(filtered_1a549)
# vol_plot=ggplot(filtered_nhbe1,aes(x=log2FoldChange, y= -log10(padj), name=ensmblID))+
#   geom_point(aes(color=test), size=1, alpha=0.3)+
#   scale_colour_manual(values = c('darkblue','red'))+
#   geom_vline(xintercept = 1, color='green',linetype=2)+
#   geom_vline(xintercept = -1, color='green', linetype=2)+
#   geom_hline(yintercept = -log10(0.05), color='darkgreen', linetype=3)+
#   xlim(-3,3)+
#   ylim(0,10)+
#   theme_classic() +
#   theme(legend.position = 'none')
# vol_plot
# #### changing the plot into interactive plot!!!
# #install('plotly')
# ggplotly(vol_plot)>>> as we have got gene id's we can do it below

##Annotation of rnaSEQ DATA DEG's
#install('biomaRt')
library(biomaRt)
listMarts()
#useEnsembl((biomart = 'ensembl'), version = 107)
ensmbl107 = useEnsembl(biomart = "ensembl", version=106)
listEnsembl(ensmbl107)
listDatasets(ensmbl107)
ensmbl107=useDataset(dataset = 'hsapiens_gene_ensembl', mart = ensmbl107)
#view attributes to know what columns we need from the emsembl website and choosing filter those
view(listAttributes(ensmbl107))
view(listFilters(ensmbl107))
#lets check for six genes. does it work, yes, luckily !!!
# getBM(attributes = c('ensembl_gene_id',
#                      'ensembl_gene_id_version',
#                      'ensembl_transcript_id',
#                      'ensembl_transcript_id_version',
#                      'external_gene_name'),
#       filters= c('ensembl_gene_id'), 
#       values = filtered_nhbe1$ensmblID[1:6], mart = ensmbl107)
# now lets run it on our whole data sets redifined our attributes that actually we need

annotation1_a549=getBM(attributes = c('ensembl_gene_id',
                                      'chromosome_name',
                                      'start_position',
                                      'end_position',
                                      'strand',
                                      'gene_biotype',
                                      'external_gene_name',
                                      'description'),
                       filters= c('ensembl_gene_id'), 
                       values = a549_filter1$ensmblID, mart = ensmbl107)
View(annotation1_a549)
View(a549_filter1)

#now we have the annotaion table we need to merge it to our expression data
#but how do we know whether the annotation has been assined to the correct row
#thats where we use dplyr using join method
annotated1_a549= left_join(a549_filter1, annotation1_a549,
                           by= c('ensmblID'='ensembl_gene_id'))
View(annotated1_a549)
vol_plot_a549=ggplot(annotated1_a549,aes(x=log2FoldChange, y= -log10(padj), 
                                    name=external_gene_name))+
  geom_point(aes(color=test), size=1, alpha=0.3)+
  scale_colour_manual(values = c('darkblue','red'))+
  geom_vline(xintercept = 1, color='green',linetype=2)+
  geom_vline(xintercept = -1, color='green', linetype=2)+
  geom_hline(yintercept = -log10(0.05), color='darkgreen', linetype=3)+
  xlim(-3,3)+
  ylim(0,10)+
  theme_classic() +
  theme(legend.position = 'none')
vol_plot_a549
#### changing the plot into interactive plot!!!
#install('plotly')
ggplotly(vol_plot)

#heatmaps
annot_df1_a549=annotated1_a549[annotated1_a549$padj < 0.05,]
annot_df2_a549=annot_df1_a549[abs(annot_df1_a549$log2FoldChange) >1,]
degs_a549=annot_df2_a549$ensmblID
vst_a549_hm = varianceStabilizingTransformation(dds_a549)
vst_a549_matrix=assay(vst_a549_hm)
View(vst_a549_matrix)
data_hm_a549=vst_a549_matrix[degs_a549,]
rownames(data_hm_a549)=annot_df2_a549$external_gene_name
heatmap(data_hm_a549)
#lets get nice pheatmap instead
library(pheatmap)
library(RColorBrewer)
pheatmap(data_hm_a549)
pheatmap(data_hm_a549, scale = 'row', fontsize = 4)
heatmap_RdPu=colorRampPalette(brewer.pal(n=9, name='Spectral' ))(100)
pheatmap(data_hm_a549, scale = 'row', fontsize_row = 4, color = heatmap_RdPu, cutree_cols = 2, cutree_rows = 2 )

##Gene ontology functional enrichment analysis !!!!
#install('clusterProfiler')
#installing the gene ontology database for human 
#install('org.Hs.eg.db')
# GO erichment

library(clusterProfiler)
library(org.Hs.eg.db)
# i was encountering this error here ######Error in .processResults(postRes, 
#mart = mart, hostURLsep = sep, fullXmlQuery = fullXmlQuery,  : #############
#Query ERROR: caught BioMart::Exception::Database: ###########################
#Error during query execution: Table #########################################
#'ensembl_mart_107.hsapiens_gene_ensembl__ox_entrezgene__dm' doesn't exist####
#'#then please check the version your are using for mart !!!!

entgene_id_a549=getBM(attributes = 'entrezgene_id',
                 filters = 'ensembl_gene_id', values = annot_df2_a549$ensmblID, mart = ensmbl107)
entgene_id_a549=entgene_id_a549$entrezgene_id
#creating the variable for universe argument which is the whole gene expression pool ####
ent_unv_a549=getBM(attributes = 'entrezgene_id',
              filters = 'ensembl_gene_id', values = annotated1_a549$ensmblID, mart = ensmbl107)
ent_unv_a549=ent_unv_a549$entrezgene_id
ent_unv_a549= as.character(ent_unv_a549)

ego_a549= enrichGO(gene = entgene_id_a549,
              OrgDb = org.Hs.eg.db,
              ont = 'BP', ######BP=biological process check
              universe = ent_unv_a549,
              readable = TRUE)
summary(ego_a549)
View(summary(ego_a549))

fold_changes_a549=annot_df2_a549$log2FoldChange
names(fold_changes_a549)= annot_df2_a549$external_gene_name
barplot(ego_a549, showCategory = 15)
dotplot(ego_a549)
dotplot(ego_a549, showCategory = 15)
cnetplot(ego_a549,showCategory = 5, foldChange=fold_changes_a549)
goplot(ego_a549)

####KEGG Enrichment ######
kegg_a549=enrichKEGG(gene = entgene_id_a549,
                universe = ent_unv_a549)
View(summary(kegg_a549))
dotplot(kegg_a549)
barplot(kegg_a549)
dotplot(kegg_a549, showCategory = 15)
##########################writing result table to file #########################
write_tsv(annot_df2_a549, 'a549_DEG.txt')
