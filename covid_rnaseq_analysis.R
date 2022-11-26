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
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
#Creating variable Sample_type to get coloums of interest.

# #sample_table = read_csv("SraRunTable.txt") %>% 
#   # selecting interested coloumns from the data table
#   ###do no load Biomart otherwise encounter error
#   #select(sample_table, `Sample Name`, source_name,
#                  Treatment, Cell_Line, 
#                  Cell_type, 
#                  Time_point) %>% 
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
samplefiles = paste0(
  pull(sample_table, `Sample Name`),
  '/quant.sf')
names(samplefiles)=pull(sample_table, `Sample Name`)
#Naming the SF files, thats where the data will be stored otherwise instead of names will be no.'s
names(samplefiles) = pull(sample_table, `Sample Name`)
#extracting the transcript and gene ids from terminal>>>>>> 
#>>>>>>out of reference and stick then in 1 file.
genemap = read_csv("gene_map.csv", col_names = c('enstid','ensgid')) #show_col_types = FALSE)
#Running tximport for abundance, amke sure all files have same names.
countdata_a549 = tximport(files = samplefiles, type = "salmon",
                     tx2gene = genemap,
                     ignoreTxVersion = TRUE)
#load library ‘DESeq2’
library(DESeq2)
#setting the colDATA
sample_table=as.data.frame(sample_table)
#changing the `Sample Name` to sample, Keeping it simple
colnames(sample_table)[1]='sample'
conditions=c('mock_nhbe', 'infected_nhbe', 'mock_a549', 'infected_a549')
conditions=rep(conditions, each=3)
#adding the new coloumn conditions to the sample table 
sample_table$conditions = conditions

dseq_data = DESeqDataSetFromTximport(txi = countdata, 
                                     colData = sample_table,
                                     design =~ conditions )
#Normalization of dseq_data, in_depth library prep. ~sizefactors
dseq_data=estimateSizeFactors(dseq_data)
normalizationFactors(dseq_data)
counts(dseq_data, normalized= TRUE)[1:6,1:3]
boxplot(counts(dseq_data, normalized= TRUE))

#PCA #DAta must be normally distributed, but log transformed rnaseq data is almost normally distributed ->>
# there is log transformation and variance stabilizing transformation(VST-here)
#applying here variance stab transformation
vst=varianceStabilizingTransformation(dseq_data)
boxplot(assay(vst))
#after this normalization, we can do principle component analysis
# intergroup is that which we prepared before and merged into the sample table!
# thats the conditions which actually is mock vs infected expriment
plotPCA(vst, intgroup='conditions')
#after looking at the plot this seems that the cell lines are completely different or
# completely 2 different experiment, thus normalising these two togathr nhbe and a549 cell line is completely waste
#thus we are making two data's out of this oe will be for nhbe and other a549 (mock vs treated)
#splitting the exp into 2 datasets
dseq_data
dds1_nhbe= dseq_data[,1:6]
#dds2=dseq_data[,7:12]
# we have splitted the data [based on coloumns]we have to run normalization again
#to nullify the effect of other data that it was splitted from
dds1_nhbe=estimateSizeFactors(dds1_nhbe)
#whether the normalisation has worked normalizationFactors(dds1 vs dseq_data)[1:6,1:3]
normalizationFactors(dds1_nhbe)
#dds2=estimateSizeFactors(dds2)
#whether the normalisation has worked normalizationFactors(dds1 vs dseq_data)
#normalizationFactors(dds2)

##checking the difference bw splitted and combined data (few rows and coloums)
counts(dseq_data, normalized=TRUE)
counts(dds1_nhbe,normalized=TRUE )
#counts(dseq_data, normalized=TRUE)
#counts(dds2, normalized=TRUE)
#creating vst normalized datasts
vst1_nhbe=varianceStabilizingTransformation(dds1_nhbe)
plotPCA(vst1_nhbe,  intgroup='conditions')

#same for 2nd datasets dds2
#vst2=varianceStabilizingTransformation(dds2)
#plotPCA(vst2,intgroup='conditions')

#clustering genes
# for hierarchical clustering we need to calculate distance matrix
#
dm1_nhbe=assay(vst1_nhbe)
dm1_nhbe=t(dm1_nhbe)
dm1_nhbe=dist(dm1_nhbe)
# or this way simple
dist(t(assay(vst1_nhbe)))
cluster1_nhbe=hclust(dm1_nhbe)
plot(cluster1_nhbe)

#kmeans clustering for 2 groups
kmns1_nhbe=kmeans(t(assay(vst1_nhbe)), centers = 2)
#plot(kmns1_nhbe)
kmns1_nhbe$cluster

#3 steps to deseq analysis
#1.normalization of data from tximport(estimate size factor)
#2. Estimate the dispersion or variance (biological variation, technical variation(shortnoise))
#3. Apply statistics(Wald test)


#counts1_nhbe=counts(dds1_nhbe)
#sample_table1_nhbe=sample_table[1:6, ]

#dds1_nhbe=estimateDispersions(dds1_nhbe)



#Above doesnt work, as we imported data to DEseq2 which was combined nhbe+a549
#split the data out before importing into DEseq2
#splitting the sample files as well 
samplefiles_nhbe=samplefiles[1:6]
countdata_nhbe=tximport(files = samplefiles_nhbe, 
                        type = 'salmon',
                        tx2gene = genemap,
                        ignoreTxVersion = TRUE)
sample_table_nhbe=sample_table[1:6, ]


sample_table_nhbe$conditions=factor(rep(c('mock', 'infected'), 
                                    each=3), 
                                levels=c('mock', 'infected'))
dds_nhbe=DESeqDataSetFromTximport(txi = countdata_nhbe,
                                  colData = sample_table_nhbe,
                                  design =~ conditions )

dds_nhbe=estimateSizeFactors(dds_nhbe)
# normalizationFactors(dseq_nhbe)
# counts(dseq_nhbe, normalized= TRUE)[1:6,1:3]
# boxplot(counts(dseq_nhbe, normalized= TRUE))
#2. Estimate the dispersion or variance (biological variation, technical variation(shortnoise))
dds_nhbe=estimateDispersions(dds_nhbe)

plotDispEsts(dds_nhbe)

#3. Apply statistics(Wald test)
dds_nhbe=nbinomWaldTest(dds_nhbe)


#DEseq2 shortcut
#dds_nhbe=DESeq(dds_nhbe)
result_nhbe = results(dds_nhbe)
summary(result_nhbe)

# DATAFrame
res_nhbe_df=as.data.frame(result_nhbe)
View(res_nhbe_df)


plotCounts(dds_nhbe, gene = 'ENSG00000265794', intgroup = 'conditions')
sum(complete.cases(res_nhbe_df))
#filter data to remove na's 
filtered_nhbe1= res_nhbe_df[complete.cases(res_nhbe_df), ]
View(filtered_nhbe1)

filtered_nhbe1$padj < 0.05
filtered_nhbe2=filtered_nhbe1[filtered_nhbe1$padj < 0.05,]
View(filtered_nhbe2)
abs(filtered_nhbe2$log2FoldChange) > 1
filtered_nhbe3 = filtered_nhbe2[abs(filtered_nhbe2$log2FoldChange) > 1, ]
View(filtered_nhbe3)
plotMA(result_nhbe)
#volcano plot, -log10 just to not create - pvalue coz log10 0f 0-1 = negative.
filtered_nhbe1$test=filtered_nhbe1$padj < 0.05 & abs(filtered_nhbe1$log2FoldChange) > 1
filtered_nhbe1= rownames_to_column(filtered_nhbe1, var = 'ensmblID')
#for colored just create the coloum which will be basis of color
View(filtered_nhbe1)
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

annotation_nbhe1=getBM(attributes = c('ensembl_gene_id',
                     'chromosome_name',
                     'start_position',
                     'end_position',
                     'strand',
                     'gene_biotype',
                     'external_gene_name',
                     'description'),
      filters= c('ensembl_gene_id'), 
      values = filtered_nhbe1$ensmblID, mart = ensmbl107)
View(annotation_nbhe1)
View(filtered_nhbe1)

#now we have the annotaion table we need to merge it to our expression data
#but how do we know whether the annotation has been assined to the correct row
#thats where we use dplyr using join method
annotated_nhbe1= left_join(filtered_nhbe1, annotation_nbhe1,
                           by= c('ensmblID'='ensembl_gene_id'))
View(annotated_nhbe1)
vol_plot=ggplot(annotated_nhbe1,aes(x=log2FoldChange, y= -log10(padj), 
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
vol_plot
#### changing the plot into interactive plot!!!
#install('plotly')
ggplotly(vol_plot)

#heatmaps
annot_df1_nhbe=annotated_nhbe1[annotated_nhbe1$padj < 0.05,]
annot_df2_nhbe=annot_df1_nhbe[abs(annot_df1_nhbe$log2FoldChange) >1,]
degs=annot_df2_nhbe$ensmblID
vst_nhbe_hm = varianceStabilizingTransformation(dds_nhbe)
vst_nhbe_matrix=assay(vst_nhbe_hm)
View(vst_nhbe_matrix)
data_hm_nhbe=vst_nhbe_matrix[degs,]
rownames(data_hm_nhbe)=annot_df2_nhbe$external_gene_name
heatmap(data_hm_nhbe)
#lets get nice pheatmap instead
pheatmap(data_hm_nhbe)

pheatmap(data_hm_nhbe, scale = 'row', fontsize_row = 4)
heatmap_RdPu=colorRampPalette(brewer.pal(n=9, name = 'Spectral'))(100)
pheatmap(data_hm_nhbe, scale = 'row', fontsize_row = 4, color = heatmap_RdPu, cutree_cols = 2, cutree_rows = 2 )

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

entgene_id_nhbe=getBM(attributes = 'entrezgene_id',
            filters = 'ensembl_gene_id', 
            values = annot_df2_nhbe$ensmblID,
            mart = ensmbl107)
entgene_id_nhbe=entgene_id_nhbe$entrezgene_id
#creating the variable for universe argument which is the whole gene expression pool ####
ent_unv_nhbe=getBM(attributes = 'entrezgene_id',
              filters = 'ensembl_gene_id', values = annotated_nhbe1$ensmblID, mart = ensmbl107)
ent_unv_nhbe=ent_unv_nhbe$entrezgene_id
ent_unv_nhbe= as.character(ent_unv_nhbe)

ego_nhbe= enrichGO(gene = entgene_id_nhbe,
              OrgDb = org.Hs.eg.db,
              ont = 'BP', ######BP=biological process check
              universe = ent_unv_nhbe)
summary(ego_nhbe)
View(summary(ego_nhbe))
barplot(ego_nhbe, showCategory = 15)
dotplot(ego_nhbe)
dotplot(ego_nhbe, showCategory = 15)
cnetplot(ego_nhbe,showCategory = 10)
goplot(ego_nhbe)

####KEGG Enrichment ######
kegg_nhbe=enrichKEGG(gene = entgene_id_nhbe,
                universe = ent_unv_nhbe)
View(summary(kegg_nhbe))
dotplot(kegg_nhbe)
barplot(kegg_nhbe)
dotplot(kegg_nhbe, showCategory = 15)
##########################writing result table to file #########################
write_tsv(annot_df2_nhbe, 'nhbe_DEG.txt')

################################################################################
##############Repeat The analysis with the Other CELL_LINE######################
################################################################################



