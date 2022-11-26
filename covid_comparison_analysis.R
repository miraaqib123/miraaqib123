#############comparison analysis bw the cell lines ############################
#####attach Library######
library(tidyverse)
library(clusterProfiler)
source('covid_rnaseq_analysis.R', echo = TRUE)
source('a549cellline_analysis.R', echo = TRUE)
#nhbe first filter results
#annotated_nhbe1
#a549 first filter results
#annotated1_a549
#### common genes in both analysis#######
nhbe_ensgene=annotated_nhbe1$ensmblID
a549_ensgene=annotated1_a549$ensmblID
common_genes=intersect(nhbe_ensgene, a549_ensgene)
################fetching ensmbl_ids, log2FoldChange #################
nhbe_common_genes=annotated_nhbe1[annotated_nhbe1$ensmblID %in% common_genes, c('ensmblID', 'log2FoldChange')]
a549_common_genes=annotated1_a549[annotated1_a549$ensmblID %in% common_genes, c('ensmblID', 'log2FoldChange')]
############getting the commom gene list into a variable#############
common_FC=left_join(nhbe_common_genes, a549_common_genes, by=c('ensmblID'='ensmblID'))
##########plotting the common genes foldchange from 2 analysis #########
ggplot(data=common_FC, aes(x=log2FoldChange.x, y=log2FoldChange.y))+
  geom_point(alpha=0.3)+
  geom_abline(slope = 1, intercept = 0, color='red')+
  xlim(min=-5, max=5)+
  ylim(min=-5, max=5)
cor.test(common_FC$log2FoldChange.x, common_FC$log2FoldChange.y, method = 'spearman')
cor.test(common_FC$log2FoldChange.x, common_FC$log2FoldChange.y, method = 'pearson')

#######REpeating the same for significantly diff expressed genes#####################
####
nhbe_difexp=annot_df1_nhbe$ensmblID
a549_difexp=annot_df1_a549$ensmblID
#####venn diagram
common_difexp=intersect(nhbe_difexp, a549_difexp)
nhbe_only=setdiff(nhbe_difexp, a549_difexp)
a549_only=setdiff(a549_difexp, nhbe_difexp)

##installing VennDiagram#######
#BiocManager::install('VennDiagram')
#ibrary(VennDiagram)
plot.new()
draw.pairwise.venn(area1 = length(nhbe_difexp), 
                   area2 = length(a549_difexp), 
                   cross.area = length(common_difexp),
                   scaled = TRUE, fill = c('Purple','blue'),
                   aplha=0.5)

######heatmaps#######
union_difexp=union(nhbe_difexp, a549_difexp)
union_fc_hm=filter(common_FC, ensmblID %in% union_difexp)
union_fc_hm_matrix=as.matrix(union_fc_hm[ ,2:3])

my_color=colorRampPalette(rev(brewer.pal(n=7, name = 'RdYlBu')))(100)
my_breaks=c(seq(-2, -0.01, length.out=50),
            0,seq(0.01, 2, length.out=50))
pheatmap(union_fc_hm_matrix, breaks = my_breaks, color = my_color)

