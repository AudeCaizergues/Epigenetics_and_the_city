#!/usr/bin/Rscript
### PURPOSE OF THE SCRIPT: test for an enrichment of GO

########################################################################
################      GO analyses with topGO       #####################
########################################################################


#### LIBRARIES #########################################################
#BiocManager::install("topGO")
#BiocManager::install("GO.db")
# BiocManager::install("biomaRt")
# BiocManager::install("Rgraphviz")
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
########################################################################


#### STEP 1: DATA PREPARATION ##########################################
setwd("/PATH/TO/06_results/genes_DMR")
exp_data= read.table('bg_genes_total_5kb.txt', header=TRUE)
bg_genes=as.character(exp_data[,1])
length(bg_genes)

# Read in genes of interest
setwd("/PATH/TO/06_results/genes_DMR/5kb")
candidate_list =read.table('liste_genes_5kb_DMR_toutesvilles.txt', header=TRUE)
candidate_list= as.character(candidate_list[,1])
length(candidate_list)

setwd("/PATH/TO/06_results/GO_enrichment analyses/total_5kb")
########################################################################
########################################################################

#### STEP2: GO annotation ##############################################
# create GO db for genes to be used using biomaRt - please note that this takes a while
# ensembl <- useMart("ensembl")
# listDatasets(ensembl)
#dataset interessant pour le GO : le poulet = ggallus_gene_ensembl, le zebra finch = tguttata_gene_ensembl
db= useMart('ENSEMBL_MART_ENSEMBL',dataset='ggallus_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg_genes, mart=db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])

# remove any candidate genes without GO annotation
keep = candidate_list %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes
########################################################################

#### STEP 3: MAKE topGO DATA OBJECT ####################################
GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO,nodeSize=3)
########################################################################

#### STEP 4: TEST FOR SIGNIFICANCE #####################################
# define test using the classic algorithm with fisher (refer to [1] if you want to understand how the different algorithms work)
classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')
# Selecting ‘algorithm=classic’ means that the GO hierarchy isn’t taken into account, 
# so each GO term is tested independently (over-representation enrichment). The limitation 
# of using this is that all genes annotated to a GO terms will be automatically annotated 
# to its parents as well, therefore a GO term might look enriched just because its children 
# are enriched. Thus, it is important that GO hierarchy is taken into account (conditional 
# enrichment) to avoid redundancy.

# define test using the weight01 algorithm (default) with fisher
weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher') 
# generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
allGO=usedGO(GOdata)
all_res=GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))
results.table.p= all_res[which(all_res$weightFisher<=0.05),]
results.table.p
write.table(results.table.p,"summary_topGO_analysis.csv",sep=",",quote=FALSE,row.names=FALSE)
showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "none", sigForAll=FALSE, firstSigNodes=13,.NO.CHAR=50)

########################################################################

# #### STEP 5: CORRECTING FOR MULTIPLE TESTING ###########################
#performing BH correction on our p values
p.adj<-round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_final=cbind(all_res,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]

#get list of significant GO before multiple testing correction
results.table.p= all_res_final[which(all_res_final$weightFisher<=0.05),]

# #get list of significant GO after multiple testing correction
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]

# #save first top 50 ontolgies sorted by adjusted pvalues
 write.table(all_res_final[1:50,],"summary_topGO_analysis_with_pval_correction.csv",sep=",",quote=FALSE,row.names=FALSE)

# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level

#pdf(file='topGOPlot_fullnames_useinfo_all.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "none", sigForAll=FALSE, firstSigNodes=13,.NO.CHAR=50)
#dev.off()

#### STEP 6: GET THE CANDIDATE GENES IN SIGNIFICANT GO TERMS
myterms <- results.table.p$GO.ID
mygenes <- genesInTerm(GOdata, myterms)
var=c()
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% candidate_list # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  var[i]=print(paste("Term",myterm,"genes:",mygenesforterm2))
}

write.table(var,"candidate_genes_in_significant_GO.txt",sep="\t",quote=F)

#### STEP 7: GET ALL (even those that are not in the candidate list) THE GENES IN SIGNIFICANT GO TERMS
myterms =results.table.p$GO.ID 
mygenes = genesInTerm(GOdata, myterms)

variable=c()
for (i in 1:length(myterms))
{
  myterm=myterms[i]
  mygenesforterm= mygenes[myterm][[1]]
  mygenesforterm=paste(mygenesforterm, collapse=',')
  var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}
write.table(variable,"all_genes_in_significant_GO.txt",sep="\t",quote=F)

