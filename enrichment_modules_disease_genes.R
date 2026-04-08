###compare clustered module of MDD and PTSD
#LOADING libraries
library("clusterProfiler")
library("enrichplot")
library(ggplot2)
library(GOSemSim)
library(stats)
library(DOSE)
library("GO.db")
library(msigdbr)
library(DOSE)
library(forcats)
library(enrichplot)
library(GOplot)
library("ReactomePA")





#orgsnism definition
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#downloading and storing the immune cells marker from Msigdb 
m_df = msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

immuno_cells = msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)




#immune cells GSEA 


mdd = read.csv("E:/Paper1_all_files/common_degs_mdd.csv", header = T)
ids_mdd = bitr(mdd$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

ptsd = read.csv("E:/Paper1_all_files/common_degs_PTSD1.csv", header = T)
ids_ptsd = bitr(ptsd$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)


smoking = read.csv("E:/Paper1_all_files/common_degs_smoking.csv", header = T)
ids_smoking = bitr(smoking$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)


obesity = read.csv("E:/Paper1_all_files/common_degs_obesity.csv", header = T)
ids_obesity = bitr(obesity$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

mdd = mdd[which(mdd$gene %in% ids_mdd$SYMBOL),]
mdd["ENTREZID"]=ids_mdd$ENTREZID
ptsd["ENTREZID"]=ids_ptsd$ENTREZID
smoking["ENTREZID"]=ids_smoking$ENTREZID
obesity["ENTREZID"]=ids_obesity$ENTREZID

mdd_lis = mdd$ENTREZID
smoking_lis = smoking$ENTREZID
ptsd_lis = ptsd$ENTREZID
obesity_lis = obesity$ENTREZID

list = list()
list = append(list,list(mdd_lis) )
list = append(list, list(ptsd_lis))
list = append(list, list(smoking_lis))
list = append(list, list(obesity_lis))

names(list) = c("MDD", "PTSD", "Smoking", "Obesity")


mdd["othegroup"] = "MDD"
ptsd["othegroup"] = "PTSD"
obesity["othegroup"] = "Obesity"
smoking["othegroup"] = "Smoking"





fc_df = rbind(mdd[,4:6], ptsd[,4:6], smoking[,4:6], obesity[,4:6])
fc_df_sort = sort(fc_df$LogFC_avg, decreasing = TRUE)
fc_df$LogFC_avg = fc_df_sort


####################MsigDB
#overreperesntation to find common terms 
msig_clus_enrich= compareCluster(ENTREZID~othegroup,data = fc_df, fun =  enricher, TERM2GENE=immuno_cells,
                                 pvalueCutoff = 1,minGSSize    = 3)
msig_result = msig_clus_enrich@compareClusterResult[msig_clus_enrich@compareClusterResult$pvalue < 0.05 & 
                                                      msig_clus_enrich@compareClusterResult$Count> 2,]
msig_clus_enrich@compareClusterResult = msig_result
noccur = data.frame(table(msig_clus_enrich@compareClusterResult$Description))
noccur = noccur[noccur$Freq >1,]
common_msig_result = msig_clus_enrich@compareClusterResult[which(msig_clus_enrich@compareClusterResult$Description %in%
                                                             noccur$Var1),]
msig_clus_enrich@compareClusterResult = common_msig_result
write.csv(msig_result,"C:/Users/dell/OneDrive - Nile University/Documents/Paper_1/comparecluster_pap1/msigdb_c7_OR.csv")
msig_clus_enrich= setReadable(msig_clus_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dotplot(msig_clus_enrich, showCategory = 15, font.size = 7.5, 
        label_format = 110,color = 'pvalue')+ggtitle("Clustered MDD,PTSD,Obesity,Smoking immunologic gene sets from MsigDB C7")

####GSEA of sets overreperesneted commonly in one mental state versus behavioural phenotype 

#mdd and obesity
mdd_msig= c(mdd$LogFC_avg)
mdd_msig = sort(mdd_msig, decreasing = TRUE)
names(mdd_msig) = c(mdd$ENTREZID)


obes_msig = c(obesity$LogFC_avg)
obes_msig = sort(obes_msig, decreasing = TRUE)
names(obes_msig) = c(obesity$ENTREZID)


smokgsea= c(smoking$LogFC_avg)
smokgsea = sort(smokgsea, decreasing = TRUE)
names(smokgsea) = c(smoking$ENTREZID)

ptsdgsea = c(ptsd$LogFC)
ptsdgsea = sort(ptsdgsea, decreasing = TRUE)
names(ptsdgsea) = c(ptsd$ENTREZID)




mdd_gsea_c7 =GSEA(mdd_msig, TERM2GENE = immuno_cells, pvalueCutoff = 1, minGSSize = 1)
mdd_gsea_c7plot = mdd_gsea_c7@result[mdd_gsea_c7@result$pvalue< 0.05 & mdd_gsea_c7@result$setSize >2 ,]
mdd_gsea_c7@result= mdd_gsea_c7plot[order(mdd_gsea_c7plot$Description),]




obes_gsea_c7 = GSEA(obes_msig,  TERM2GENE = immuno_cells, pvalueCutoff = 1, minGSSize = 1)
obes_gsea_c7plot = obes_gsea_c7@result[obes_gsea_c7@result$pvalue < 0.05 & obes_gsea_c7@result$setSize >2, ]
obes_gsea_c7@result = obes_gsea_c7plot[order(obes_gsea_c7plot$Description),]

inter_obes_mdd = intersect(obes_gsea_c7plot$Description, mdd_gsea_c7plot$Description)


mdd_gsea_c7@result = mdd_gsea_c7plot[which(mdd_gsea_c7plot$Description %in% inter_obes_mdd),]
obes_gsea_c7@result = obes_gsea_c7plot[which(obes_gsea_c7plot$Description %in% inter_obes_mdd),]


r1=ridgeplot(mdd_gsea_c7, fill = "pvalue", label_format = 100, showCategory = 20, orderBy= "ID", decreasing = TRUE)+ggtitle("Enrcihed immunologic gene sets of MDD common with Obesity transcriptome signature ")
r2= ridgeplot(obes_gsea_c7, fill = "pvalue", label_format = 100, showCategory = 20, orderBy= "ID", decreasing = TRUE)+ggtitle("Enriched immunologic gene sets of Obesity common with MDD transcriptome signature ")
cowplot::plot_grid(r1,r2, ncol=1)





##PTSD and smoking 

mdd_gsea_c7 =GSEA(mdd_msig, TERM2GENE = immuno_cells, pvalueCutoff = 1, minGSSize = 1)
mdd_gsea_c7plot = mdd_gsea_c7@result[mdd_gsea_c7@result$pvalue< 0.05 & mdd_gsea_c7@result$setSize>2,]
mdd_gsea_c7@result= mdd_gsea_c7plot[order(mdd_gsea_c7plot$Description),]


smok_msig = c(smoking$LogFC_avg)
smok_msig = sort(smok_msig, decreasing = TRUE)
names(smok_msig) = c(smoking$ENTREZID)

smok_gsea_c7 = GSEA(smok_msig,  TERM2GENE = immuno_cells, pvalueCutoff = 1, minGSSize = 1)
smok_gsea_c7plot = smok_gsea_c7@result[smok_gsea_c7@result$pvalue < 0.05 & mdd_gsea_c7@result$setSize>2, ]
smok_gsea_c7@result = smok_gsea_c7plot[order(smok_gsea_c7plot$Description),]

inter_smok_mdd = intersect(smok_gsea_c7plot$Description, mdd_gsea_c7plot$Description)

ptsd_msig = c(ptsd$LogFC)
ptsd_msig = sort(ptsd_msig, decreasing = TRUE)
names(ptsd_msig) = c(ptsd$ENTREZID)

ptsd_gsea_c7 = GSEA(ptsd_msig,  TERM2GENE = immuno_cells, pvalueCutoff = 1, minGSSize = 3)
ptsd_gsea_c7plot = ptsd_gsea_c7@result[ptsd_gsea_c7@result$pvalue < 0.05 & ptsd_gsea_c7@result$setSize >2, ]
ptsd_gsea_c7@result = ptsd_gsea_c7plot[order(ptsd_gsea_c7plot$Description),]

inter_ptsd_mdd = intersect(ptsd_gsea_c7plot$Description, mdd_gsea_c7plot$Description)

ptsd_gsea_c7@result = ptsd_gsea_c7plot[which(ptsd_gsea_c7plot$Description %in% inter_ptsd_mdd),]
mdd_gsea_c7@result = mdd_gsea_c7plot[which(mdd_gsea_c7plot$Description %in% inter_ptsd_mdd),]


r1=ridgeplot(ptsd_gsea_c7, fill = "pvalue", label_format = 100, showCategory = 15, orderBy= "ID", decreasing = TRUE)+ggtitle("Enrcihed immunologic gene sets of PTSD common with Obesity transcriptome signature ")
r2= ridgeplot(mdd_gsea_c7, fill = "pvalue", label_format = 100, showCategory = 15, orderBy= "ID", decreasing = TRUE)+ggtitle("Enriched immunologic gene sets of Obesity common with PTSD transcriptome signature ")
cowplot::plot_grid(r1,r2, ncol=1)



#########DisGeNet

dgn = compareCluster(ENTREZID~othegroup, data = fc_df, fun = enrichDGN, pvalueCutoff = 1)
result_do = dgn@compareClusterResult[dgn@compareClusterResult$pvalue< 0.05 & dgn@compareClusterResult$Count >2,]
result_do_ptsd = result_do[result_do$Cluster== "PTSD",]
dgn@compareClusterResult = result_do
write.csv(result_do, "C:/Users/dell/OneDrive - Nile University/Documents/Paper_1/comparecluster_pap1/disgenet_cluster.csv")
noccur = data.frame(table(dgn@compareClusterResult$Description))
noccur= noccur[noccur$Freq >1,]
common_result = rbind(dgn@compareClusterResult[which(dgn@compareClusterResult$Description %in%
                                                                   noccur$Var1),], result_do_ptsd)
dgn@compareClusterResult = common_result

dgn = pairwise_termsim(dgn)
emapplot(dgn, showCategory = 40, font.size = 10, color = "pvalue", group_category = TRUE, group_legend = TRUE)


go_bp = compareCluster(list, enrichGO, OrgDb=org.Hs.eg.db, ont = "BP", minGSSize =3, readable=TRUE, pvalueCutoff= 1)
go_bp_res = go_bp@compareClusterResult[go_bp@compareClusterResult$pvalue < 0.05 & 
                                                 go_bp@compareClusterResult$Count >2, ]
go_bp@compareClusterResult = go_bp_res
go_bp@keytype = "ENTREZID"

go_bp = pairwise_termsim(go_bp)
treeplot(go_bp,showCategory = 15, font.size = 2, cex_category = 6,
         offset = rel(5),color ="pvalue", geneClusterPanel = "pie",offset_tiplab = rel(3),
         group_color = c("gray", "orange", "lightblue", "lightgreen", "green")
         )+ggtitle( "Clustering MDD,PTSD,Obesity and Smoking enriched GO terms")



####################
setwd("C:/Users/dell/OneDrive - Nile University/Documents/Paper_1/")

a_reac = gsePathway(ptsdgsea, organism = "human", minGSSize = 1, pvalueCutoff = 1)
result_a_reac = a_reac@result[a_reac@result$pvalue < 0.05 & a_reac@result$setSize >2,]
write.csv(result_a_reac, "reactome_gsea_smok.csv")

a_reac@result = result_a_reac



reac = read.csv('reactome_gsea_mdd.csv', header = T)
a_reac@result = reac
dotplot(a_reac, x = "NES", showCategory = 30, color = "pvalue") + 
  geom_point(aes(shape = Phenotype), size = 3) + 
  scale_shape_manual(values = c(15, 16, 17,18)) 
a_reac = setReadable(a_reac, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
a_reac = pairwise_termsim(a_reac)


d = as.data.frame(a_reac@termsim)
emapplot(a_reac)


library(ggplot2)
ggplot(reac, aes(reac$Description,reac$NES, 
                              color=reac$X, size = setSize)) +
  geom_point(stat = "identity") + 
  scale_color_brewer(type="qual", palette = 6) +
  labs(y = "enricment score", x='')  +
  coord_flip()


gseabpm = gseGO(mdd_msig, ont = "BP", minGSSize = 1, pvalueCutoff = 1, OrgDb = org.Hs.eg.db )
gseabpresultm = gseabpm@result[gseabpm@result$pvalue < 0.05 & gseabpm@result$setSize >2,]
gseabpm@result= gseabpresultm


gseabpp = gseGO(smokgsea, ont = "BP", minGSSize = 1, pvalueCutoff = 1, OrgDb = org.Hs.eg.db )
gseabpresultp = gseabpp@result[gseabpp@result$pvalue < 0.05 & gseabpp@result$setSize >2,]
gseabpp@result= gseabpresultp



gseabp = gseGO(ptsdgsea, ont = "BP", minGSSize = 1, pvalueCutoff = 1, OrgDb = org.Hs.eg.db )
gseabpresult = gseabp@result[gseabp@result$pvalue < 0.05 & gseabp@result$setSize >2,]
gseabp@result= gseabpresult


gseabp = gseGO(obes_msig, ont = "BP", minGSSize = 1, pvalueCutoff = 1, OrgDb = org.Hs.eg.db )
gseabpresult = gseabp@result[gseabp@result$pvalue < 0.05 & gseabp@result$setSize >2,]
gseabp@result= gseabpresult
















####################################################################################################################
#2nd paper 


alldegs = read.csv("E:/Paper1_all_files/common_degs_mdd.csv", header = TRUE)
allptsddegs = read.csv("E:/Paper1_all_files/common_degs_ptsds.csv", header = TRUE)





###getting the entrezid instead of sybmols 
ids_mdd = bitr(alldegs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
ids_ptsd =  bitr(allptsddegs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

mdd = alldegs[which(alldegs$gene %in% ids_mdd$SYMBOL),]
mdd["ENTREZID"] = ids_mdd$ENTREZID
mdd["othergroup"] = "MDD"

allptsddegs = allptsddegs[which(allptsddegs$gene %in% ids_ptsd$SYMBOL),]
allptsddegs["ENTREZID"]= ids_ptsd$ENTREZID
allptsddegs["othergroup"] = "PTSD"

fcdf = rbind(mdd[,3:5],allptsddegs[,3:5])


list = list(ids_mdd$ENTREZID)
list = append(list, list(ids_ptsd$ENTREZID))
names(list)= c("MDD Modules","PTSD Modules")

# mddcomunities = read.delim("E:/MDD_network/mdd_disease_genes.txt", sep = "\t", header = F)
# mdddd = alldegs[which(alldegs$gene %in% mddcomunities$V1),]
logfc_mdd = mdd$LogFC_avg
names(logfc_mdd) = mdd$ENTREZID
logfc_mdd = sort(logfc_mdd, decreasing = T)

logfc_pstd = allptsddegs$LogFC
names(logfc_pstd) = allptsddegs$ENTREZID
logfc_pstd = sort(logfc_pstd, decreasing = T)



###MDD and PTSD modules clustering and Comparsion for GO and immune cells in Msig, PATHWAYS


go_bp = compareCluster(list, enrichGO, OrgDb=org.Hs.eg.db, ont = "BP", minGSSize =3, readable=TRUE, pvalueCutoff= 1)
go_bp_res = go_bp@compareClusterResult[go_bp@compareClusterResult$pvalue < 0.05 & 
                                         go_bp@compareClusterResult$Count >2, ]
go_bp@compareClusterResult = go_bp_res
go_bp@keytype = "ENTREZID"

go_bp = pairwise_termsim(go_bp)
treeplot(go_bp,showCategory = 30, font.size = 2, cex_category = 6,
         offset = rel(5),color ="pvalue", geneClusterPanel = "pie",offset_tiplab = rel(3),
         group_color = c("gray", "orange", "lightblue", "lightgreen", "green")
)+ggtitle( "Clustering MDD,PTSD enriched GO terms")


write.csv(go_bp_res ,"C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/gobp_mdd_ptsd.csv")

# mddcomunities = read.delim("E:/MDD_network/mdd_disease_genes.txt", sep = "\t", header = F)


gseabp_mdd = gseGO(logfc_mdd, ont = "BP", minGSSize = 1, pvalueCutoff = 1, OrgDb = org.Hs.eg.db )
gseabp_mdddf = gseabp_mdd@result[gseabp_mdd@result$pvalue < 0.05 & 
                                         gseabp_mdd@result$setSize >2, ]  

gseabp_mdd@result = gseabp_mdddf
gseabp_mdd = setReadable(gseabp_mdd, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(gseabp_mdd@result ,"C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/gsea_bp_mdd_readable.csv")





####PTSD
gseabp_ptsd = gseGO(logfc_pstd, ont = "BP", minGSSize = 1, pvalueCutoff = 1, OrgDb = org.Hs.eg.db)
gseabp_ptsddf = gseabp_ptsd@result[gseabp_ptsd@result$pvalue < 0.05 & 
                                   gseabp_ptsd@result$setSize >2, ] 
gseabp_ptsd@result = gseabp_ptsddf
gseabp_ptsd = setReadable(gseabp_ptsd, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(gseabp_ptsd@result ,"C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/gsea_bp_ptsd_readable.csv")


interbp= intersect(gseabp_mdddf$Description, gseabp_ptsddf$Description)


gseabp_mdd@keytype = "ENTREZID"

bpmdd = pairwise_termsim(gseabp_mdd)
dev.new()
treeplot(bpmdd,showCategory = 60, font.size = 2, cex_category = 0.8,
         offset = rel(1),color ="NES", geneClusterPanel = "pie",offset_tiplab = rel(3),
         group_color = c("#52688F","#ECD5BB","#BDF0D6", "#FFFFFF", "#710117"))

bpptsd = pairwise_termsim(gseabp_ptsd)

write.csv(gseabp_mdddf,"bp_gseamdd.csv")
write.csv(gseabp_ptsddf,"bp_gseaptsd.csv")





sub = c("gene", "LogFC")
ptsddegs = allptsddegs[,which(colnames(allptsddegs) %in% sub)]
colnames(ptsddegs) = c("ID", "logFC")
write.csv(ptsddegs, "ptsd_degs_for_chord.csv")

ptsd_gse = data.frame(Category= allpaths_enriched$Category,
                       ID = allpaths_enriched$ID, 
                       Term = allpaths_enriched$Description, 
                       Count = allpaths_enriched$setSize,
                       Genes = allpaths_enriched$core_enrichment,
                       adj_pval=as.numeric(allpaths_enriched$pvalue))
ptsd_gse = as.data.frame(lapply(ptsd_gse, function(x) gsub("/", ",", x)))
write.csv(allpaths_enriched, "C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/pathwyas_pstd.csv")



library(GOplot)
library(ggrepel)
library(ggrepel)

mdd_chord = read.csv("C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/mdd_pathways.csv")
mddchord = data.frame(Category= mdd_chord$database,
                      ID = mdd_chord$ID, 
                      Term = mdd_chord$Description, 
                      Count = mdd_chord$setSize,
                      Genes = mdd_chord$core_enrichment,
                      adj_pval=as.numeric(mdd_chord$pvalue))




sub = c("gene", "LogFC_avg")
alldegsmdd = alldegs[,which(colnames(alldegs) %in% sub)]
colnames(alldegsmdd) = c("ID", "logFC")
write.csv(alldegsmdd, "mdd_degs_for_chord.csv")
circ = circle_dat(mddchord,alldegsmdd)

GOCircle(circ, zsc.col = c("green","red", "black"), lfc.col = c("#BDF0D6", "#F652A0")) 

genes = unname(EC@gene2Symbol)
chord = chord_dat(circ, genes, mdd_gse$Term)

# Generate the matrix with a list of selected genes
chord = chord_dat(data = circ, genes = genes)
# Generate the matrix with selected processes
Logfc= alldegs[which(alldegs$gene %in% rownames(chord)),]
Logfc = Logfc[ order(match(Logfc$gene, rownames(chord))), ]

chord = as.data.frame(chord_dat(data = circ, process = mdd_gse$Term))
chord$logFC = Logfc$LogFC_avg

GOChord(chord, gene.order = "logFC",gene.space = 0.4, gene.size = 3, 
        space = 0.01, lfc.col = c("red", "white", "cyan"), nlfc = 1, border.size = 0.1)


GOHeat(chord, nlfc = 1, fill.col = c('red', 'white', 'cyan'))
library(GSEA)
GSEA.barplot(EC@result, category = EC@result$Description, score = EC@result$NES,
             pvalue = EC@result$pvalue, sort = EC@result$NES, decreasing = TRUE)

#######msig: immunologic gene sets enrichement
mdd_gsea_c7 =GSEA(logfc_mdd, TERM2GENE = immuno_cells, pvalueCutoff = 1, minGSSize = 1)
mdd_gsea_c7plot = mdd_gsea_c7@result[mdd_gsea_c7@result$pvalue< 0.05 & mdd_gsea_c7@result$setSize>2,]
mdd_gsea_c7@result= mdd_gsea_c7plot[order(mdd_gsea_c7plot$Description),]
write.csv(mdd_gsea_c7plot ,"C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/gsea_mdd.csv")

ptsd_gsea_c7 = GSEA(logfc_pstd,  TERM2GENE = immuno_cells, pvalueCutoff = 1, minGSSize = 3)
ptsd_gsea_c7plot = ptsd_gsea_c7@result[ptsd_gsea_c7@result$pvalue < 0.05 & ptsd_gsea_c7@result$setSize >2, ]
ptsd_gsea_c7@result = ptsd_gsea_c7plot[order(ptsd_gsea_c7plot$Description),]
write.csv(ptsd_gsea_c7plot , "C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/gsea_ptsd.csv")


inter_ptsd_mdd = intersect(ptsd_gsea_c7plot$Description, mdd_gsea_c7plot$Description)

ptsd_gsea_c7@result = ptsd_gsea_c7plot[which(ptsd_gsea_c7plot$Description %in% inter_ptsd_mdd),]
mdd_gsea_c7@result = mdd_gsea_c7plot[which(mdd_gsea_c7plot$Description %in% inter_ptsd_mdd),]

setwd("C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/")
r1=ridgeplot(ptsd_gsea_c7, fill = "pvalue", label_format = 100, showCategory = 17, orderBy= "ID", decreasing = TRUE)+ggtitle("Enrcihed immunologic gene sets of PTSD common with MDD DEGs")
r2= ridgeplot(mdd_gsea_c7, fill = "pvalue", label_format = 100, showCategory = 17, orderBy= "ID", decreasing = TRUE)+ggtitle("Enriched immunologic gene sets of MDD common with PTSD DEGs")
cowplot::plot_grid(r1,r2, ncol=1)




#for pathways

paths= gsePathway(logfc_pstd,organism = "human", pvalueCutoff = 1, minGSSize=1,pAdjustMethod = "BH")
paths = setReadable(paths, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

pathway_raec_res = paths@result[paths@result$pvalue < 0.05 &
                                                     paths@result$setSize >2, ]
# path_kegg@result = pathway_kegg_res
# path_kegg@keytype = "ENTREZID"
# path_kegg@readable = TRUE
# heatplot(path_kegg)


paths = gseWP(logfc_pstd,organism  = "Homo sapiens", pvalueCutoff = 1, minGSSize=1,pAdjustMethod = "BH")
paths = setReadable(paths, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

pathway_wiki_result = paths@result[paths@result$pvalue < 0.05 &
                                     paths@result$setSize >2, ]

paths = gseKEGG(logfc_pstd,organism  = "hsa", pvalueCutoff = 1, minGSSize=1,pAdjustMethod = "BH")
paths = setReadable(paths, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

pathway_kegg_result = paths@result[paths@result$pvalue < 0.05 &
                                     paths@result$setSize >2, ]


pathway_raec_res$Category = "Reactome"
pathway_kegg_result$Category = "KEGG"
pathway_wiki_result$Category = "Wikipathways"
allpaths_enriched = rbind(pathway_raec_res,pathway_kegg_result,pathway_wiki_result)

# path_kegg@result = pathway_kegg_res
# path_kegg@keytype = "ENTREZID"
# path_kegg@readable = TRUE
# heatplot(path_kegg)
# 
write.csv(allpaths_enriched, "C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/ptsd_pathways.csv")



library(ggplot2)
ggplot(allpaths_enriched, aes(allpaths_enriched$Description,allpaths_enriched$enrichmentScore, 
                             color=allpaths_enriched$database, size = setSize)) +
  geom_point(stat = "identity") + 
  scale_color_brewer(type="qual", palette = 6) +
  labs(y = "enricment score", x='')  +
  coord_flip()

paths@result = allpaths_enriched

write.csv(allpaths_enriched, "C:/Users/dell/OneDrive - Nile University/Documents/Paper_2/pathwyas_pstd.csv")


pathway_rep = compareCluster(ENTREZID~othegroup,data =fc_df, fun = enrichPathway , pvalueCutoff = 1, minGSSize= 1,pAdjustMethod = "BH")
pathway_res = pathway_rep@compareClusterResult[pathway_rep@compareClusterResult$pvalue < 0.05 &
                                                  pathway_rep@compareClusterResult$Count >2, ]
pathway_rep@compareClusterResult = pathway_res
pathway_rep@keytype = "ENTREZID"
pathway_rep@readable = TRUE
path_mdd =enrichPathway(mdd_lis, organism = "human", pvalueCutoff = 0.05, readable = TRUE,pAdjustMethod = "BH", minGSSize = 3)
res_mdd = path_mdd@result[path_mdd@result$pvalue < 0.05 & path_mdd@result$Count >2,]
res_mdd["Cluster"] = "MDD"
res_mdd["othegroup"] = "MDD"
path_mdd@result= res_mdd
pathway_rep@compareClusterResult = rbind(pathway_rep@compareClusterResult, res_mdd)



msig_clus_enrich= setReadable(msig_clus_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dgn_clus_enrich = setReadable(dgn_clus_enrich , OrgDb = org.Hs.eg.db, keyType="ENTREZID")

pathsem = pairwise_termsim(paths)
emapplot(pathsem, showCategory = 30, label_format = 10, cex_label_gene = 0.6)+ggtitle("Comparison between MDD and PTSD GSEA\n for MsigDB immune cells marker ")
dotplot(dgn_clus_enrich, showCategory = 10, font.size = 10,color = "pvalue", by= "GeneRatio")+ggtitle("Comparison between MDD and PTSD DisGeNET GSEA")

write.csv(dgn_clus_enrich@compareClusterResult ,"C:/Users/dell/OneDrive - Nile University/Documents/Paper_1/DISGENET_MDD_PTSD_GSEA.csv")
write.csv(msig_clus_enrich@compareClusterResult ,"C:/Users/dell/OneDrive - Nile University/Documents/Paper_1/msigdb_MDD_PTSD_overrep.csv")





library(GOplot)
setwd("C:/Users/dell/OneDrive - Nile University/Documents/Paper_2")
ec = read.csv("pathwyas_pstd.csv", row.names = 1L)
organism = "org.Hs.eg.db"

edc = setReadable(ec, OrgDb = organism, keyType="ENTREZID")

circ <- circle_dat(EC$david, EC$genelist)

# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = EC$genes)
# Generate the matrix with selected processes
chord <- chord_dat(data = circ, process = EC$process)

GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
GOChord(chord, limit = c(3, 0), gene.order = 'logFC')

























