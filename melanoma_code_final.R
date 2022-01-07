library(AUCell)
library(data.table)
library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(rvcheck)
library(ggplot2)
library(openxlsx)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggnewscale)
library(monocle)
library(multcomp)
library(Rmisc)
library(readxl)
library(corrplot)
library(survival)
library(survminer)
library(limma)
library(edgeR)
library(NMF)
library(ggalluvial)

##Data 
melanoma_raw_data <- read_csv("melanoma_raw_data.csv")
rownames(melanoma_raw_data)<-melanoma_raw_data$Tumor
celltype<-melanoma_raw_data[1,]
celltype<-celltype[,-1]
celltype2<-data.frame(t(celltype))
melanoma<-melanoma_raw_data[-1,]
rownames<-melanoma$Tumor
melanoma<-melanoma[,-1]
rownames(melanoma)<-rownames
mydata <- CreateSeuratObject(counts = melanoma, min.cells = 3, project = 'mydata_scRNAseq')
mydata@meta.data$celltype <-celltype2
mydata[["RNA"]]@meta.features<- data.frame(row.names = rownames(mydata[["RNA"]]))

##Normalize
mydata <- NormalizeData(mydata, verbose = FALSE)
mydata<- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)

##Filter
mydata[["percent.mito"]] <- PercentageFeatureSet(mydata, pattern = "^MT",assay = 'RNA')
mito.genes <- grep(pattern = "^MT", x = rownames(x = mydata@assays$RNA), 
                   value = TRUE)
percent.mito <- Matrix::colSums( mydata@assays$RNA[mito.genes, ])/
  Matrix::colSums(mydata@assays$RNA)

mydata<- AddMetaData(object = mydata, metadata = percent.mito,
                     col.name = "percent.mito")
VlnPlot(object = mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)#绘制小提琴图
mydata <-subset(mydata, subset =  percent.mito>0&percent.mito < 0.1&nCount_RNA<20000&nFeature_RNA< 10000&nFeature_RNA> 100)
VlnPlot(object = mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)#绘制小提琴图

##RunPCA
combined_filter <- ScaleData(mydata, vars.to.regress = "percent.mito")
combined_filter <- RunPCA(combined_filter, npcs = 30, verbose = FALSE)
ElbowPlot(combined_filter)

##RunUMAP
combined_filter <- RunUMAP(combined_filter, reduction = "pca", dims = 1:10)
combined_filter <- FindNeighbors(combined_filter, reduction = "pca", dims = 1:10)
combined_filter <- FindClusters(combined_filter, resolution = 0.4)
combined_filter_marker<- FindAllMarkers(combined_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(combined_filter_marker, file = "cluster_marker.xlsx",rowNames=TRUE)

DimPlot(combined_filter,cols =palette1,reduction = "umap",label=TRUE)
table(combined_filter$seurat_clusters)

##Marker
palette1<-c('#8dd3c7','#fee391','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
marker<-c("CD3","CD8A","CD8B","GZMA","CCL5","GZMK","CD69","GPR171","ZFP36L2","KLRB1","IL7R","PRF1","KLRD1","CCL3","GNLY","CD4","CD19","HLA-DRA","C1QA","C1QB","C1QC",'GZMB','GPR183','IRF8','IRF7',"IGKC","IGLC2","JCHAIN","IGHG3",'MYL9','ACTA2','TIMP3',"ENPP2","ESM1","HBA1","SAA1","HILPDA","ANGPTL4","NDUFA4L2","MITF")
DotPlot(combined_filter,features = marker,cols = c('#f6eff7','#1c9099'),assay="RNA")+RotatedAxis()
combined_filter<- RenameIdents(combined_filter, "0"="CD8","1"="CD4","2"="B","3"="Malignant1","4"="Malignant2","5"="Malignant3","6"="Malignant4","7"="MacrophageDC","8"="Malignant5","9"="CD8","10"="EC","11"='CAF')
levels(combined_filter)<-c("CD8","CD4","B","MacrophageDC","EC","CAF","Malignant1","Malignant2","Malignant3","Malignant4","Malignant5")

##group
group<-Idents(combined_filter)
table(group)
group<-gsub("CD8", "non-malignant", group)
group<-gsub("CD4", "non-malignant", group)
group<-gsub("B", "non-malignant", group)
group<-gsub("MacrophageDC", "non-malignant", group)
group<-gsub("EC", "non-malignant", group)
group<-gsub("CAF", "non-malignant", group)
group<-gsub("Malignant1", "malignant", group)
group<-gsub("Malignant2", "malignant", group)
group<-gsub("Malignant3", "malignant", group)
group<-gsub("Malignant4", "malignant", group)
group<-gsub("Malignant5", "malignant", group)

combined_filter@meta.data$"type"=group
table(combined_filter$type)
save(combined_filter,file="combined_filter.Rda")

combined_filter<- RenameIdents(combined_filter, "0"="CD8","1"="CD4","2"="B","3"="Malignant","4"="Malignant","5"="Malignant","6"="Malignant","7"="MacrophageDC","8"="Malignant","9"="CD8","10"="EC","11"='CAF')
levels(combined_filter)<-c("CD8","CD4","B","MacrophageDC","EC","CAF","Malignant")
table(combined_filter@active.ident)
save(combined_filter,file="combined_filter.Rda")
load("combined_filter.Rda")

##Pathway score
WP_FERROPTOSIS<-getGmt("/home/syh/R_project/WP_FERROPTOSIS.gmt")
WP_FERROPTOSIS_gene<-read.gmt("/home/syh/R_project/WP_FERROPTOSIS.gmt")
WP_FERROPTOSIS_gene<-WP_FERROPTOSIS_gene[,2]

HYPOXIA<-getGmt("/data/scRNA_mix/zrh_mix/HALLMARK_HYPOXIA.gmt")
HYPOXIA_gene<-read.gmt("/data/scRNA_mix/zrh_mix/HALLMARK_HYPOXIA.gmt")
HYPOXIA_gene<-HYPOXIA_gene[,2]

oxidative<-getGmt("/data/scRNA_mix/zrh_mix/GO_RESPONSE_TO_OXIDATIVE_STRESS.gmt")
oxidative_gene<-read.gmt("/data/scRNA_mix/zrh_mix/GO_RESPONSE_TO_OXIDATIVE_STRESS.gmt")
oxidative_gene<-oxidative_gene[,2]

GLYCOLYSIS<-getGmt("/home/syh/R_project/HALLMARK_GLYCOLYSIS.gmt")
GLYCOLYSIS_gene<-read.gmt("/home/syh/R_project/HALLMARK_GLYCOLYSIS.gmt")
GLYCOLYSIS_gene<-GLYCOLYSIS_gene[,2]

OXPHOS<-getGmt("/home/syh/R_project/HALLMARK_OXIDATIVE_PHOSPHORYLATION.gmt")
OXPHOS_gene<-read.gmt("/home/syh/R_project/HALLMARK_OXIDATIVE_PHOSPHORYLATION.gmt")
OXPHOS_gene<-OXPHOS_gene[,2]

metastasis<-getGmt("WINNEPENNINCKX_MELANOMA_METASTASIS_UP.gmt")
metastasis_gene<-read.gmt("WINNEPENNINCKX_MELANOMA_METASTASIS_UP.gmt")
metastasis_gene<-metastasis_gene[,2]

HYPOXIA_gene2<-as.data.frame(HYPOXIA_gene)
names(HYPOXIA_gene2)<-c("gene")
WP_FERROPTOSIS_gene2<-as.data.frame(WP_FERROPTOSIS_gene)
names(WP_FERROPTOSIS_gene2)<-c("gene")
oxidative_gene2<-as.data.frame(oxidative_gene)
names(oxidative_gene2)<-c("gene")
GLYCOLYSIS_gene2<-as.data.frame(GLYCOLYSIS_gene)
names(GLYCOLYSIS_gene2)<-c("gene")
OXPHOS_gene2<-as.data.frame(OXPHOS_gene)
names(OXPHOS_gene2)<-c("gene")
metastasis_gene2<-as.data.frame(metastasis_gene)
names(metastasis_gene2)<-c("gene")

gene<-rbind(HYPOXIA_gene2,WP_FERROPTOSIS_gene2,oxidative_gene2,GLYCOLYSIS_gene2,OXPHOS_gene2,metastasis_gene2)
gene<-gene[,1]

exp<-subset(combined_filter,features=gene)
exp<-as.matrix(exp@assays$RNA@data)
cells_rankings_scale<- AUCell_buildRankings(exp, nCores=7, plotStats=TRUE) 

genesets <- list(HYPOXIA=HYPOXIA_gene,
                 FERROPTOSIS=WP_FERROPTOSIS_gene,
                 GLYCOLYSIS=GLYCOLYSIS_gene,
                 OXPHOS=OXPHOS_gene,
                 metastasis=metastasis_gene)

cells_AUC <- AUCell_calcAUC(genesets, cells_rankings_scale,nCores = 5,normAUC = TRUE,aucMaxRank = ceiling(1 * nrow(cells_rankings_scale)))
auc<-getAUC(cells_AUC)
auc_assay<-CreateAssayObject(auc)

combined_filter@assays$gene<-auc_assay
combined_filter<-ScaleData(combined_filter,assay ="gene")

##Cell proportion
cell.prop<-as.data.frame(prop.table(table(combined_filter$type,Idents(combined_filter))))
names(cell.prop)<-c("celltype","cluster","proportion")
ggplot(cell.prop,aes(celltype,proportion,fill=cluster))+
  
  geom_bar(stat="identity",position="fill")+
  
  ggtitle("")+
  
  theme_bw()+
  
  theme(axis.ticks.length=unit(0.5,'cm'))+
  
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = palette1)

## malignant and nonmalignant cell
Idents(combined_filter)<-combined_filter$celltype
table(combined_filter$celltype)
malignant<-subset(combined_filter,idents =c("Malignant1","Malignant2","Malignant3","Malignant4","Malignant5"))
malignant1<-subset(combined_filter,idents =c("Malignant1"))
malignant2<-subset(combined_filter,idents =c("Malignant2"))
malignant3<-subset(combined_filter,idents =c("Malignant3"))
malignant4<-subset(combined_filter,idents =c("Malignant4"))
malignant5<-subset(combined_filter,idents =c("Malignant5"))

nonmalignant<-subset(combined_filter,idents =c("CD8","CD4","B","MacrophageDC","EC","CAF"))

malignant_ave<-AverageExpression(malignant,assay="gene") 
cor<-t(malignant_ave[["gene"]])
cor2<-cor(cor)

nonmalignant_avr<-AverageExpression(nonmalignant,assay="gene") 
cor3<-t(nonmalignant_avr[["gene"]])
cor4<-cor(cor3)

col <- colorRampPalette(c("#4477AA","#77AADD", "#FFFFFF","#EE9988","#BB4444"))
corrplot(cor2, method="color", col=col(200),  
         type="lower", order="hclust", 
         addCoef.col = "black", #添加相关系数
         tl.col="black", tl.srt=30, #修改字体
         insig = "blank", #显著性筛选
         diag=FALSE 
)


##monocle
library(monocle)
a<-malignant@meta.data
a<-cbind(rownames(a),a)
names(a)[1]<-c("cellname")
b<-data.frame(malignant@active.ident)
b<-cbind(rownames(b),b)
names(b)<-c("cellname","active.ident")
c<-merge(a,b,by="cellname",sort=FALSE)
rownames(c)<-c$cellname
c$celltype

c<-c%>% 
  dplyr::select(-celltype) 


pd <- new('AnnotatedDataFrame', data = c)
data <- as(as.matrix(malignant@assays$RNA@counts), 'sparseMatrix')
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(data,
                            phenoData = pd,
                            featureData = fd,
                            lowerDetectionLimit = 0.5,
                            expressionFamily = negbinomial.size())
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"

cds<- estimateSizeFactors(cds)
cds<- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(tail(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))

diff_test_res2<- differentialGeneTest(cds,fullModelFormulaStr = "~active.ident") 

write.table(diff_test_res2,file="monocle.DEG.txt",col.names=T,row.names=T,sep="\t",quote=F)
diff_test_res2<-read.table("monocle.DEG.txt",header=T,row.names = 1,sep="\t")

ordering_genes <- row.names (subset(diff_test_res2, qval < 0.01))
oredering2<-row.names(diff_test_res2)[order(diff_test_res2$qval)][1:2000] 

cds <- setOrderingFilter(cds, oredering2)
cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
cds <- orderCells(cds,root_state = 7)

p1<-plot_cell_trajectory(cds, color_by = "active.ident",show_branch_points = FALSE)+ scale_color_manual(values=palette1)
p2<-plot_cell_trajectory(cds, color_by = "State")+ scale_color_manual(values=palette1)
p3<- plot_cell_trajectory(cds, color_by = "Pseudotime")

p1
p2
p3

diff_test_res3<- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.table(diff_test_res3,file="Pseudotime.monocle.DEG.txt",col.names=T,row.names=T,sep="\t",quote=F)

Time_genes <- top_n(diff_test_res3, n = 250, desc(qval)) %>% pull(gene_short_name) %>% as.character()
write.table(Time_genes,file="Time_genes250.txt",col.names=T,row.names=T,sep="\t",quote=F)

##top 250 genes and ferroptosis 
##TF TP53 MAP1LC3A CP
cds_mel<-cds
s.genes <- c("TF","TP53","MAP1LC3A","CP")
cds2<-cds[s.genes,]
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "active.ident", color_by = "State")+ scale_color_manual(values=palette1)
p2 <- plot_genes_violin(cds[s.genes,], grouping = "active.ident", color_by = "Pseudotime")+ scale_color_manual(values=palette1)
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")+ scale_color_manual(values=palette1)
plotc <- p1|p3

to_be_tested <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("TF","TP53","MAP1LC3A","CP")))
cds_subset <- cds[to_be_tested,]


plot_genes_in_pseudotime(cds_subset, color_by = "State")

to_be_tested <- row.names(subset(fData(cds),
                                 gene_short_name %in% c( "ACSL1","ACSL3","ACSL4","ACSL5","ACSL6","ALOX15","ATG5","ATG7",
                                                         "CP","CYBB","FTH1","FTL","FTMT","GCLC","GCLM","GPX4",   
                                                         "GSS","HMOX1","LPCAT3","MAP1LC3A","MAP1LC3B","MAP1LC3C","NCOA4","PCBP1",
                                                         "PCBP2","PRNP","SAT1","SAT2","SLC11A2","SLC39A14","SLC39A8","SLC3A2",  
                                                         "SLC40A1","SLC7A11","STEAP3","TF","TFRC","TP53","VDAC2","VDAC3" )))
cds_subset <- cds[to_be_tested,]
color=colorRampPalette(c("#80b1d3","#bebada",'#cab2d6'))(100)

plot_pseudotime_heatmap(cds_subset,
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T,
                        hmcols=color)

##cellchat
data.input<-  combined_filter@assays$RNA@data
meta <-  combined_filter@meta.data
meta<-cbind(rownames(meta),meta)
names(meta)[1]<-c("cellname")
active.ident<-data.frame(combined_filter@active.ident)
active.ident<-cbind(rownames(active.ident),active.ident)
names(active.ident)[1]<-c("cellname")
meta<-merge(meta,active.ident,by="cellname")
rownames(meta)<-meta$cellname
meta<-meta[,-1]
names(meta)[8]<-c("active.ident")
library(CellChat)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "active.ident")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "active.ident")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# simply use the default 
##unique(CellChatDB$interaction$annotation)
##CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

##Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

##Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat,population.size = TRUE)
cellchat<- filterCommunication(cellchat, min.cells = 5)
##Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
##Calculate the aggregated cell-cell communication network
cellchat<- aggregateNet(cellchat)
cellchat<- netAnalysis_computeCentrality(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

graphics.off()
##Identify signals contributing most to outgoing or incoming signaling of certain cell groups
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 10,
                                         height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 10,
                                         height = 15)
ht1 + ht2

##Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

##Compute the contribution of each ligand-receptor pair 
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(7:11), remove.isolate = FALSE)

#pattern
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,width = 10,
                                          height = 15)

netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,width = 10,
                                          height = 15)

netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(7:11), remove.isolate = FALSE)

cellchat@netP$pathways

pathways.show <- c("MIF") 
vertex.receiver = seq(1,6)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)


##TCGA
logCPM<-read.table('logCPM_exp_duplicatedelete.txt', header=T,sep='\t',
                   fill=T)
logCPM<-data.frame(t(logCPM),stringsAsFactors = F)
logCPM<-cbind(rownames(logCPM),logCPM)
names(logCPM)[1]<-c("X")
library(dplyr)
exp_single<-logCPM %>% 
  dplyr::select(TF,TP53,MAP1LC3A,CP)

exp_single<-cbind(rownames(exp_single),exp_single)
names(exp_single)[1]<-c("sample")

clinical<-read.table('clinical_used.txt', header=T,sep='\t',
                     fill=T,row.names=1)
clinical<-clinical %>% 
  dplyr::select(TCGA_id,followup_daystodeath,status,metastasis_num,stage_num,age_at_index,gender_num)
names(clinical)[1]<-c("sample")

exp_single_clinical<-merge(exp_single,clinical,by="sample")
rownames(exp_single_clinical)<-exp_single_clinical$sample

exp_single_clinical<-exp_single_clinical[,-1]

write.table(exp_single_clinical,file="exp_single_clinical.txt",col.names=T,row.names=T,sep="\t",quote=F)

# coxph
model <-coxph(Surv(followup_daystodeath, status) ~ TF+TP53+MAP1LC3A+CP,data=exp_single_clinical)
summary(model)

exp_single_clinical2<-exp_single_clinical
exp_single_clinical2$score<- (( -0.02367)*(exp_single_clinical2$TF)+(0.03474)*(exp_single_clinical2$TP53)+(-0.05943 )*(exp_single_clinical2$MAP1LC3A)
                              +(-0.09765)*(exp_single_clinical2$CP))

##
cutoffpoint<- surv_cutpoint(exp_single_clinical2,time="followup_daystodeath",event = "status",variables = "score")
plot(cutoffpoint,"score")

phe_er <- surv_categorize(cutoffpoint)
head(phe_er)
fit <- survfit(Surv(followup_daystodeath, status) ~score, data = phe_er)

ggsurvplot(fit, data = phe_er,risk.table =TRUE,check.names = FALSE,pval =TRUE,
           xlab ="Time in days",
           ggtheme =theme_light()) +
  labs(title = "TCGA")

exp_single_clinical3<-exp_single_clinical2
stage<-exp_single_clinical3$stage_num
table(stage)
stage<-gsub("0", "0", stage)
stage<-gsub("1", "0", stage)
stage<-gsub("2", "1", stage)
stage<-gsub("3", "1", stage)
stage<-gsub("4", "1", stage)
exp_single_clinical3$stage<-stage

exp_single_clinical3$age<-exp_single_clinical3$age_at_index
table(exp_single_clinical3$age)
median(exp_single_clinical3$age_at_index)
exp_single_clinical3$age <- ifelse(exp_single_clinical3$age< median(exp_single_clinical3$age),'0','1')

model <-coxph(Surv(followup_daystodeath, status) ~ stage,data=exp_single_clinical3)
model <-coxph(Surv(followup_daystodeath, status)  ~ metastasis_num,data=exp_single_clinical3)
model <-coxph(Surv(followup_daystodeath, status)  ~ gender_num,data=exp_single_clinical3)
model <-coxph(Surv(followup_daystodeath, status)  ~ age_at_index,data=exp_single_clinical3)
model <-coxph(Surv(followup_daystodeath, status)  ~ score,data=exp_single_clinical3)
model <-coxph(Surv(followup_daystodeath, status)  ~ age,data=exp_single_clinical3)
summary(model)

model <-coxph(Surv(followup_daystodeath, status) ~ stage+age+gender_num+score+metastasis_num,data=exp_single_clinical3)
summary(model)

##GSE59455
gse59455<-read.table('GSE59455_GPL8432_sample_exp.txt', header=T,sep='\t',
                     fill=T,row.names=1)

dge <- DGEList(counts = gse59455)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

gse59455<-data.frame(t(logCPM),stringsAsFactors = F)
gse59455<-cbind(rownames(gse59455),gse59455)
names(gse59455)[1]<-c("X")
gse59455_single<-gse59455 %>% 
  dplyr::select(TF,TP53,MAP1LC3A,CP,MITF)
gse59455_single<-cbind(rownames(gse59455_single),gse59455_single)
names(gse59455_single)[1]<-c("X")

gse59455_clinical<-read.table('GSE59455_sample2.txt', header=T,sep='\t',
                              fill=T,row.names=1)
gse59455_clinical<-cbind(rownames(gse59455_clinical),gse59455_clinical)
names(gse59455_clinical)[1]<-c("X")

exp59455_single_clinical<-merge(gse59455_clinical,gse59455_single,by="X")
rownames(exp59455_single_clinical)<-exp59455_single_clinical$X

exp59455_single_clinical<-exp59455_single_clinical[,-1]
write.table(exp59455_single_clinical,file="exp59455_single_clinical.txt",col.names=T,row.names=T,sep="\t",quote=F)

exp59455_single_clinical<-read.table('exp59455_single_clinical.txt', header=T,sep='\t',
                                     fill=T,row.names=1)

cutoffpoint<- surv_cutpoint(exp59455_single_clinical,time="survival_days",event = "status",variables = "score")
plot(cutoffpoint,"score")

phe_er <- surv_categorize(cutoffpoint)
head(phe_er)
fit <- survfit(Surv(survival_days, status) ~score, data = phe_er)

ggsurvplot(fit, data = phe_er,risk.table =TRUE,check.names = FALSE,pval =TRUE,
           xlab ="Time in days",
           ggtheme =theme_light()) +
  labs(title = "GSE59455")

model_GSE <-coxph(Surv(survival_days, status) ~ TF+TP53+MAP1LC3A+CP,data=exp59455_single_clinical)

summary(model_GSE)

##GSE22153
GSE22153<-read.table('GSE22153.csv', header=T,sep=',',
                     fill=T,row.names=1)
dge <- DGEList(counts = GSE22153)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

GSE22153<-data.frame(t(logCPM),stringsAsFactors = F)
GSE22153<-cbind(rownames(GSE22153),GSE22153)
names(GSE22153)[1]<-c("X")


GSE22153_clinical<-read.table('GSE22153_sample2.txt', header=T,sep='\t',
                              fill=T,row.names=1)
GSE22153_clinical<-cbind(rownames(GSE22153_clinical),GSE22153_clinical)
names(GSE22153_clinical)[1]<-c("X")

GSE22153_single_clinical<-merge(GSE22153_clinical,GSE22153,by="X")
rownames(GSE22153_single_clinical)<-GSE22153_single_clinical$X

GSE22153_single_clinical<-GSE22153_single_clinical[,-1]
write.table(GSE22153_single_clinical,file="GSE22153_single_clinical.txt",col.names=T,row.names=T,sep="\t",quote=F)

summary(model)

GSE22153_single_clinical<-read.table('GSE22153_single_clinical.txt', header=T,sep='\t',
                                     fill=T,row.names=1)

cutoffpoint<- surv_cutpoint(GSE22153_single_clinical,time="survival_days",event = "status",variables = "score")
plot(cutoffpoint,"score")

phe_er <- surv_categorize(cutoffpoint)
head(phe_er)
fit <- survfit(Surv(survival_days, status) ~score, data = phe_er)

ggsurvplot(fit, data = phe_er,risk.table =TRUE,check.names = FALSE,pval =TRUE,
           xlab ="Time in days",
           ggtheme =theme_light()) +
  labs(title = "GSE22153")













