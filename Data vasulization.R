library(tidyverse)
library(ggtern)
library(hrbrthemes)
library(ggtext)
library(ggplot2)
library(ggpubr)
library(patchwork)
#install.packages("ggrepel")
#install.packages("matrixStats")
library(matrixStats)
library(ggrepel)
require(pals)

####RNA-seq
setwd("D:/胚乳/差异表达基因")
filePath <- "D:/LYQ/RNA/RNA-seq/v1/"
s2c <- read.table("D:/LYQ/RNA/RNA-seq/design_matrix.txt", header = TRUE, sep='\t',stringsAsFactors=FALSE)
sampleNames <- s2c[,1]
countData.list <- sapply(sampleNames, function(x) read.table(file=paste0(filePath, x, "/abundance.tsv"), header=T, sep="\t"), simplify=F)
countData.df <- do.call("cbind", countData.list)
colsToKeep <- c(1,grep("est_count", names(countData.df)))
ct <- countData.df[,colsToKeep]
names(ct) <- c("transcript_id", sampleNames)
ct[,2:13] <- round(ct[,2:13])
write.csv(ct,"count_matrix.csv")
colsToKeep <- c(1,grep("tpm", names(countData.df)))
ct <- countData.df[,colsToKeep]
names(ct) <- c("transcript_id", sampleNames)
write.csv(ct,"TPM_matrix.csv",row.names=F)
base_dir <- "D:/LYQ/RNA/RNA-seq/v1/"  #R1
sample_id <- dir(file.path(base_dir))
sample_id
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs
s2c <- read.table("D:/LYQ/RNA/RNA-seq/design_matrix.txt", header = TRUE, sep='\t',stringsAsFactors=FALSE)
s2c
s2c <- dplyr::mutate(s2c, path = kal_dirs[-13])
print(s2c)
t2g <- read.table("D:/ZLH/t2gv1.txt", header = TRUE, stringsAsFactors=FALSE)
library(sleuth)
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, num_cores = 1,extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE,gene_mode=T,aggregation_column = 'gene')
#DEG
library(stringr)
a <- list(c("7DAF","14DAF"))
for (xxx in a) {
  s2b <- dplyr::filter(s2c, condition == xxx[1] | condition == xxx[2])
  s2o <- sleuth_prep(s2b, ~ condition, 
                     num_cores = 1,target_mapping = t2g,
                     aggregation_column = 'gene', 
                     extra_bootstrap_summary=TRUE,
                     read_bootstrap_tpm=TRUE,
                     gene_mode = TRUE)
  s2o2 <- sleuth_fit(s2o)
  s2o2 <- sleuth_fit(s2o2, formula = ~ 1, fit_name = "reduced")
  s2o_lrt <- sleuth_lrt(s2o2, "reduced", "full")
  models(s2o_lrt)
  sleuth_table <- sleuth_results(s2o_lrt, 'reduced:full', 'lrt', show_all = FALSE)
  table(sleuth_table[,"qval"] < 0.05)
  sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
  head(sleuth_significant)
  write.csv(sleuth_table,str_c(xxx[1], "vs", xxx[2],"_sleuth_gene_level.csv"),row.names=TRUE,quote=TRUE)
  write.csv(sleuth_significant,str_c(xxx[1], "vs", xxx[2],"_sleuth_significant_gene_level.csv"),row.names=TRUE,quote=TRUE)
  sleuth_matrix <- sleuth_to_matrix(s2o_lrt, 'obs_norm', 'tpm')
  head(sleuth_matrix)
  write.csv(sleuth_matrix, str_c(xxx[1],"vs", xxx[2],"_sleuth_tpm_norm_gene_level.csv"),row.names=TRUE,quote=TRUE)
}

####express gene venn graph
setwd("D:/胚乳/表达的基因")
inter1<-read.table("venn4_inter+un0815.csv",header=T,row.names=1,sep ="\t")
inter<-inter1
inter<-inter[,6]
inter<-strsplit(inter,', ')
totalc <- data.frame()
for (i in 1:15) {
  c <- as.data.frame(inter[[i]])
  colnames(c) <- "X"
  c$type <- rep(paste0("c",i),nrow(c))
  totalc <- rbind(totalc,c)
}
library(ggnewscale)
library(clusterProfiler)
library(enrichplot)
library(org.Twheat.eg.db)
for(i in 1:15)
{
  m<-totalc[which(totalc$type==paste0("c",i)),]
  cg1<-enrichGO(unique(m$X), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  #cg1 <- simplify(cg1,cutoff=0.01,by="p.adjust",select_fun=min)
  #cg1 <- simplify(cg1,cutoff=0.7,by="p.adjust",select_fun=min)  #只有BP CC MF单个才能去冗余
  
  write.csv(cg1,paste0("c",i,"_GO.csv"))
}
totag <- data.frame()
for (i in 1:15){
  fe <- read.csv(paste0("c",i,"_GO.csv"),header=T)
  fe$type <- rep(paste("C",i),nrow(fe))
  totag <- rbind(totag,fe)
}
head(totag)
write.csv(totag,"GO_c15_total.csv")
data <- totag[,c(3,4,8,10,11)]
data <- separate(data = data, col = GeneRatio, into = c("value1", "value2"), sep = "/")
data$value <- as.numeric(data$value1)/as.numeric(data$value2)
data <- data[,c(-2,-3)]
head(data)
top<- data %>% group_by(type) %>% top_n(n = 5, wt = -log10(qvalue))
top <- top[order(top$type),]
top$Description <- factor(top$Description,levels=unique(top$Description),ordered=TRUE)
data$Description <- factor(data$Description,levels=unique(data$Description),ordered=TRUE)
pdf("GO_c15.pdf")
ggplot(top,aes(x=type,y=Description))+
  geom_point(aes(size=`Count`,
                 color=`qvalue`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient(low="pink",high="blue") +
  labs(x=NULL,y=NULL)
dev.off()
write.csv(data,"GO_c15_total.csv")
write.csv(top,"GO_top14.csv")

##14 DAP appear disapper
app <- c(inter[[6]],inter[[7]],inter[[8]],inter[[5]])
dis <- c(inter[[6]],inter[[4]],inter[[3]],inter[[7]])
cg1<-enrichGO(unique(app), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
              pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cg2<-enrichGO(unique(dis), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
              pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
write.csv(cg1,"14_appear.csv")
write.csv(cg2,"14_disappear.csv")
##4 DAP disappear
dis4 <- c(inter[[15]],inter[[7]],inter[[5]],inter[[13]])
cg3<-enrichGO(unique(dis4), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
              pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
write.csv(cg3,"4_disappear.csv")

##7 DAP appear disappear
app7 <- c(inter[[12]],inter[[4]],inter[[2]],inter[[10]])
dis7 <- c(inter[[12]],inter[[11]],inter[[9]],inter[[10]])
cg4<-enrichGO(unique(app7), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
              pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cg5<-enrichGO(unique(dis7), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
              pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
write.csv(cg4,"7_appear.csv")
write.csv(cg5,"7_disappear.csv")

##18 DAP appear disappear
app18 <- c(inter[[14]],inter[[10]],inter[[9]],inter[[13]])
cg6<-enrichGO(unique(app18), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
              pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
write.csv(cg6,"18_appear.csv")
name <- c("4_disappear","14_appear","14_disappear","7_appear","7_disappear","18_appear")
totag <- data.frame()
for (i in name){
  fe <- read.csv(paste0(i,".csv"),header=T)
  fe$type <- rep(paste("C",i),nrow(fe))
  totag <- rbind(totag,fe)
}
head(totag)
data <- totag[,c(3,4,8,10,11)]
data <- separate(data = data, col = GeneRatio, into = c("value1", "value2"), sep = "/")
data$value <- as.numeric(data$value1)/as.numeric(data$value2)
data <- data[,c(-2,-3)]
head(data)
top<- data %>% group_by(type) %>% top_n(n = 5, wt = -log10(qvalue))
top <- top[order(top$type),]
top$Description <- factor(top$Description,levels=unique(top$Description),ordered=TRUE)
data$Description <- factor(data$Description,levels=unique(data$Description),ordered=TRUE)
pdf("GO_appear_disappear.pdf")
ggplot(top,aes(x=type,y=Description))+
  geom_point(aes(size=`Count`,
                 color=`qvalue`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient(low="orange",high="yellow") +
  labs(x=NULL,y=NULL)
dev.off()

####Heatmap for clustering of temporal-specific gene expressions
setwd("D:/胚乳/时序变化基因")
sx <- read.csv("dynamic_gene_expression_filter.csv",row.names=1)
head(sx)
sx <- sx[,-5]
scaled_mat = t(scale(t(sx)))
pdf("dynamic_time_heatmap.pdf")
pheatmap(scaled_mat,
         #annotation_row = annotation_row, 
         cluster_row = F, show_rownames = F, 
         #scale="row",
         cluster_col = FALSE,color = c("#A8C0FF","white","red"))
dev.off()

#########GO
setwd("D:/胚乳/时序变化基因")
library(ggnewscale)
library(clusterProfiler)
library(enrichplot)
library(org.Twheat.eg.db)
for(i in 1:4)
{
  m<-sx[which(sx$cluster==i),]
  cg1<-enrichGO(unique(row.names(m)), OrgDb=org.Twheat.eg.db,ont='BP', keyType="GID",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  #cg1 <- simplify(cg1,cutoff=0.01,by="p.adjust",select_fun=min)
  #cg1 <- simplify(cg1,cutoff=0.7,by="p.adjust",select_fun=min)  #只有BP CC MF单个才能去冗余
  
  write.csv(cg1,paste0("cluster",i,"_GO.csv"))
}
totag <- data.frame()
for (i in 1:4){
  fe <- read.csv(paste0("cluster",i,"_GO.csv"),header=T)
  fe$type <- rep(paste0("Cluster",i),nrow(fe))
  totag <- rbind(totag,fe)
}
head(totag)
data <- totag[,c(3,4,5,8,10,11)]
data <- separate(data = data, col = GeneRatio, into = c("value1", "value2"), sep = "/")
data$value <- as.numeric(data$value1)/as.numeric(data$value2)
head(data)
#top5,即p.adjust最小的ko
top<- data %>% group_by(type) %>% top_n(n = 10, wt = -log10(qvalue))
top <- top[order(top$type),]
top$Description <- factor(top$Description,levels=unique(top$Description),ordered=TRUE)
data$Description <- factor(data$Description,levels=unique(data$Description),ordered=TRUE)
#G0 气泡图
pdf("时序变化4个clusterGOR1.pdf")
ggplot(top,aes(x=type,y=Description))+
  geom_point(aes(size=`Count`,
                 color=`qvalue`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient(low="darkorange",high="goldenrod1") +
  labs(x=NULL,y=NULL)
dev.off()

####KEGG
head(sx)
term2gene <- read.csv("D:/小麦信息/KEGG/KO_Infov1.1.csv",header=T)
term2name <- term2gene[,c(2,3)]
term2gene <- term2gene[,c(2,1)] 
colnames(term2name) <- c("ko_term","Description")
colnames(term2gene) <- c("ko_term","gene_ID")
head(term2gene)
head(term2name)
for(i in 1:4){
  genelist <- row.names(sx[which(sx$cluster == i),])
  gene <- as.factor(genelist)
  fe <- enricher(gene = gene,TERM2GENE = term2gene,TERM2NAME =term2name)
  write.csv(fe,file = paste0("cluster",i,"_KEGGR1.csv"))
}
tota <- data.frame()
for (i in 1:4){
  fe <- read.csv(paste0("cluster",i,"_KEGGR1.csv"),header=T,row.names=1)
  fe$type <- rep(paste0("cluster",i),nrow(fe))
  tota <- rbind(tota,fe)
}
head(tota)
library(tidyverse)
data <- tota[,c(2,6,9,10,3)]
head(data)
data <- separate(data = data, col = GeneRatio, into = c("value1", "value2"), sep = "/")
data$value <- as.numeric(data$value1)/as.numeric(data$value2)
head(data)
data <- data[,c(-5,-6)]
head(data)
top<- data %>% group_by(type) %>% top_n(n = 10, wt = -log10(p.adjust))
top <- top[order(top$type),]
top$Description <- factor(top$Description,levels=unique(top$Description),ordered=TRUE)
pdf("KEGG_top10R1.pdf")
ggplot(top,aes(x=type,y=Description))+
  geom_point(aes(size=`Count`,
                 color=`p.adjust`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="#3F2B96",high="#A8C0FF")+
  labs(x=NULL,y=NULL)
dev.off()
write.csv(data,"KEGG_total.csv")
write.csv(top,"KEGG_top10.csv")

###表达量与组蛋白修饰相关性热图
#####与RNA-seq进行相关性分析
library(ggplot2)
library(dplyr)
library(viridis) # 使用viridis提供的翠绿色标度：scale_fill_viridis()
library(ggpointdensity) # 绘制密度散点图
library(cowplot) #组图
###将所有数据合并chip rna-seq
setwd("D:/胚乳/ChIP/bed/定量")
rna<-read.csv("D:/胚乳/表达的基因/ave_tpm.csv",row.names=1)
head(rna)
rna$gene <- rownames(rna)
#导入chip数据
#先只做gene区域
#先做一个框架
c<-read.csv("0709_18DAF-H3K9ac_gene_peak_tpms.csv",row.names=1)
c <- c[which(c$Tyep == "intergeric" | c$Tyep == "downstream"),]
cc<-as.data.frame(cbind(c$Gene,c$TPM))
colnames(cc)<-c("gene","18DAF-H3K9ac")
cc$"18DAF-H3K9ac"<-as.numeric(cc$"18DAF-H3K9ac")
cc<-tapply(cc$"18DAF-H3K9ac",INDEX=cc$gene,FUN=mean)   #求相同基因的平均peak value
cc<-as.data.frame(cc)
colnames(cc)<-"18DAF-H3K9ac"
cc<-merge(rna,cc,by="row.names",all.x=TRUE)
head(cc)
#循环往里面加数据
name <- c("4DAF-H3K9ac","7DAF-H3K9ac","14DAF-H3K9ac",
          "4DAF-H3K4me3","7DAF-H3K4me3","14DAF-H3K4me3","18DAF-H3K4me3",
          "4DAF-H3K27me3","7DAF-H3K27me3","14DAF-H3K27me3","18DAF-H3K27me3")
for (i in name){
  x<-read.csv(paste0("0709_",i,"_gene_peak_tpms.csv"),row.names=1)
  x <- x[which(c$Tyep == "intergeric" | c$Tyep == "downstream"),]
  head(x)
  xx<-as.data.frame(cbind(x$Gene,x$TPM),stringsAsFactors=F)
  colnames(xx)<-c("gene",i)
  xx[,2]<-as.numeric(xx[,2])
  xx<-tapply(xx[,2],INDEX=xx[,1],FUN=mean)
  xx<-as.data.frame(xx)
  colnames(xx)<-i
  xx$gene<-row.names(xx)
  cc<-merge(cc,xx,by="gene",all.x=TRUE)
}
#保存数据
write.csv(cc,"RNA_chip_intergenericCorrelation_data.csv")
#cc <- read.csv("RNA_chip_genebodyCorrelation_data.csv",row.names=1)
rownames(cc)<-cc$gene
head(cc)
cc<-cc[,c(-1,-2)]
cc[is.na(cc)] <- 0
#绘制chip-atac-rna相关性热图
library('corrplot')
library(stringr)
#cc<-read.csv("RNA-chip-atac-peak-gene-valueR1.csv",header=T, row.names)
ms <- as.data.frame(colnames(cc))
#ac = data.frame(group=str_split(ms,'', simplify = T)[,1])
rownames(ms) = colnames(cc)
M=cor(log2(cc+1), use= "complete.obs", method = 'spearman')
pdf("Corrlation between RNA-chip-intergeneric.pdf",width=12)
pheatmap::pheatmap(M, annotation_col = ms,display_numbers=T)
#pheatmap::pheatmap(M,display_numbers=F)
dev.off()

###k43-only, k27-only, k9ac, k43-k27, k43-k9ac, k9ac-k273, k9ac-k273-k43的基因的表达量箱式图
cd /public/home/chaohe/ChIP/endosperm
module load BEDTools/2.27
awk '{print $0"\t"0}' 4DAF-H3K27me3.peaks.bed > 4DAF-H3K27me3R1.peaks.bed
awk '{print $0"\t"0}' 7DAF-H3K27me3.peaks.bed > 7DAF-H3K27me3R1.peaks.bed
awk '{print $0"\t"0}' 14DAF-H3K27me3.peaks.bed > 14DAF-H3K27me3R1.peaks.bed
awk '{print $0"\t"0}' 18DAF-H3K27me3.peaks.bed > 18DAF-H3K27me3R1.peaks.bed
cat 4DAF-H3K4me3.peaks.bed 4DAF-H3K27me3R1.peaks.bed 4DAF-H3K9ac.peaks.bed | sort -k1,1 -k2n,2 > 4dap.bed
bedtools merge -i - -c 4 -o collapse 4dap.bed > 4dap_merged.peaks.bed
cat 7DAF-H3K4me3.peaks.bed 7DAF-H3K27me3R1.peaks.bed 7DAF-H3K9ac.peaks.bed | sort -k1,1 -k2n,2 > 7dap.bed
bedtools merge -i - -c 4 -o collapse 7dap.bed > 7dap_merged.peaks.bed
cat 14DAF-H3K4me3.peaks.bed 14DAF-H3K27me3R1.peaks.bed 14DAF-H3K9ac.peaks.bed | sort -k1,1 -k2n,2 > 14dap.bed
bedtools merge -i - -c 4 -o collapse 14dap.bed > 14dap_merged.peaks.bed
cat 18DAF-H3K4me3.peaks.bed 18DAF-H3K27me3R1.peaks.bed 18DAF-H3K9ac.peaks.bed | sort -k1,1 -k2n,2 > 18dap.bed
bedtools merge -i - -c 4 -o collapse 18dap.bed > 18dap_merged.peaks.bed

###correlation bewteen DEG and DEP
setwd("D:/胚乳/ChIP/bed/定量/")
deg <- read.csv("D:/胚乳/差异表达基因/7DAFvs14DAF_sleuth_log2FC.csv")
head(deg)
deg <- deg[,c(1,3,4)]
j=1
plot_list=list()
name <- c("H3K9ac","H3K4me3")
for (i in name) {
  d7 <- read.csv(paste0("0709_7DAF-",i,"_gene_peak_tpms.csv"))
  d14 <- read.csv(paste0("0709_14DAF-",i,"_gene_peak_tpms.csv"))
  d7 <- d7[which(d7$Tyep == "promoter"),c(4,15)]
  d14 <- d14[which(d14$Tyep == "promoter"),c(4,15)]
  d7 <- aggregate(d7$TPM,list(d7$Gene),sum)
  d14 <- aggregate(d14$TPM,list(d14$Gene),sum)
  dt <- merge(d7,d14,by="Group.1",all.x=T)
  head(dt)
  #将NA替换为0
  dt[is.na(dt)] <- 1
  #求差异倍数
  dt$log2FC <- log2(dt$x.y/dt$x.x)
  colnames(dt) <- c("target_id","TPM1","TPM2","log2FC")
  head(dt)
  #合并差异表达基因
  deg_his <- merge(deg,dt,by="target_id",all.x=F)
  head(deg_his)
  #画图
  up <- deg_his[which((deg_his$b > 1 & deg_his$log2FC > 1) | (deg_his$b < -1 & deg_his$log2FC < -1)),]
  down <- deg_his[which((deg_his$b > 1 & deg_his$log2FC < -1) | (deg_his$b < -1 & deg_his$log2FC > 1)),]
  up$type <- rep("Up",nrow(up))
  down$type <- rep("Down",nrow(down))
  x <- list(c(up$target_id,down$target_id))
  y<- list(deg_his$target_id)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "target_id"
  no <- merge(deg_his,mm,by="target_id",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  no[no=="#VALUE!"] <- 0
  ###将no分为3类，no1是rna倍数未知，no2是hisone倍数未知，no3是都已知
  tn <- unique(rbind(up,down,no))
  #up <-  tn[which(tn$RNA_log2FC > 1 & tn$histone_log2FC > 1),]
  #down <- tn[which(tn$RNA_log2FC < -1 & tn$histone_log2FC < -1),]
  #up$type <- rep("Up",nrow(up))
  #down$type <- rep("Down",nrow(down))
  #up <- tn[which(tn$RNA_log2FC > 1 & tn$histone_log2FC < -1),]
  #down <- tn[which(tn$RNA_log2FC < -1 & tn$histone_log2FC > 1),]
  up$type <- rep("positive",nrow(up))
  down$type <- rep("negative ",nrow(down))
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no) 
  ttt$b <- as.numeric(ttt$b)
  this_tile <- paste0('positive: ',nrow(up),',negative: ',nrow(down),",R value is ",round(cor.test(ttt$b,ttt$log2FC)$estimate[[1]],digits=3))
  #绘图
  head(ttt)
  library(ggpointdensity)
  p1 <- ggplot(data = ttt, mapping = aes(x = b, y = log2FC)) + 
    #xlim(0,4)+ylim(0,40)+theme_bw()+theme(panel.grid=element_blank()) +
    geom_pointdensity() + #密度散点图（geom_pointdensity???
    #geom_point(aes(colour = factor(context))) +
    scale_color_viridis() + 
    #scale_colour_manual(values=c("#4E5BF4","lightgrey","lightgrey","#E7298A")) +
    geom_smooth(method = lm,se = F,color='red',size=1) +  ##省略拟合曲线
    stat_cor(method = "pearson") + 
    ggtitle("7 DAP vs 14 DAP") +
    theme(text = element_text(family = ,face='bold'),
          axis.text = element_text(family = ,size = 12,face = 'bold'),
          axis.ticks.length=unit(-0.22, "cm"), 
          #加宽图边???
          #panel.border = element_rect(size=1),
          axis.line = element_line(size = .8),
          axis.ticks = element_line(size = .8),
          #去除图例标题
          #legend.title = element_blank(),
          #设置刻度label的边???
          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    labs(x ='expression_log2FC',y= paste0(i,"_log2FC")) +  
    geom_hline(aes(yintercept=0), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed")
plot_list[[j]] <- p1
j=j+1
}
#H3K27me3
i="H3K27me3"
d7 <- read.csv(paste0("0709_7DAF-",i,"_gene_peak_tpms.csv"))
d14 <- read.csv(paste0("0709_14DAF-",i,"_gene_peak_tpms.csv"))
d7 <- d7[which(d7$Tyep == "promoter" | d7$Tyep == "genebody" ),c(4,15)]
d14 <- d14[which(d14$Tyep == "promoter" | d14$Tyep == "genebody"),c(4,15)]
d7 <- aggregate(d7$TPM,list(d7$Gene),sum)
d14 <- aggregate(d14$TPM,list(d14$Gene),sum)
dt <- merge(d7,d14,by="Group.1",all.x=T)
head(dt)
#将NA替换为0
dt[is.na(dt)] <- 1
#求差异倍数
dt$log2FC <- log2(dt$x.y/dt$x.x)
colnames(dt) <- c("target_id","TPM1","TPM2","log2FC")
head(dt)
#合并差异表达基因
head(deg)
deg_his <- merge(deg,dt,by="target_id",all.x=F)
head(deg_his)
deg_his[deg_his=="#VALUE!"] <- 0
deg_his$b <- as.numeric(deg_his$b)
#画图
up <- deg_his[which(((deg_his$b > 1) & (deg_his$log2FC) > 1) | ((deg_his$b < -1) & (deg_his$log2FC < -1))),]
down <- deg_his[which((deg_his$b > 1 & deg_his$log2FC < -1) | (deg_his$b < -1 & deg_his$log2FC > 1)),]
up$type <- rep("Positive",nrow(up))
down$type <- rep("Negativ",nrow(down))
x <- list(c(up$target_id,down$target_id))
y<- list(deg_his$target_id)
mm <- as.data.frame(setDT(y)[!x, on = names(y)])
colnames(mm) <- "target_id"
no <- merge(deg_his,mm,by="target_id",is.all=FALSE)
no$type <- rep("other",nrow(no))
no[no=="#VALUE!"] <- 0
ttt <- rbind(up,down,no) 
ttt$b <- as.numeric(ttt$b)
this_tile <- paste0('positive: ',nrow(up),',negative: ',nrow(down),",R value is ",round(cor.test(ttt$b,ttt$log2FC)$estimate[[1]],digits=3))
#绘图
head(ttt)
#library(ggpointdensity)
p2 <- ggplot(data = ttt, mapping = aes(x = b, y = log2FC)) + 
  #xlim(0,4)+ylim(0,40)+theme_bw()+theme(panel.grid=element_blank()) +
  geom_pointdensity() + #密度散点图（geom_pointdensity???
  #geom_point(aes(colour = factor(context))) +
  scale_color_viridis() + 
  #scale_colour_manual(values=c("#4E5BF4","lightgrey","lightgrey","#E7298A")) +
  geom_smooth(method = lm,se = F,color='red',size=1) +  ##省略拟合曲线
  stat_cor(method = "pearson") + 
  ggtitle("7 DAP vs 14 DAP") +
  theme(text = element_text(family = ,face='bold'),
        axis.text = element_text(family = ,size = 8,face = 'bold'),
        axis.ticks.length=unit(-0.22, "cm"), 
        #加宽图边???
        #panel.border = element_rect(size=1),
        axis.line = element_line(size = .4),
        axis.ticks = element_line(size = .4),
        #去除图例标题
        #legend.title = element_blank(),
        #设置刻度label的边???
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
  labs(x ='expression_log2FC',y= paste0(i,"_log2FC")) +  
  geom_hline(aes(yintercept=0), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed")
#组图
grid.arrange(plot_list[[1]],plot_list[[2]],p2,
             nrow=1,ncol=3)     %>%  ggsave("D7VS14D_deg_dep.pdf",.,width=300,height=100, units="mm")

#####another visualization
setwd("D:/胚乳/ChIP/bed/定量/")
deg <- read.csv("D:/胚乳/差异表达基因/7DAFvs14DAF_sleuth_log2FC.csv")
head(deg)
deg <- deg[,c(1,3,4)]
j=1
plot_list=list()
name <- c("H3K9ac","H3K4me3")
for (i in name) {
  d7 <- read.csv(paste0("0709_7DAF-",i,"_gene_peak_tpms.csv"))
  d14 <- read.csv(paste0("0709_14DAF-",i,"_gene_peak_tpms.csv"))
  d7 <- d7[which(d7$Tyep == "promoter"),c(4,15)]
  d14 <- d14[which(d14$Tyep == "promoter"),c(4,15)]
  d7 <- aggregate(d7$TPM,list(d7$Gene),sum)
  d14 <- aggregate(d14$TPM,list(d14$Gene),sum)
  dt <- merge(d7,d14,by="Group.1",all.x=T)
  head(dt)
  #将NA替换为0
  dt[is.na(dt)] <- 1
  #求差异倍数
  dt$log2FC <- log2(dt$x.y/dt$x.x)
  colnames(dt) <- c("target_id","TPM1","TPM2","log2FC")
  head(dt)
  #合并差异表达基因
  deg_his <- merge(deg,dt,by="target_id",all.x=F)
  deg_his[deg_his=="#VALUE!"] <- 0
  deg_his$b <- as.numeric(deg_his$b)
  head(deg_his)
  #画图
  up <- deg_his[which((deg_his$b > 1 & deg_his$log2FC > 1) | (deg_his$b < -1 & deg_his$log2FC < -1)),]
  down <- deg_his[which((deg_his$b > 1 & deg_his$log2FC < -1) | (deg_his$b < -1 & deg_his$log2FC > 1)),]
  up$type <- rep("Up",nrow(up))
  down$type <- rep("Down",nrow(down))
  x <- list(c(up$target_id,down$target_id))
  y<- list(deg_his$target_id)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "target_id"
  no <- merge(deg_his,mm,by="target_id",is.all=FALSE)
  no$type <- rep("other",nrow(no))
  no[no=="#VALUE!"] <- 0
  ###将no分为3类，no1是rna倍数未知，no2是hisone倍数未知，no3是都已知
  tn <- unique(rbind(up,down,no))
  #up <-  tn[which(tn$RNA_log2FC > 1 & tn$histone_log2FC > 1),]
  #down <- tn[which(tn$RNA_log2FC < -1 & tn$histone_log2FC < -1),]
  #up$type <- rep("Up",nrow(up))
  #down$type <- rep("Down",nrow(down))
  #up <- tn[which(tn$RNA_log2FC > 1 & tn$histone_log2FC < -1),]
  #down <- tn[which(tn$RNA_log2FC < -1 & tn$histone_log2FC > 1),]
  up$type <- rep("positive",nrow(up))
  down$type <- rep("negative ",nrow(down))
  no$type <- rep("other",nrow(no))
  ttt <- rbind(up,down,no) 
  ttt$b <- as.numeric(ttt$b)
  this_tile <- paste0('positive: ',nrow(up),',negative: ',nrow(down),",R value is ",round(cor.test(ttt$b,ttt$log2FC)$estimate[[1]],digits=3))
  #绘图
  head(ttt)
  library(ggpointdensity)
  p1 <- ggplot(ttt, aes(b, log2FC)) + 
    geom_point(aes(colour = factor(type))) +
    #geom_point() +
    scale_colour_manual(values=c("#4E5BF4","lightgrey","#E7298A")) +
    theme(text = element_text(family = ,face='bold'),
          axis.text = element_text(size = 12,face = 'bold'),
          axis.ticks.length=unit(-0.22, "cm"), 
          #加宽图边???
          #panel.border = element_rect(size=1),
          axis.line = element_line(size = .8),
          #axis.ticks = element_line(size = .8),
          #去除图例标题
          legend.title = element_blank(),
          #设置刻度label的边???
          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    theme_bw() + 
    labs(x ="log2_gene_expression_FC",y=paste0("log2_",i,"FC"), title = this_tile) +
    geom_hline(aes(yintercept=1), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=1), colour="#BB0000", linetype="dashed") +
    geom_hline(aes(yintercept=-1), colour="#BB0000", linetype="dashed") +
    geom_vline(aes(xintercept=-1), colour="#BB0000", linetype="dashed")
  plot_list[[j]] <- p1
  j=j+1
}
#H3K27me3
i="H3K27me3"
d7 <- read.csv(paste0("0709_7DAF-",i,"_gene_peak_tpms.csv"))
d14 <- read.csv(paste0("0709_14DAF-",i,"_gene_peak_tpms.csv"))
d7 <- d7[which(d7$Tyep == "promoter" | d7$Tyep == "genebody" ),c(4,15)]
d14 <- d14[which(d14$Tyep == "promoter" | d14$Tyep == "genebody"),c(4,15)]
d7 <- aggregate(d7$TPM,list(d7$Gene),sum)
d14 <- aggregate(d14$TPM,list(d14$Gene),sum)
dt <- merge(d7,d14,by="Group.1",all.x=T)
head(dt)
#将NA替换为0
dt[is.na(dt)] <- 1
#求差异倍数
dt$log2FC <- log2(dt$x.y/dt$x.x)
colnames(dt) <- c("target_id","TPM1","TPM2","log2FC")
head(dt)
#合并差异表达基因
head(deg)
deg_his <- merge(deg,dt,by="target_id",all.x=F)
head(deg_his)
deg_his[deg_his=="#VALUE!"] <- 0
deg_his$b <- as.numeric(deg_his$b)
#画图
up <- deg_his[which(((deg_his$b > 1) & (deg_his$log2FC) > 1) | ((deg_his$b < -1) & (deg_his$log2FC < -1))),]
down <- deg_his[which((deg_his$b > 1 & deg_his$log2FC < -1) | (deg_his$b < -1 & deg_his$log2FC > 1)),]
up$type <- rep("Positive",nrow(up))
down$type <- rep("Negativ",nrow(down))
x <- list(c(up$target_id,down$target_id))
y<- list(deg_his$target_id)
mm <- as.data.frame(setDT(y)[!x, on = names(y)])
colnames(mm) <- "target_id"
no <- merge(deg_his,mm,by="target_id",is.all=FALSE)
no$type <- rep("other",nrow(no))
no[no=="#VALUE!"] <- 0
ttt <- rbind(up,down,no) 
ttt$b <- as.numeric(ttt$b)
this_tile <- paste0('positive: ',nrow(up),',negative: ',nrow(down),",R value is ",round(cor.test(ttt$b,ttt$log2FC)$estimate[[1]],digits=3))
#绘图
head(ttt)
p3<- ggplot(ttt, aes(b, log2FC)) + 
  geom_point(aes(colour = factor(type))) +
  #geom_point() +
  scale_colour_manual(values=c("#4E5BF4","lightgrey","#E7298A")) +
  theme(text = element_text(family = ,face='bold'),
        axis.text = element_text(size = 12,face = 'bold'),
        axis.ticks.length=unit(-0.22, "cm"), 
        #加宽图边???
        #panel.border = element_rect(size=1),
        axis.line = element_line(size = .8),
        #axis.ticks = element_line(size = .8),
        #去除图例标题
        legend.title = element_blank(),
        #设置刻度label的边???
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
  theme_bw() + 
  labs(x ="log2_gene_expression_FC",y=paste0("log2_",i,"FC"), title = this_tile) +
  geom_hline(aes(yintercept=1), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=1), colour="#BB0000", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=-1), colour="#BB0000", linetype="dashed")
#组图
grid.arrange(plot_list[[1]],plot_list[[2]],p3,
             nrow=1,ncol=3)     %>%  ggsave("D7VS14D_deg_dep_number.pdf",.,width=300,height=100, units="mm")

####生成TGT分析文件
name <- c("H3K9ac","H3K4me3")
his <- data.frame()
for (i in name) {
  d7 <- read.csv(paste0("0709_7DAF-",i,"_gene_peak_tpms.csv"))
  d14 <- read.csv(paste0("0709_14DAF-",i,"_gene_peak_tpms.csv"))
  d7 <- d7[which(d7$Tyep == "promoter"),c(4,15)]
  d14 <- d14[which(d14$Tyep == "promoter"),c(4,15)]
  d7 <- aggregate(d7$TPM,list(d7$Gene),sum)
  d14 <- aggregate(d14$TPM,list(d14$Gene),sum)
  dt <- merge(d7,d14,by="Group.1",all.x=T)
  head(dt)
  #将NA替换为0
  dt[is.na(dt)] <- 1
  #求差异倍数
  dt$log2FC <- log2(dt$x.y/dt$x.x)
  colnames(dt) <- c("target_id","TPM1","TPM2","log2FC")
  head(dt)
  #合并差异表达基因
  deg_his <- merge(deg,dt,by="target_id",all.x=F)
  deg_his[deg_his=="#VALUE!"] <- 0
  deg_his$b <- as.numeric(deg_his$b)
  head(deg_his)
  #画图
  up <- deg_his[which((deg_his$b > 1 & deg_his$log2FC > 1) | (deg_his$b < -1 & deg_his$log2FC < -1)),]
  down <- deg_his[which((deg_his$b > 1 & deg_his$log2FC < -1) | (deg_his$b < -1 & deg_his$log2FC > 1)),]
  up$type <- rep(paste0(i,"_positive"),nrow(up))
  down$type <- rep(paste0(i,"_negitive"),nrow(down))
  x <- list(c(up$target_id,down$target_id))
  y<- list(deg_his$target_id)
  mm <- as.data.frame(setDT(y)[!x, on = names(y)])
  colnames(mm) <- "target_id"
  no <- merge(deg_his,mm,by="target_id",is.all=FALSE)
  no$type <- rep(paste0(i,"_none"),nrow(no))
  no[no=="#VALUE!"] <- 0
  ta <- dplyr::bind_rows(as.data.frame(up$target_id))
  colnames(ta) <- c(paste0(i,"_positive"))
  his <- dplyr::bind_rows(his,ta)
}

i="H3K27me3"
d7 <- read.csv(paste0("0709_7DAF-",i,"_gene_peak_tpms.csv"))
d14 <- read.csv(paste0("0709_14DAF-",i,"_gene_peak_tpms.csv"))
d7 <- d7[which(d7$Tyep == "promoter" | d7$Tyep == "genebody" ),c(4,15)]
d14 <- d14[which(d14$Tyep == "promoter" | d14$Tyep == "genebody"),c(4,15)]
d7 <- aggregate(d7$TPM,list(d7$Gene),sum)
d14 <- aggregate(d14$TPM,list(d14$Gene),sum)
dt <- merge(d7,d14,by="Group.1",all.x=T)
head(dt)
#将NA替换为0
dt[is.na(dt)] <- 1
#求差异倍数
dt$log2FC <- log2(dt$x.y/dt$x.x)
colnames(dt) <- c("target_id","TPM1","TPM2","log2FC")
head(dt)
#合并差异表达基因
head(deg)
deg_his <- merge(deg,dt,by="target_id",all.x=F)
head(deg_his)
deg_his[deg_his=="#VALUE!"] <- 0
deg_his$b <- as.numeric(deg_his$b)
#画图
up <- deg_his[which(((deg_his$b > 1) & (deg_his$log2FC) > 1) | ((deg_his$b < -1) & (deg_his$log2FC < -1))),]
down <- deg_his[which((deg_his$b > 1 & deg_his$log2FC < -1) | (deg_his$b < -1 & deg_his$log2FC > 1)),]
up$type <- rep(paste0(i,"_positive"),nrow(up))
down$type <- rep(paste0(i,"_negitive"),nrow(down))
x <- list(c(up$target_id,down$target_id))
y<- list(deg_his$target_id)
mm <- as.data.frame(setDT(y)[!x, on = names(y)])
colnames(mm) <- "target_id"
no <- merge(deg_his,mm,by="target_id",is.all=FALSE)
no$type <- rep("other",nrow(no))
no[no=="#VALUE!"] <- 0
ta <- dplyr::bind_rows(as.data.frame(down$target_id))
colnames(ta) <- c(paste0(i,"_negitive"))
his <- dplyr::bind_rows(his,ta)
write.csv(his,"TGT_histone_regulation.csv")

#####TGT GO注释气泡图
setwd("D:/胚乳/ChIP/bed/定量")
go <- read.csv("TGT_GOTable20221028151152.csv")
#go <- go[which(go$Sample == "Balance_unblance_unblance" | go$Sample == "Balance_unblance_Balance" | go$Sample == "unalance_unblance_Balance" | go$Sample == "Balance_unblance_Balance") ,]
go <- go[,c(1,3,4,7)]
go <- separate(data = go, col = Ratio.in.foreground.list, into = c("value1", "value2"), sep = "/")
go$value <- as.numeric(go$value1)/as.numeric(go$value2)
go$value1 <- as.numeric(go$value1)
go$value2 <- as.numeric(go$value2)
go$Description <- factor(go$Description,levels=unique(go$Description),ordered=TRUE)
head(go)
top<- go %>% group_by(Sample) %>% top_n(n = 10, wt = -log10(FDR))
#write.csv(top,"top.csv")
#top <- read.csv("top.csv")
top <- top[order(top$Sample),]  ##跑两次这个和下面两个，点的排队趋势就对了
top$Description <- factor(top$Description,levels=unique(top$Description),ordered=TRUE)
#G0 气泡图.

pdf("GO_histone_modification.pdf",width=10,family="ArialMT",height = 10)
ggplot(top,aes(x=Sample,y=Description))+
  geom_point(aes(size=`value1`,
                 color=`FDR`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient(low="#61C0BA",high="#EEF8B7") +
  labs(x=NULL,y=NULL)
dev.off()
