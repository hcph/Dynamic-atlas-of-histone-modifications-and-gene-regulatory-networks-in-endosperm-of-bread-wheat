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

####RNA-seq datasets overlap with public embryo RNA-seq datasets
mar <- read.csv("D:/生殖发育/胚的数据/embryo_TPM_average.csv",row.names=1)
colnames(mar) <- paste0("embyro_",colnames(mar))
head(mar)
ct <- read.csv("D:/胚乳/差异表达基因/TPM_matrix.csv",row.names=1)
ct$gene <- row.names(ct)
ct$gene <- gsub("\\.[0-9]","",ct$gene)
ct_tt <- aggregate(ct[,1:12],list(ct$gene),sum)
head(ct_tt)
row.names(ct_tt) <- ct_tt[,1]
ct_tt <- ct_tt[,-1]
library(matrixStats)
ct_tt$DPA4 <-  rowMeans(ct_tt[,7:9])
ct_tt$DPA7 <-  rowMeans(ct_tt[,10:12])
ct_tt$DPA14 <-  rowMeans(ct_tt[,1:3])
ct_tt$DPA18 <-  rowMeans(ct_tt[,4:6])
ct_ave <- ct_tt[,13:16]
head(ct_ave)
colnames(ct_ave) <- paste0("endosperm_",colnames(ct_ave))
#合并
total <- merge(mar,ct_ave,by="row.names",all=F)
head(total)
row.names(total) <- total[,1]
total <- total[,-1]
head(total)
dim(total)
whole <- read.table("L:/wheat/expression/wheat_development_expression.txt",sep="\t",row.names = 1,header=T)
#aa <- whole[,c(104,108)]
#aa <- whole[,c(104,108,119)]
whole[1:4,1:4]
head(aa)
#colnames(aa) <- c("DPA6_starchy_endosperm","DPA9_starchy_endosperm","Dough_embryo_proper","Dough_endosperm")
colnames(aa) <- c("DPA6_starchy_endosperm","DPA9_starchy_endosperm")
total2 <- merge(total,aa,by="row.names",all=F)
row.names(total2) <- total2[,1]
total2 <- total2[,-1]
head(total2)
#PCA
library(FactoMineR)
library(ggsci)
library("ggrepel")
ct <- t(total2)
gene.pca <- PCA(ct, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca) 
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )
group <- read.csv('D:/生殖发育/RNA-seq/groupR1.csv',row.names=1)
group <- group[rownames(pca_sample), ]
pca_sample <- cbind(pca_sample, group)
pca_sample$samples <- rownames(pca_sample)
pca_sample$line <- rep("RNA",nrow(pca_sample))
head(pca_sample)
library(ggplot2)
pdf("D:/胚乳/大修/与肖军胚的数据的overlap/胚乳和胚的PCA聚类图R1.pdf",width=8)
ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2, label = samples, shape = samples)) +
  geom_point(aes(color = samples), size = 3) +  #根据样本坐标绘制二维散点图
  geom_text_repel(aes(label = samples,color=samples),
                  #color = "gray20",
                  #data = subset(pca_sample, samples %in% pointsToLabel),
                  force = 10) +
  scale_shape_manual(values = c(1:60)) +
  scale_color_d3("category20") +
  #scale_linetype_manual(values = c('twodash', 'longdash', 'dashed')) 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中
dev.off()
#heatmap
library(matrixStats)
ac = data.frame(group=str_split(colnames(total2),'_',simplify = T)[,1])
rownames(ac) = colnames(total2)
M=cor(log(total2+1),method = "pearson")
pdf("D:/胚乳/大修/与肖军胚的数据的overlap/heatmap_tpm_spearman.pdf",family="ArialMT")
pheatmap::pheatmap(M, annotation_col = ac,display_numbers=T)
dev.off()

###expression of embryo specific gene
head(total)
spe <- read.csv("D:/胚乳/胚乳0728/specific_gene.csv")
head(spe)
row.names(spe) <- spe[,1]
spe <- spe[,c(1,3)]
spe_total <- merge(spe,total,by="row.names",all=F)
head(spe_total)
row.names(spe_total) <- spe_total[,1]
spe_total <- spe_total[,-1:-3]
library(pheatmap)
head(spe_total)
spe_total <- as.numeric(spe_total)
pdf("D:/胚乳/大修/与肖军胚的数据的overlap/胚和胚乳中特异表达的基因.pdf",family="ArialMT")
pheatmap(spe_total[,1:8],
         #annotation_row = annotation_row, 
         cluster_row = F, show_rownames = T, 
         display_numbers = T,
         #scale="row",
         cluster_cols = F)
dev.off()

###Data on embryo-specific expression in embryo and expression in endosperm
total <- merge(mar,ct_ave,by="row.names",all.x=F)
head(total)
row.names(total) <- total[,1]
total <- total[,-1]
#中位数标准化
library("preprocessCore") 
colname <- colnames(total)
rowname <- row.names(total)
total <- normalize.quantiles(as.matrix(total))
row.names(total) <- rowname
colnames(total) <- colname
head(total)

#Embryo-specific expressed genes, homologous
library(openxlsx)
em <- read.xlsx("D:/胚乳/胚乳0728/胚乳特异基因R1.xlsx",sheet = 4)
head(em)
em <- unique(em)
row.names(em) <- em[,3]
em_total <- merge(em,total,by="row.names",all.x=TRUE)
row.names(em_total) <- paste(em_total[,3],em_total[,1],sep="_")
em_total <- em_total[,5:16]
head(em_total)
em_total <- log2(em_total+1)
library(pheatmap)
pdf("D:/胚乳/大修/与肖军胚的数据的overlap/胚特异表达的基因.pdf",family = "ArialMT")
pheatmap(em_total,
         #annotation_row = annotation_row, 
         cluster_row = T, show_rownames = T, 
         display_numbers = T,
         #scale="row",
         cluster_cols = F)
dev.off()
###marker
#AtS2
em <- unique(em)
at <- em[which(em$symbol == "AtS2"),]
head(at)
row.names(at) <- at[,3]
em_total <- merge(at,total,by="row.names",all.x=TRUE)
row.names(em_total) <- paste(em_total[,3],em_total[,1],sep="_")
em_total <- em_total[,5:16]
head(em_total)
em_total <- log2(em_total+1)
library(pheatmap)
pdf("D:/胚乳/大修/与肖军胚的数据的overlap/AtS2特异表达的基因.pdf",family = "ArialMT")
pheatmap(em_total,
         #annotation_row = annotation_row, 
         cluster_row = T, show_rownames = T, 
         display_numbers = T,
         scale="row",
         cluster_cols = F)
dev.off()
AtS2_171200
AtS2_168900
AtS2_186100
AtS2_367200
AtS2_364100

##Genes specifically expressed in embryo and endosperm
setwd("D:/胚乳/胚乳0728")
gene <- read.xlsx("D:/胚乳/胚乳0728/胚乳特异基因R1.xlsx",sheet = 5)
#gene <- gene[which(gene$symbol != "Zm.66589" & gene$symbol != "Zm.3896"),]
gene <- unique(gene)
row.names(gene) <- gene[,1]
tpm <- read.csv("yq_tpm_norm_gene_levelR1.csv",header=T,row.names=1,check.names = F)
head(tpm)
spe_tpm <- merge(gene,tpm,by="row.names",all.x=F)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm<- spe_tpm %>%
  rownames_to_column(var = 'sample') %>%
  pivot_longer( cols =  c("1_Four_DAF-1":"Dough(embryoproper)_3"),
                names_to = 'stage',
                values_to = 'expr')
spe_tpm$log2expr<-log2(spe_tpm$expr+1)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm$gene_id <- paste(spe_tpm$type,spe_tpm$gene_id,sep="_")
head(spe_tpm)
pdf("markerR1.pdf",family = "ArialMT")
ggplot(spe_tpm,aes(x=stage,y=gene_id))+
  geom_point(aes(size= log2expr,
             color= log2expr)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient2(low="#323996",high="#B12241",mid = "white")
  #scale_color_gradient2(midpoint =mid,low="#323996",high="#B12241",mid = "white")
dev.off()
gene <- read.xlsx("D:/胚乳/胚乳0728/胚乳特异基因R1.xlsx",sheet = 5)
gene <- gene[which(gene$symbol != "Zm.66589" & gene$symbol != "Zm.3896"),]
gene <- unique(gene)
row.names(gene) <- gene[,1]
tpm <- read.csv("yq_tpm_norm_gene_levelR1.csv",header=T,row.names=1,check.names = F)
head(tpm)
spe_tpm <- merge(gene,tpm,by="row.names",all.x=F)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm<- spe_tpm %>%
  rownames_to_column(var = 'sample') %>%
  pivot_longer( cols =  c("1_Four_DAF-1":"Dough(embryoproper)_3"),
                names_to = 'stage',
                values_to = 'expr')
spe_tpm$log2expr<-log2(spe_tpm$expr+1)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm$gene_id <- paste(spe_tpm$type,spe_tpm$gene_id,sep="_")
head(spe_tpm)
#pdf("marker.pdf",family = "ArialMT")
p1 <- ggplot(spe_tpm,aes(x=stage,y=gene_id))+
  geom_point(aes(size= expr,
             color= expr)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient2(low="#323996",high="#B12241",mid = "white")
  #scale_color_gradient2(midpoint =mid,low="#323996",high="#B12241",mid = "white")
  scale_color_gradient2(low = "white", mid = 'white', high = "red")

  gene <- read.xlsx("D:/胚乳/胚乳0728/胚乳特异基因R1.xlsx",sheet = 5)
gene <- gene[which(gene$symbol != "Zm.66589" & gene$symbol != "Zm.3896"),]
gene <- unique(gene)
row.names(gene) <- gene[,1]
tpm <- read.csv("yq_tpm_norm_gene_levelR1.csv",header=T,row.names=1,check.names = F)
head(tpm)
spe_tpm <- merge(gene,tpm,by="row.names",all.x=F)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm<- spe_tpm %>%
  rownames_to_column(var = 'sample') %>%
  pivot_longer( cols =  c("1_Four_DAF-1":"Dough(embryoproper)_3"),
                names_to = 'stage',
                values_to = 'expr')
spe_tpm$log2expr<-log2(spe_tpm$expr+1)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm$gene_id <- paste(spe_tpm$type,spe_tpm$gene_id,sep="_")
head(spe_tpm)
#pdf("marker.pdf",family = "ArialMT")
p1 <- ggplot(spe_tpm,aes(x=stage,y=gene_id))+
  geom_point(aes(size= expr,
             color= expr)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  #scale_color_gradient2(low="#323996",high="#B12241",mid = "white")
  #scale_color_gradient2(midpoint =mid,low="#323996",high="#B12241",mid = "white")
  scale_color_gradient2(low = "white", mid = 'white', high = "red")
##Embryo-specific maize homologous genes drawn separately
gene <- read.xlsx("D:/胚乳/胚乳0728/胚乳特异基因R1.xlsx",sheet = 5)
gene <- gene[which(gene$symbol == "Zm.66589" | gene$symbol == "Zm.3896"),]
gene <- unique(gene)
row.names(gene) <- gene[,1]
tpm <- read.csv("yq_tpm_norm_gene_levelR1.csv",header=T,row.names=1,check.names = F)
head(tpm)
spe_tpm <- merge(gene,tpm,by="row.names",all.x=F)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm<- spe_tpm %>%
  rownames_to_column(var = 'sample') %>%
  pivot_longer( cols =  c("1_Four_DAF-1":"Dough(embryoproper)_3"),
                names_to = 'stage',
                values_to = 'expr')
spe_tpm$log2expr<-log2(spe_tpm$expr+1)
spe_tpm <- spe_tpm[order(spe_tpm$type),]
spe_tpm$gene_id <- paste(spe_tpm$type,spe_tpm$gene_id,sep="_")
head(spe_tpm)
pdf("marker2.pdf",family = "ArialMT")
ggplot(spe_tpm,aes(x=stage,y=gene_id))+
  geom_point(aes(size= expr,
                 color= expr)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient2(low = "white", high = "red")
dev.off()  
#scale_color_gradient2(low="#323996",high="#B12241",mid = "white")
#scale_color_gradient2(midpoint =mid,low="#323996",high="#B12241",mid = "white")

# Define the min-max normalization function
min_max_normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
spe_tpm_normalized <- spe_tpm %>%
  group_by(gene_id) %>%
  mutate(normalized_log2expr = min_max_normalize(log2expr)) %>%
  ungroup()
ggplot(spe_tpm_normalized, aes(x = stage, y = gene_id)) +
  geom_point(aes(size = normalized_log2expr, color = normalized_log2expr)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_color_gradient2(
    low = "white",   # Color for the minimum normalized value
    mid = "white",   # Color for the midpoint of the gradient
    high = "red",    # Color for the maximum normalized value
    midpoint = 0.5   # Midpoint for the normalized data (0.5 in case of min-max normalization)
  ) +
  scale_size_continuous(range = c(2, 10))  # Adjust size scale as needed
pdf("normalized_marker.pdf",family="ArialMT")
ggplot(spe_tpm_normalized, aes(x = stage, y = gene_id)) +
  geom_point(aes(size = normalized_log2expr, color = normalized_log2expr)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_color_gradient2(
    low = "white",   # Color for the minimum normalized value
    high = "red",    # Color for the maximum normalized value
  ) +
  scale_size_continuous(range = c(2, 10))  # Adjust size scale as needed
dev.off()


#####Significantly Differentially Expressed Genes Correlate with Histone Modification Malignancy, 4DAF vs. 14DAF
###Correlation of differentially expressed and differentially changed genes
#Differential Modification and Differential Expression Correlation
#7 vs 14
setwd("D:/胚乳/ChIP/bed/定量/")
deg <- read.csv("D:/胚乳/差异表达基因/7DAFvs14DAF_sleuth_log2FC.csv")
#deg <- deg[which((deg$b > 1 | deg$b < -1) & deg$qval < 0.05),]
head(deg)
deg <- deg[,c(1,3,4)]
j=1
plot_list=list()
name <- c("H3K9ac","H3K4me3")
for (i in name) {
  d7 <- read.csv(paste0("D:/胚乳/ChIP/bed/定量/","0709_7DAF-",i,"_gene_peak_tpms.csv"))
  d14 <- read.csv(paste0("D:/胚乳/ChIP/bed/定量/","0709_14DAF-",i,"_gene_peak_tpms.csv"))
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


#####The correlation between significantly differentially expressed genes and significantly differentially modified peak
##Get the value of each repeat
#count value
cd /public/home/chaohe/endosperm
#H3K27me3
cat /public/home/stbi/task/yqchip/chip1/macs2/final/4DAF-H3K27me3.peaks.bed /public/home/stbi/task/yqchip/chip1/macs2/final/7DAF-H3K27me3.peaks.bed \
/public/home/stbi/task/yqchip/chip1/macs2/final/14DAF-H3K27me3.peaks.bed /public/home/stbi/task/yqchip/chip1/macs2/final/18DAF-H3K27me3.peaks.bed | sort -k1,1 -k2n,2 | bedtools merge -i - >merged_H3K27me3.bed
bsub  -J h3k27me3 -n 2 -o h3k27me3.out -e h3k27me3.err -q smp "bedtools multicov -bams /public/home/stbi/task/yqchip/chip1/align/Q17.final.bam /public/home/stbi/task/yqchip/chip1/align/Q18.final.bam \
/public/home/stbi/task/yqchip/chip1/align/Q19.final.bam /public/home/stbi/task/yqchip/chip1/align/Q20.final.bam \
/public/home/stbi/task/yqchip/chip1/align/Q21.final.bam /public/home/stbi/task/yqchip/chip1/align/Q22.final.bam  \
/public/home/stbi/task/yqchip/chip1/align/Q23.final.bam /public/home/stbi/task/yqchip/chip1/align/Q24.final.bam  \
 -bed  merged_H3K27me3.bed >H3K27me3_mergedpeak_rep_count_total.bed"
#H3K4me3
cat /public/home/stbi/task/yqchip/chip2/macs2/final/4DAF-H3K4me3.peaks.bed /public/home/stbi/task/yqchip/chip2/macs2/final/7DAF-H3K4me3.peaks.bed \
/public/home/stbi/task/yqchip/chip2/macs2/final/14DAF-H3K4me3.peaks.bed /public/home/stbi/task/yqchip/chip2/macs2/final/18DAF-H3K4me3.peaks.bed | sort -k1,1 -k2n,2 | bedtools merge -i - >merged_H3K4me3.bed
bsub  -J h3k4me3 -n 2 -o h3k4e3.out -e h3k4me3.err -q smp "bedtools multicov -bams /public/home/stbi/task/yqchip/chip2/align/Q1.final.bam /public/home/stbi/task/yqchip/chip2/align/Q2.final.bam \
/public/home/stbi/task/yqchip/chip2/align/Q3.final.bam /public/home/stbi/task/yqchip/chip2/align/Q4.final.bam \
/public/home/stbi/task/yqchip/chip2/align/Q5.final.bam /public/home/stbi/task/yqchip/chip2/align/Q6.final.bam  \
/public/home/stbi/task/yqchip/chip2/align/Q11.final.bam /public/home/stbi/task/yqchip/chip2/align/Q12.final.bam  \
-bed  merged_H3K4me3.bed >H3K4me3_mergedpeak_rep_count_total.bed"
#H3K9ac
cat /public/home/stbi/task/yqchip/chip3/macs2/final/4DAF-H3K9ac.peaks.bed /public/home/stbi/task/yqchip/chip3/macs2/final/7DAF-H3K9ac.peaks.bed \
/public/home/stbi/task/yqchip/chip3/macs2/final/14DAF-H3K9ac.peaks.bed /public/home/stbi/task/yqchip/chip3/macs2/final/18DAF-H3K9ac.peaks.bed | sort -k1,1 -k2n,2 | bedtools merge -i - >merged_H3K9ac.bed
bsub  -J h3k9ac -n 2 -o h3k9ac.out -e h3k9ac.err -q smp "bedtools multicov -bams /public/home/stbi/task/yqchip/chip2/align/Q7.final.bam /public/home/stbi/task/yqchip/chip2/align/Q13.final.bam \
/public/home/stbi/task/yqchip/chip2/align/Q8.final.bam /public/home/stbi/task/yqchip/chip2/align/Q14.final.bam \
/public/home/stbi/task/yqchip/chip2/align/Q9.final.bam /public/home/stbi/task/yqchip/chip2/align/Q15.final.bam  \
/public/home/stbi/task/yqchip/chip2/align/Q10.final.bam /public/home/stbi/task/yqchip/chip2/align/Q16.final.bam  \
-bed  merged_H3K9ac.bed >H3K9ac_mergedpeak_rep_count_total.bed"


#####The multiplicity of differences in simultaneous protein modification with differentially expressed genes was 0.75
setwd("D:/胚乳/ChIP/bed/定量/")
###DEP
name <- c("H3K4me3","H3K9ac","H3K27me3")
for (i in name) {
  countdata <- read.table(paste0(i,"_mergedpeak_rep_countR1.bed"), row.names = 4,header=F)
  countdata <- countdata[,-1:-5]
  colnames(countdata) <- c("V0_rep1","V0_rep2","V26N0_rep1","V26N0_rep2","V26N6_rep1","V26N6_rep2")
  head(countdata)
  ## 过滤在所有重复样本中小于1的基因
  countdata = countdata[rowMeans(countdata) > 1,]
  #导入样本注释信息
  coldata  <- read.csv("../deseq2/coldata.csv",row.names = 1)
  coldata$sample <- coldata$condition
  head(coldata)
  #检查数据Counts文件与coldata数据是否匹配
  all(rownames(coldata) %in% colnames(countdata))  
  all(rownames(coldata) == colnames(countdata))
  ##差异分析
  # 制作差异矩阵
  ##V0vsV26N0
  a1 <- countdata[,1:4]
  b1 <- coldata[1:4,]
  head(b1)
  dds <-  DESeqDataSetFromMatrix(countData = a1,colData = b1,design = ~ condition) 
  # 过滤
  dds <- dds[rowSums(counts(dds)) > 1,]  
  nrow(dds)  
  # 差异比较
  dep <- DESeq(dds)
  res <- results(dep)
  diff = res
  diff <- na.omit(diff)  ## 去除NA
  #加上tpm值
  t <- read.csv(paste0(i,"_peak_tpms_pos.csv"),row.names = 1)
  t <- t[,-1:-3]
  head(t)
  t1 <- t[,1:4]
  head(t1)
  diff <- as.data.frame(diff)
  diff_tpm <- merge(diff,t1,by="row.names",all.x=F)
  head(diff_tpm)
  write.csv(diff_tpm,paste0(i,"V0vsV26N0_diff.csv"))  # 导出所有的差异文件
  ##V26N0vsV26N6
  a2 <- countdata[,-1:-2]
  head(a2)
  b2 <- coldata[c(-1,-2),]
  head(b2)
  dds <-  DESeqDataSetFromMatrix(countData = a2,colData = b2,design = ~ condition) 
  # 过滤
  dds <- dds[rowSums(counts(dds)) > 1,]  
  nrow(dds)  
  # 差异比较
  dep <- DESeq(dds)
  res <- results(dep)
  diff = res
  diff <- na.omit(diff)  ## 去除NA
  #加上tpm值
  t2 <- t[,-1:-2]
  diff <- as.data.frame(diff)
  diff_tpm <- merge(diff,t2,by="row.names",all.x=F)
  head(diff_tpm)
  write.csv(diff_tpm,paste0(i,"V26N0vsV26N6_diff.csv"))  # 导出所有的差异文件
  ##V0vsV26N6
  a3 <- countdata[,-3:-4]
  head(a3)
  b3 <- coldata[c(-3,-4),]
  head(b3)
  dds <-  DESeqDataSetFromMatrix(countData = a3,colData = b3,design = ~ condition) 
  # 过滤
  dds <- dds[rowSums(counts(dds)) > 1,]  
  nrow(dds)  
  # 差异比较
  dep <- DESeq(dds)
  res <- results(dep)
  diff = res
  diff <- na.omit(diff)  ## 去除NA
  #加上tpm值
  t3 <- t[,-3:-4]
  diff <- as.data.frame(diff)
  diff_tpm <- merge(diff,t3,by="row.names",all.x=F)
  head(diff_tpm)
  write.csv(diff_tpm,paste0(i,"V0vsV26N6_diff.csv"))  # 导出所有的差异文件
  ##加上位置信息
  pos <- read.csv(paste0(i,"_peak_tpms_pos.csv"))
  a1 <- read.csv(paste0(i,"V0vsV26N0_diff.csv"))
  a1$peak <- a1$Row.names
  b1 <- merge(a1,pos,by="peak",all.x=F)
  b1 <- b1[which(b1$padj < 0.05),]
  write.csv(b1,paste0(i,"V0vsV26N0_sig_DEP.csv"))
  #b1 <- b1[,c(12:14,3)]
  #write.table(b1,paste0(i,"V0vsV26N0_sig_DEP.bed"))
  #V26N0vsV26N6
  pos <- read.csv(paste0(i,"_peak_tpms_pos.csv"))
  a1 <- read.csv(paste0(i,"V26N0vsV26N6_diff.csv"))
  a1$peak <- a1$Row.names
  b2 <- merge(a1,pos,by="peak",all.x=F)
  b2 <- b2[which(b2$padj < 0.05),]
  write.csv(b2,paste0(i,"V26N0vsV26N6_sig_DEP.csv"))
  #b2 <- b2[,c(12:14,3)]
  #write.table(b2,paste0(i,"V26N0vsV26N6_sig_DEP.bed"))
  #V0vsV26N6
  pos <- read.csv(paste0(i,"_peak_tpms_pos.csv"))
  a1 <- read.csv(paste0(i,"V0vsV26N6_diff.csv"))
  a1$peak <- a1$Row.names
  b3 <- merge(a1,pos,by="peak",all.x=F)
  b3 <- b3[which(b3$padj < 0.05),]
  write.csv(b3,paste0(i,"V0vsV26N6_sig_DEP.csv"))
  #b3 <- b3[,c(12:14,3)]
  #write.table(b3,"V0vsV26N6_sig_DEP.bed")
  ###韦恩图看差异
  library(VennDiagram)
  venn_list <- list(V0vsV26N0 = b1$peak, V0vsV26N6 = b2$peak, V26N0vsV26N6 = b3$peak)
  venn.plot<-venn.diagram(venn_list, filename = NULL, 
                          fill = c('#7FC97F', '#BEAED4', '#FDC086'), alpha = 0.50, 
                          cat.col = c('#7FC97F', '#BEAED4', '#FDC086'), cat.cex = 1.5, cat.fontfamily = 'serif',
                          col = c('#7FC97F', '#BEAED4', '#FDC086'), cex = 1.5, fontfamily = 'serif')
  pdf(file=paste0(i,"venn3.pdf"))
  grid.draw(venn.plot)
  dev.off()
  ###生成所有DEP文件
  b1 <- read.csv(paste0(i,"V0vsV26N0_sig_DEP.csv"),row.names=1)
  b2 <- read.csv(paste0(i,"V26N0vsV26N6_sig_DEP.csv"),row.names=1)
  b3 <- read.csv(paste0(i,"V0vsV26N6_sig_DEP.csv"),row.names=1)
  b1$type <- rep("V0vsV26N0",nrow(b1))
  b2$type <- rep("V26N0vsV26N6",nrow(b2))
  b3$type <- rep("V0vsV26N6",nrow(b3))
  colnames(b1) <- colnames(b2)
  colnames(b3) <- colnames(b1)
  #total <- rbind(b1[,-8:-11],b2[,-8:-11],b3[,-8:-11])
  total <- rbind(b1,b2,b3)
  #total <- total[which(total$log2FoldChange >= 0.75 | total$log2FoldChange =< -0.75),]
  #total <- total[,c(8:10,3,6,7,17,1,2)]
  head(total)
  ##绘制热图,用最宽的peak来做
  #按行标准化
  pea <- total[,c(1,17:22)]
  colnames(pea) <- c("peak","TPM1","TPM2","TPM3","TPM4","TPM5","TPM6")
  pea <- unique(pea)
  row.names(pea) <- pea[,1]
  pea <- pea[,-1]
  head(pea)
  scaled_mat = t(scale(t(pea)))
  head(scaled_mat)
  library(circlize)
  p1 <- Heatmap(scaled_mat, name = i, km = 7, 
                #top_annotation = ha, 
                cluster_columns = FALSE,
                #top_annotation_height = unit(4, "mm"),
                show_row_names = FALSE, show_column_names = TRUE,
                col = colorRamp2(c(-2,0,2), c("blue","#EEEEEE", "red")))
  #Heatmap(expr$length, name = "length", width = unit(5, "mm"), col = colorRamp2(c(0, 100000), c("white", "orange"))) +
  #Heatmap(expr$type, name = "type", width = unit(5, "mm")) +
  #Heatmap(expr$chr, name = "chr", width = unit(5, "mm"), col = rand_color(length(unique(expr$chr))))
  png(paste0(i,"_DEG_heatmapRfR1.png"))
  ht = draw(p1)
  dev.off()
  saveRDS(ht,paste0(i,"_ht.rds"))
  pdf(paste0(i,"_DEG_heatmapRfR1.pdf"))
  ht = draw(p1)
  dev.off()
  cluster1<-row_dend(ht)
  set.seed(nrow(scaled_mat))   
  cls <-plyr::ldply(row_order(ht),data.frame)
  names(cls) <- c("id","rowid")
  cls <- dplyr::arrange(cls,rowid)
  pea$cluster <- cls$id
  write.csv(pea,paste0(i,"_row_order.csv"))
}

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
  up <- deg_his[which((deg_his$b > 1 & deg_his$log2FC > 0.75) | (deg_his$b < -1 & deg_his$log2FC < -0.75)),]
  down <- deg_his[which((deg_his$b > 1 & deg_his$log2FC < -0.75) | (deg_his$b < -1 & deg_his$log2FC > 0.75)),]
  #deg_his <- deg_his[which(deg_his$qval < 0.05),]
  #up <- deg_his[which((deg_his$log2FC > 1) | (deg_his$log2FC < -1)),]
  #down <- deg_his[which((deg_his$log2FC < -1) | (deg_his$log2FC > 1)),]
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
  #ttt <- rbind(up,down,no) 
  ttt <- rbind(up,down)
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
p3<- ggplot(data = ttt, mapping = aes(x = b, y = log2FC)) + 
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
  geom_hline(aes(yintercept=1), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=1), colour="#BB0000", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="#BB0000", linetype="dashed") +
  geom_vline(aes(xintercept=-1), colour="#BB0000", linetype="dashed")
#组图
grid.arrange(plot_list[[1]],plot_list[[2]],p3,
             nrow=1,ncol=3)     %>%  ggsave("D7VS14D_deg1_dep0.7_corr.pdf",.,width=300,height=100, units="mm")


####haplotype
cd /public/home/chaohe/prgwas/haptype
awk '{print $1"\tIWGSC_v1.1\tgene\t"$2"\t"$3"\t.\t"$6"\t.\t"$4}' /public/home/chaohe/db/geneR1.bed >geneR1.txt  
perl -p -i -e 's/chr1A/1/g' geneR1.txt
perl -p -i -e 's/chr1B/2/g' geneR1.txt
perl -p -i -e 's/chr1D/3/g' geneR1.txt
perl -p -i -e 's/chr2A/4/g' geneR1.txt
perl -p -i -e 's/chr2B/5/g' geneR1.txt
perl -p -i -e 's/chr2D/6/g' geneR1.txt
perl -p -i -e 's/chr3A/7/g' geneR1.txt
perl -p -i -e 's/chr3B/8/g' geneR1.txt
perl -p -i -e 's/chr3D/9/g' geneR1.txt
perl -p -i -e 's/chr4A/10/g' geneR1.txt
perl -p -i -e 's/chr4B/11/g' geneR1.txt
perl -p -i -e 's/chr4D/12/g' geneR1.txt
perl -p -i -e 's/chr5A/13/g' geneR1.txt
perl -p -i -e 's/chr5B/14/g' geneR1.txt
perl -p -i -e 's/chr5D/15/g' geneR1.txt
perl -p -i -e 's/chr6A/16/g' geneR1.txt
perl -p -i -e 's/chr6B/17/g' geneR1.txt
perl -p -i -e 's/chr6D/18/g' geneR1.txt
perl -p -i -e 's/chr7A/19/g' geneR1.txt
perl -p -i -e 's/chr7B/20/g' geneR1.txt
perl -p -i -e 's/chr7D/21/g' geneR1.txt
perl -p -i -e 's/chrUN/22/g' geneR1.txt
#提取候选基因上有3.5kb和基因内部上的snp
awk '{if($7 == "+") print $1"\t"$4-3500"\t"$5"\t"$9; else if($7 == "-") print $1"\t"$4"\t"$5+3500"\t"$9}' geneR1.txt >target.txt
plink --bfile /public/home/chaohe/prgwas/HC_genotype --extract range target.txt --make-bed --out target_snp
#plink --bfile target_snp --recode vcf-iid --out target_snp_vcf
##做ttest
#计算p-value
#提前准备好基因型文件，gene.txt 表型文件 和  NL18r2.txt
#准备NL18r2.txt
#cd /public/home/tllu/haptype
cd /public/home/chaohe/prgwas/haptype
mkdir tgw
mkdir gw
mkdir gl
mkdir glgw
awk '{print "tgw\t.\t.\t.\t.\t.\t"$9"\t."}' geneR1.txt >tgw/NL18r0.txt
awk '{print "gw\t.\t.\t.\t.\t.\t"$9"\t."}' geneR1.txt >gw/NL18r0.txt
awk '{print "gl\t.\t.\t.\t.\t.\t"$9"\t."}' geneR1.txt >gl/NL18r0.txt
awk '{print "glgw\t.\t.\t.\t.\t.\t"$9"\t."}' geneR1.txt >glgw/NL18r0.txt
#准备表型数据
cp /public/home/chaohe/prgwas/tgw.txt tgw/tgw.txt
cp /public/home/chaohe/prgwas/gw.txt gw/gw.txt
cp /public/home/chaohe/prgwas/gl.txt  gl/gl.txt 
cp /public/home/chaohe/prgwas/glgw.txt glgw/glgw.txt
#运行程序
#cd /public/home/tllu/haptype
cd /public/home/chaohe/prgwas/haptype
module load Python/3.8.6
module load plink/1.9
for i in tgw gw gl glgw;
do
cd /public/home/chaohe/prgwas/haptype/"$i"
bsub  -J "$i" -n 1 -o "$i".out -e "$i".err -q high -R "rusage[mem=80GB]" "python NL18r0.py"
done
###整合所有的单倍型显著的基因
setwd("D:/胚乳/大修/单倍型")
name <- c("GL","GW","GWGL","TGW")
tt_total <- data.frame()
for (i in name) {
  tt <- read.table(paste0(i,"_result.txt"),sep="\t")
  tt <- unique(tt[,c(1,2)])
  tt_total <- rbind(tt_total,tt)
}
#合并相同的基因
hy <- tt_total
colnames(hy) <- c("trait","gene")
hy$trait <- gsub("glgw.txt","GLGW",hy$trait)
hy$trait <- gsub("tgw.txt","TGW",hy$trait)
hy$trait <- gsub("gw.txt","GW",hy$trait)
hy$trait <- gsub("gl.txt","GL",hy$trait)
head(hy)
hy$value <- 1
hy_wide <- spread(hy,key=trait,value=value)
hy_wide[is.na(hy_wide)] <- 0
head(hy_wide)
hy_wide$GW <- gsub("1","GW",hy_wide$GW)
hy_wide$GL <- gsub("1","GL",hy_wide$GL)
hy_wide$GLGW <- gsub("1","GWGL",hy_wide$GLGW)
hy_wide$TGW <- gsub("1","TGW",hy_wide$TGW)
hy_wide$triat <- paste(hy_wide$GW,hy_wide$GL,hy_wide$GLGW,hy_wide$TGW,sep=",")
hy_wide$triat <- gsub("0,","",hy_wide$triat)
hy_wide$triat <- gsub(",0","",hy_wide$triat)
write.csv(hy_wide[,c(1,6)],"haplatype_genev1.1.csv",row.names = F)
##计算单倍型显著基因箱式图
#在excel中筛选出单倍型显著的source
setwd("D:/胚乳/大修/单倍型")
net <- read.csv("单倍型筛选的网络.csv")
source <- unique(net[c(1,1)])
#留下pvalue最大的snp
tt_total <- data.frame()
for (i in name) {
  tt <- read.table(paste0(i,"_result.txt"),sep="\t")
  tt <- tt[order(tt$V2,tt$V4),]
  topsig <- tt[!duplicated(tt[,c(1,2)]),]
  tt_total <- rbind(tt_total,tt)
}
head(tt_total)
tt_total <- tt_total[!duplicated(tt_total[,c(1,2),]),]
#挑选出source
colnames(source) <- c("V2","source")
tt_total_source <- merge(tt_total,source,by="V2",all=FALSE)
head(tt_total_source)
write.table(tt_total_source,"source_haplotype.txt",row.names=F,col.names = F,quote = F,sep="\t")
###调取筛选SNP的表型数据
cd /public/home/chaohe/prgwas
mkdir source_haplotype
cut -f 3 source_haplotype.txt | uniq | while read i;
do
echo "$i" > "$i"_SNP.txt
plink --bfile HC_genotype --extract "$i"_SNP.txt  --recode --out "$i"_topsig_SNP
paste -d'\t' "$i"_topsig_SNP.ped gl.txt >source_haplotype/"$i"_topsig_SNP_gl.txt
paste -d'\t' "$i"_topsig_SNP.ped gw.txt >source_haplotype/"$i"_topsig_SNP_gw.txt
paste -d'\t' "$i"_topsig_SNP.ped glgw.txt >source_haplotype/"$i"_topsig_SNP_glgw.txt
paste -d'\t' "$i"_topsig_SNP.ped tgw.txt >source_haplotype/"$i"_topsig_SNP_tgw.txt
done
perl -p -i -e 's/ /\t/g' *_topsig_SNP_*.txt

#####补充调取所有source的单倍型
#在excel中筛选出单倍型显著的source
setwd("D:/胚乳/大修/单倍型")
net <- read.csv("")
source <- unique(net[c(1,1)])
#留下pvalue最大的snp
tt_total <- data.frame()
for (i in name) {
  tt <- read.table(paste0(i,"_result.txt"),sep="\t")
  tt <- tt[order(tt$V2,tt$V4),]
  topsig <- tt[!duplicated(tt[,c(1,2)]),]
  tt_total <- rbind(tt_total,tt)
}
head(tt_total)
tt_total <- tt_total[!duplicated(tt_total[,c(1,2),]),]
#挑选出source
colnames(source) <- c("V2","source")
tt_total_source <- merge(tt_total,source,by="V2",all=FALSE)
head(tt_total_source)
write.table(tt_total_source,"source_haplotype.txt",row.names=F,col.names = F,quote = F,sep="\t")
###调取筛选SNP的表型数据
cd /public/home/chaohe/prgwas
mkdir source_haplotype
cut -f 3 source_haplotype.txt | uniq | while read i;
do
echo "$i" > "$i"_SNP.txt
plink --bfile HC_genotype --extract "$i"_SNP.txt  --recode --out "$i"_topsig_SNP
paste -d'\t' "$i"_topsig_SNP.ped gl.txt >source_haplotype/"$i"_topsig_SNP_gl.txt
paste -d'\t' "$i"_topsig_SNP.ped gw.txt >source_haplotype/"$i"_topsig_SNP_gw.txt
paste -d'\t' "$i"_topsig_SNP.ped glgw.txt >source_haplotype/"$i"_topsig_SNP_glgw.txt
paste -d'\t' "$i"_topsig_SNP.ped tgw.txt >source_haplotype/"$i"_topsig_SNP_tgw.txt
done
perl -p -i -e 's/ /\t/g' *_topsig_SNP_*.txt


####source haplotype visualization
#TraesCS3B02G400200
setwd("D:/胚乳/大修/单倍型/显著性")
setwd("D:/胚乳/大修/单倍型/显著性")
#读入snp列表
perl -p -i -e 's/ /\t/g' *txt
list <- read.table("list.txt")
list <- unique(unlist(list))
head(list)
#将对应的基因找出来
net <- read.table("../source_haplotype.txt")
source <- unique(net[c(1,3)])
colnames(source) <- c("gene","snp")
head(source)
#zat_snp <- merge(topsig,zat,by="gene",all=F)
#zat_snp <- zat_snp[which(zat_snp$SNP != "SNP-00153525"),]
#导入基因表达量
tpm <- read.csv("TPM_means.csv")
##grain length
plot_list <- list()
for (i in list) {
  a1 <- read.table(paste0(i,"_topsig_SNP","_gl.txt"),sep="\t",header=F)
  a1 <- a1[,c(1,7,11)]
  a2 <- source[which(source$snp == i),]
  a1$gene <- a2[,1]
  colnames(a1) <- c("stage","genotype","value","gene")
  a1 <- a1[which(a1$genotype != 0),]
  xm <- unique(a1$genotype)
  a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
  a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
  my_comparisons <- list(c("H1", "H2"))
  p1 <- ggboxplot(a1, x = "genotype", y = "value", 
                  color = "genotype", palette = "jco",
                  add = c("jitter"),
                  #bxp.errorbar=T,width = 0.5,
                  add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
    #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
    #scale_fill_manual(values=brewer.pal(12,"Paired")) +
    #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
    #ggtitle("total_triad") + 
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
    xlab(NULL) + ylab("Grain length") +
    stat_compare_means(comparisons = my_comparisons,
                       #aes(group = Haplotype),
                       #label = "p.signif",
                       method = "t.test") 
  plot_list[[i]] <- p1
}
##Grain width
plot_list2 <- list()
for (i in list) {
  a1 <- read.table(paste0(i,"_topsig_SNP","_gw.txt"),sep="\t",header=F)
  a1 <- a1[,c(1,7,11)]
  a2 <- source[which(source$snp == i),]
  a1$gene <- a2[,1]
  colnames(a1) <- c("stage","genotype","value","gene")
  a1 <- a1[which(a1$genotype != 0),]
  xm <- unique(a1$genotype)
  a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
  a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
  my_comparisons <- list(c("H1", "H2"))
  p1 <- ggboxplot(a1, x = "genotype", y = "value", 
                  color = "genotype", palette = "jco",
                  add = c("jitter"),
                  #bxp.errorbar=T,width = 0.5,
                  add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
    #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
    #scale_fill_manual(values=brewer.pal(12,"Paired")) +
    #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
    #ggtitle("total_triad") + 
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
    xlab(NULL) + ylab("Grain width") +
    stat_compare_means(comparisons = my_comparisons,
                       #aes(group = Haplotype),
                       #label = "p.signif",
                       method = "t.test") 
  plot_list2[[i]] <- p1
}
##Grain length/ Grain width
plot_list3 <- list()
for (i in list) {
  a1 <- read.table(paste0(i,"_topsig_SNP","_glgw.txt"),sep="\t",header=F)
  a1 <- a1[,c(1,7,11)]
  a2 <- source[which(source$snp == i),]
  a1$gene <- a2[,1]
  colnames(a1) <- c("stage","genotype","value","gene")
  a1 <- a1[which(a1$genotype != 0),]
  xm <- unique(a1$genotype)
  a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
  a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
  my_comparisons <- list(c("H1", "H2"))
  p1 <- ggboxplot(a1, x = "genotype", y = "value", 
                  color = "genotype", palette = "jco",
                  add = c("jitter"),
                  #bxp.errorbar=T,width = 0.5,
                  add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
    #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
    #scale_fill_manual(values=brewer.pal(12,"Paired")) +
    #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
    #ggtitle("total_triad") + 
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
    xlab(NULL) + ylab("Grain length / grain width") +
    stat_compare_means(comparisons = my_comparisons,
                       #aes(group = Haplotype),
                       #label = "p.signif",
                       method = "t.test") 
  plot_list3[[i]] <- p1
}
##TGW
plot_list4 <- list()
for (i in list) {
  a1 <- read.table(paste0(i,"_topsig_SNP","_tgw.txt"),sep="\t",header=F)
  a1 <- a1[,c(1,7,11)]
  a2 <- source[which(source$snp == i),]
  a1$gene <- a2[,1]
  colnames(a1) <- c("stage","genotype","value","gene")
  a1 <- a1[which(a1$genotype != 0),]
  xm <- unique(a1$genotype)
  a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
  a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
  my_comparisons <- list(c("H1", "H2"))
  p1 <- ggboxplot(a1, x = "genotype", y = "value", 
                  color = "genotype", palette = "jco",
                  add = c("jitter"),
                  #bxp.errorbar=T,width = 0.5,
                  add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
    #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
    #scale_fill_manual(values=brewer.pal(12,"Paired")) +
    #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
    #ggtitle("total_triad") + 
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
    xlab(NULL) + ylab("Thousand grain weight") +
    stat_compare_means(comparisons = my_comparisons,
                       #aes(group = Haplotype),
                       #label = "p.signif",
                       method = "t.test") 
  plot_list4[[i]] <- p1
}
####组图
grid.arrange(plot_list[[1]],plot_list2[[1]],plot_list3[[1]],plot_list4[[1]],
             plot_list[[2]],plot_list2[[2]],plot_list3[[2]],plot_list4[[2]],
             plot_list[[3]],plot_list2[[3]],plot_list3[[3]],plot_list4[[3]],
             plot_list[[4]],plot_list2[[4]],plot_list3[[4]],plot_list4[[4]],
             plot_list[[5]],plot_list2[[5]],plot_list3[[5]],plot_list4[[5]],
             plot_list[[6]],plot_list2[[6]],plot_list3[[6]],plot_list4[[6]],
             plot_list[[7]],plot_list2[[7]],plot_list3[[7]],plot_list4[[7]],
             plot_list[[8]],plot_list2[[8]],plot_list3[[8]],plot_list4[[8]],
             plot_list[[9]],plot_list2[[9]],plot_list3[[9]],plot_list4[[9]],
             plot_list[[10]],plot_list2[[10]],plot_list3[[10]],plot_list4[[10]],
             plot_list[[11]],plot_list2[[11]],plot_list3[[11]],plot_list4[[11]],
             plot_list[[12]],plot_list2[[12]],plot_list3[[12]],plot_list4[[12]],
             plot_list[[13]],plot_list2[[13]],plot_list3[[13]],plot_list4[[13]],
             plot_list[[14]],plot_list2[[14]],plot_list3[[14]],plot_list4[[14]],
             plot_list[[15]],plot_list2[[15]],plot_list3[[15]],plot_list4[[15]],
             plot_list[[16]],plot_list2[[16]],plot_list3[[16]],plot_list4[[16]],
             plot_list[[17]],plot_list2[[17]],plot_list3[[17]],plot_list4[[17]],
             plot_list[[18]],plot_list2[[18]],plot_list3[[18]],plot_list4[[18]],
             plot_list[[19]],plot_list2[[19]],plot_list3[[19]],plot_list4[[19]],
             plot_list[[20]],plot_list2[[20]],plot_list3[[20]],plot_list4[[20]],
             plot_list[[21]],plot_list2[[21]],plot_list3[[21]],plot_list4[[21]],
             plot_list[[22]],plot_list2[[22]],plot_list3[[22]],plot_list4[[22]],
             plot_list[[23]],plot_list2[[23]],plot_list3[[23]],plot_list4[[23]],
             plot_list[[24]],plot_list2[[24]],plot_list3[[24]],plot_list4[[24]],
             plot_list[[25]],plot_list2[[25]],plot_list3[[25]],plot_list4[[25]],
             plot_list[[26]],plot_list2[[26]],plot_list3[[26]],plot_list4[[26]],
             plot_list[[27]],plot_list2[[27]],plot_list3[[27]],plot_list4[[27]],
             plot_list[[28]],plot_list2[[28]],plot_list3[[28]],plot_list4[[28]],
             plot_list[[29]],plot_list2[[29]],plot_list3[[29]],plot_list4[[29]],
             plot_list[[30]],plot_list2[[30]],plot_list3[[30]],plot_list4[[30]],
             plot_list[[31]],plot_list2[[31]],plot_list3[[31]],plot_list4[[31]],
             plot_list[[32]],plot_list2[[32]],plot_list3[[32]],plot_list4[[32]],
             plot_list[[33]],plot_list2[[33]],plot_list3[[33]],plot_list4[[33]],
             plot_list[[34]],plot_list2[[34]],plot_list3[[34]],plot_list4[[34]],
             plot_list[[35]],plot_list2[[35]],plot_list3[[35]],plot_list4[[35]],
             plot_list[[36]],plot_list2[[36]],plot_list3[[36]],plot_list4[[36]],
             plot_list[[37]],plot_list2[[37]],plot_list3[[37]],plot_list4[[37]],
             plot_list[[38]],plot_list2[[38]],plot_list3[[38]],plot_list4[[38]],
             plot_list[[39]],plot_list2[[39]],plot_list3[[39]],plot_list4[[39]],
             plot_list[[40]],plot_list2[[40]],plot_list3[[40]],plot_list4[[40]],
             plot_list[[41]],plot_list2[[41]],plot_list3[[41]],plot_list4[[41]],
             plot_list[[42]],plot_list2[[42]],plot_list3[[42]],plot_list4[[42]],
             plot_list[[43]],plot_list2[[43]],plot_list3[[43]],plot_list4[[43]],
             plot_list[[44]],plot_list2[[44]],plot_list3[[44]],plot_list4[[44]],
             plot_list[[45]],plot_list2[[45]],plot_list3[[45]],plot_list4[[45]],
             plot_list[[46]],plot_list2[[46]],plot_list3[[46]],plot_list4[[46]],
             plot_list[[47]],plot_list2[[47]],plot_list3[[47]],plot_list4[[47]],
             plot_list[[48]],plot_list2[[48]],plot_list3[[48]],plot_list4[[48]],
             plot_list[[49]],plot_list2[[49]],plot_list3[[49]],plot_list4[[49]],
             plot_list[[50]],plot_list2[[50]],plot_list3[[50]],plot_list4[[50]],
             plot_list[[51]],plot_list2[[51]],plot_list3[[51]],plot_list4[[51]],
             plot_list[[52]],plot_list2[[52]],plot_list3[[52]],plot_list4[[52]],
             plot_list[[53]],plot_list2[[53]],plot_list3[[53]],plot_list4[[53]],
             plot_list[[54]],plot_list2[[54]],plot_list3[[54]],plot_list4[[54]],
             plot_list[[55]],plot_list2[[55]],plot_list3[[55]],plot_list4[[55]],
             plot_list[[56]],plot_list2[[56]],plot_list3[[56]],plot_list4[[56]],
             plot_list[[57]],plot_list2[[57]],plot_list3[[57]],plot_list4[[57]],
             plot_list[[58]],plot_list2[[58]],plot_list3[[58]],plot_list4[[58]],
             plot_list[[59]],plot_list2[[59]],plot_list3[[59]],plot_list4[[59]],
             plot_list[[60]],plot_list2[[60]],plot_list3[[60]],plot_list4[[60]],
             plot_list[[61]],plot_list2[[61]],plot_list3[[61]],plot_list4[[61]],
             plot_list[[62]],plot_list2[[62]],plot_list3[[62]],plot_list4[[62]],
             plot_list[[63]],plot_list2[[63]],plot_list3[[63]],plot_list4[[63]],
             plot_list[[64]],plot_list2[[64]],plot_list3[[64]],plot_list4[[64]],
             plot_list[[65]],plot_list2[[65]],plot_list3[[65]],plot_list4[[65]],
             plot_list[[66]],plot_list2[[66]],plot_list3[[66]],plot_list4[[66]],
             plot_list[[67]],plot_list2[[67]],plot_list3[[67]],plot_list4[[67]],
             plot_list[[68]],plot_list2[[68]],plot_list3[[68]],plot_list4[[68]],
             plot_list[[69]],plot_list2[[69]],plot_list3[[69]],plot_list4[[69]],
             plot_list[[70]],plot_list2[[70]],plot_list3[[70]],plot_list4[[70]],
             plot_list[[71]],plot_list2[[71]],plot_list3[[71]],plot_list4[[71]],
             plot_list[[72]],plot_list2[[72]],plot_list3[[72]],plot_list4[[72]],
             plot_list[[73]],plot_list2[[73]],plot_list3[[73]],plot_list4[[73]],
             plot_list[[74]],plot_list2[[74]],plot_list3[[74]],plot_list4[[74]],
             plot_list[[75]],plot_list2[[75]],plot_list3[[75]],plot_list4[[75]],
             plot_list[[76]],plot_list2[[76]],plot_list3[[76]],plot_list4[[76]],
             plot_list[[77]],plot_list2[[77]],plot_list3[[77]],plot_list4[[77]],
             plot_list[[78]],plot_list2[[78]],plot_list3[[78]],plot_list4[[78]],
             plot_list[[79]],plot_list2[[79]],plot_list3[[79]],plot_list4[[79]],
             plot_list[[80]],plot_list2[[80]],plot_list3[[80]],plot_list4[[80]],
             nrow=40,ncol=8)     %>%  ggsave("单倍型分析.pdf",.,width=300,height=5000, units="mm", limitsize = FALSE)
###单独绘制ERF5的单倍型结果
i = "IND-155856731"
j = "IND-155856781"
a1 <- read.table(paste0(i,"_topsig_SNP","_gl.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p1 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Grain length") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")
a1 <- read.table(paste0(i,"_topsig_SNP","_gw.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p2 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Grain width") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")

a1 <- read.table(paste0(i,"_topsig_SNP","_glgw.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p3 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Grain length / grain width") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")
a1 <- read.table(paste0(i,"_topsig_SNP","_tgw.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p4 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Thousand grain weight") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")
i = "IND-155856781"
a1 <- read.table(paste0(i,"_topsig_SNP","_gl.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p5 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Grain length") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")
a1 <- read.table(paste0(i,"_topsig_SNP","_gw.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p6 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Grain width") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")

a1 <- read.table(paste0(i,"_topsig_SNP","_glgw.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p7 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Grain length / grain width") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")
a1 <- read.table(paste0(i,"_topsig_SNP","_tgw.txt"),sep="\t",header=F)
a1 <- a1[,c(1,7,11)]
a2 <- source[which(source$snp == i),]
a1$gene <- a2[,1]
colnames(a1) <- c("stage","genotype","value","gene")
a1 <- a1[which(a1$genotype != 0),]
xm <- unique(a1$genotype)
a1$genotype <- gsub(xm[[1]],"H1",a1$genotype)
a1$genotype <- gsub(xm[[2]],"H2",a1$genotype)
my_comparisons <- list(c("H1", "H2"))
p8 <- ggboxplot(a1, x = "genotype", y = "value",
                color = "genotype", palette = "jco",
                add = c("jitter"),
                #bxp.errorbar=T,width = 0.5,
                add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggtitle(unique(paste(unique(a4$gene.y),a2[,2],sep="_"))) +
  xlab(NULL) + ylab("Thousand grain weight") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")
grid.arrange(p5, p6, p8,
             nrow=1,ncol=3)     %>%  ggsave("ERF5_haplotype.pdf",.,width=210,height=60, units="mm")

######20240630胚乳相关基因的组蛋白修饰情况
####淀粉合成途径所有相关基因的表达水平变化和组蛋白修饰变化
#表达水平
setwd("D:/胚乳/大修/组蛋白与表达量相关性/淀粉相关")
star <- read.csv("最终使用淀粉合成相关基因20240630.csv")
#导入表达量和cluster
tpm <- read.csv("D:/胚乳/时序变化基因/dynamic_gene_expression_filter.csv")
head(tpm)
colnames(star) <- c("X","cluster")
star_tpm <- merge(star,tpm,by="X",all.x=FALSE)
head(star_tpm)
star <- star_tpm
colnames(star) <- c("id","symbol","4DAF","7DAF","14DAF","18DAF","cluster")
head(star)
row.names(star) <- star[,1]
#annotation_row = data.frame(GeneClass = factor(rep(c("SuSy", "RuBisCO", "PGK","FBA", "AGPase", "BT1", "Waxy","SSS", "SBE", "DBE","GBSSII"), c(4,34,18,18,11,1,3,5,1,1,6))))
#rownames(annotation_row) = star$id
test <- log2(star_tpm[,3:6])
#test$max <- rowMax(as.matrix(test))
#test <- test[,-5]
#pheatmap(test,
#         #annotation_row = annotation_row, 
#         cluster_row = F, show_rownames = T, 
#         scale="row",
#         cluster_col = FALSE,color = c("#A8C0FF","red"))
####组蛋白修饰水平
setwd("D:/胚乳/ChIP/bed/定量/时序比较定量/")
#H3K9ac
name1 <- c("4DPA_H3K9ac","7DPA_H3K9ac","14DPA_H3K9ac","18DPA_H3K9ac",
           "4DPA_H3K4me3","7DPA_H3K4me3","14DPA_H3K4me3","18DPA_H3K4me3")
fram <- data.frame()
for (i in name1) {
  dt <- read.csv(paste0("0902_",i,"_merge_gene_peak_tpms.csv"))
  dt <- dt[which(dt$Tyep == "promoter" | dt$Tyep == "genebody"),]
  dt <- aggregate(dt$TPM,list(dt$Gene),sum)
  dt$type <- rep(i,nrow(dt))
  fram <- rbind(fram,dt)
}
name2 <- c("4DPA_H3K27me3","7DPA_H3K27me3","14DPA_H3K27me3","18DPA_H3K27me3")
for (i in name2) {
  dt <- read.csv(paste0("0902_",i,"_merge_gene_peak_tpms.csv"))
  dt <- dt[which(dt$Tyep == "promoter" | dt$Tyep == "genebody"),]
  dt <- aggregate(dt$TPM,list(dt$Gene),sum)
  dt$type <- rep(i,nrow(dt))
  fram <- rbind(fram,dt)
}
colnames(fram) <- c("id","TPM","type")
##留下淀粉合成相关基因
fram_exp <- merge(star, fram,  by="id",all.x=F)
#长变宽
data_new = spread(fram_exp,key=type,value=TPM)
#data_new <- data_new[,c(1,2,3,4,5,6,15,18,9,12,16,19,10,13,17,20,11,14)]
data_new <- data_new[,c(1:6,14,17,8,11,15,18,9,12,16,19,10,13)]
data_new[is.na(data_new)] <- 0
head(data_new)
#write.csv(data_new,"total_starch_gene.csv")  #excel中挑顺序，将NA变成0
#write.csv(data_new,"D:/胚乳/大修/组蛋白与表达量相关性/淀粉相关/total_starch_geneR1.csv")   ##EXCEL筛选组蛋白修饰与表达量趋势相一致的
#重新读入
#data_new <- read.csv("total_starch_gene.csv",row.names=1)
#删除全是NA的行或列的函数
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
row.names(data_new) <- data_new$id
head(data_new)
name <- unique(data_new$symbol)
plot_list <- list()
for (i in name) {
  mm <- data_new[which(data_new$symbol == i),]
  RNA <- mm[,c(3:6)]
  k273 <- mm[,c(7:10)]
  k43 <- mm[,c(11:14)]
  k9ac <- mm[,c(15:18)]
  RNA[RNA == 0] <- NA
  k273[k273 == 0] <- NA
  k43[k43 == 0] <- NA
  k9ac[k9ac == 0] <- NA
  #RNA <- removeRowsAllNa(RNA)
  #k273 <- removeRowsAllNa(k273)
  #k43 <- removeRowsAllNa(k43)
  #k9ac <- removeRowsAllNa(k9ac)
  colnames(k273) <- colnames(RNA)
  colnames(k43) <- colnames(RNA)
  colnames(k9ac) <- colnames(RNA)
  annotation_row = data.frame(GeneClass = factor(rep(c("RNA", "K273", "K43","K9ac"), c(nrow(RNA),nrow(k273),nrow(k43),nrow(k9ac)))))
  #Group = factor(rep(c("RNA-seq","H3K27me3","H3K4me3","H3K9ac"),times = c(4,4,4,4)))#分组信息，用于热图分割
  #Group = factor(Group,levels = c("RNA-seq","H3K27me3","H3K4me3","H3K9ac"))
  #tt <- mm[,3:18]
  p1 <- Heatmap(t(scale(t(mm[,3:6]))),#表达矩阵
                #col = colorRampPalette(c("#A8C0FF","white","red"))(100),#颜色定义
                show_row_names = F,#不展示行名
                #top_annotation = top_annotation,#顶部分组信息
                #column_split = Group,#用group信息将热土分开，以group聚类
                cluster_columns = FALSE,
                column_title = i ,#不显示列标题
                show_column_names = T)#不显示列名
  p2 <- Heatmap(t(scale(t(mm[,7:10]))),#表达矩阵
                #col = colorRampPalette(c("#A8C0FF","white","red"))(100),#颜色定义
                show_row_names = F,#不展示行名
                #top_annotation = top_annotation,#顶部分组信息
                #column_split = Group,#用group信息将热土分开，以group聚类
                cluster_columns = FALSE,
                column_title = i ,#不显示列标题
                show_column_names = T)#不显示列名
  p3 <- Heatmap(t(scale(t(mm[,11:14]))),#表达矩阵
                #col = colorRampPalette(c("#A8C0FF","white","red"))(100),#颜色定义
                show_row_names = F,#不展示行名
                #top_annotation = top_annotation,#顶部分组信息
                #column_split = Group,#用group信息将热土分开，以group聚类
                cluster_columns = FALSE,
                column_title = i ,#不显示列标题
                show_column_names = T)#不显示列名
  p4 <- Heatmap(t(scale(t(mm[,15:18]))),#表达矩阵
                #col = colorRampPalette(c("#A8C0FF","white","red"))(100),#颜色定义
                show_row_names = F,#不展示行名
                #top_annotation = top_annotation,#顶部分组信息
                #column_split = Group,#用group信息将热土分开，以group聚类
                cluster_columns = FALSE,
                column_title = i ,#不显示列标题
                show_column_names = T)#不显示列名
  
  h1 <- draw(p1 + p2 + p3 + p4, newpage = F, auto_adjust = FALSE, column_title = i, column_title_gp = gpar(fontsize = 15, fontface = "bold"), heatmap_legend_side = "right")
  plot_list[[i]] <- h1
  pdf(paste0("D:/胚乳/大修/组蛋白与表达量相关性/淀粉相关/",i,".pdf"))
  draw(h1)
  dev.off()
}

for (i in name) {
  pdf(paste0("D:/胚乳/大修/组蛋白与表达量相关性/淀粉相关/",i,".pdf"))
  plot_list[[i]]
  dev.off()
}

#####统计G1P前和G1P后的基因的数量
#分cluster
setwd("D:/胚乳/大修/组蛋白与表达量相关性/淀粉相关")
star <- read.csv("最终使用淀粉合成相关基因20240630.csv")
#导入表达量和cluster
tpm <- read.csv("D:/胚乳/时序变化基因/dynamic_gene_expression_filter.csv")
head(tpm)
colnames(star) <- c("X","cluster")
star_tpm <- merge(star,tpm,by="X",all.x=FALSE)
head(star_tpm)
star <- star_tpm
colnames(star) <- c("id","symbol","4DAF","7DAF","14DAF","18DAF","cluster")
head(star)
#G1P之前的酶
be <- star[which(star$symbol == "SuSy" | star$symbol == "RuBisCO" | star$symbol == "FBA" |
                   star$symbol == "PGK"),]
af <-  star[which(star$symbol == "AGPase" | star$symbol == "BT1" | star$symbol == "Waxy" |
                    star$symbol == "DBE" | star$symbol == "SBE" | star$symbol == "SS"),]
be_number <- aggregate(be$symbol,list(be$cluster),length)
af_number <- aggregate(af$symbol,list(af$cluster),length)
##绘制立体3D饼图
be_number$Freq <- round(be_number$x/sum(be_number$x)*100,2)
af_number$Freq <- round(af_number$x/sum(af_number$x)*100,2)
#绘图
#install.packages("plotrix")
library("plotrix")
y=be_number
y$Group.1 <- paste0("Cluster",y$Group.1)
percent <- paste(y$Freq,"%",sep="")  #add % to labels
lbls <- paste(y$Group.1,y$x,sep="\n") #换行
lbls <- paste0(lbls," (",percent,")") #换行

pdf("D:/胚乳/大修/组蛋白与表达量相关性/淀粉相关/G1P_before.pdf",family="ArialMT")
pie3D(y$Freq,radius=0.8,height=0.1,labels=lbls,main="Before G1P was generated",
      col = c("Cluster1"="#eb1e2c","Cluster2"="#fd6f30","Cluster3"="#f9a729","Cluster4"="#f9d23c"))
dev.off()

y=af_number
y$Group.1 <- paste0("Cluster",y$Group.1)
percent <- paste(y$Freq,"%",sep="")  #add % to labels
lbls <- paste(y$Group.1,y$x,sep="\n") #换行
lbls <- paste0(lbls," (",percent,")") #换行
pdf("D:/胚乳/大修/组蛋白与表达量相关性/淀粉相关/G1P_after.pdf",family="ArialMT")
pie3D(y$Freq,radius=0.8,height=0.1,labels=lbls,main="Before G1P was generated",
      col = c("Cluster1"="#eb1e2c","Cluster2"="#fd6f30","Cluster3"="#f9a729","Cluster4"="#f9d23c"))
dev.off()

####DAP-seq验证网络的调控关系
setwd("/public/home/chaohe/prdap")
dap <- read.csv("/public/home/chaohe/AG_chip/dap_peak_annotation.csv")
link <- read.table("dap_tf_link.txt")
dap$V4 <- gsub("GSE[0-9]*_","",dap$V4)
dap$V4 <- gsub("GSM[0-9]*_","",dap$V4)
dap$V4 <- gsub("_PE_peaks.filter.p10.bed","",dap$V4)
dap$V4 <- gsub(".final.bed","",dap$V4)
dap$V4 <- gsub("macs_","",dap$V4)
dap$V4 <- gsub("_PE_peaks.bed","",dap$V4)
dap$V4 <- gsub("_rep2","",dap$V4)  #AP2-DREB-1A-1,AP2-DREB-1B-1,AP2-DREB-1D-1,AP2_1_A,AP2-ERF-6A-1两个重复
dap$V4 <- gsub("_DAP-seq_Rep2_peaks.bed","",dap$V4)  #AP2-DREB-1A-1,AP2-DREB-1B-1,AP2-DREB-1D-1,AP2_1_A,AP2-ERF-6A-1两个重复
dap <- dap[which(dap$V4 != "Halo-merge_PE_peaks.bed" & dap$V4 != "control.pks.bed"),]
dap <- unique(dap)
name <- read.table("DAP-seq_id.txt",sep="\t",header=TRUE)
dap$symbol <- dap$V4
dap_name <- merge(dap,name,by="symbol",all=F)
#提取dap中的source
write.csv(dap_name,"dap_name.csv",row.names=F)
#linux
awk -F[,] '{print $15"\t"$1"\t"$10}' dap_name.csv >dap_name.bed
perl -p -i -e 's/\"//g' dap_name.bed
awk 'NR==FNR {ids[$1]; next} $1 in ids' DAP-seq_id.txt dap_name.bed > matched_output.bed
##绘制overlap 比例
setwd("D:/胚乳/TF印记网络/DAP")
dap <- read.table("matched_output.bed",sep="\t")
link <- read.table("dap_tf_link.txt",header=T,sep="\t")
link <- link[,c(1,3)]
head(dap)
head(link)
colnames(dap) <- c("source","dap_source_symbol","target")
#name <- unique(link$source)
name <- c("TraesCS3A02G432500","TraesCS3B02G468400","TraesCS3D02G425800","TraesCS5B02G075300",
          "TraesCS7B02G142200", "TraesCS6D02G225700")
aabb_num_total <- data.frame()
aabb_total <- data.frame()
for (i in name) {
  aa <- unique(dap[which(dap$source == i),])
  bb <- unique(link[which(link$source == i),])
  aabb <- merge(bb,aa,by="target",all.x=TRUE)
  aabb$type <- aabb$source.y
  aabb[is.na(aabb)] <- "Other"
  xx <- aabb[which(aabb$dap_source_symbol != "Other"),4]
  aabb$dap_source_symbol <- unique(xx)
  aabb$type <- gsub(i,"DAP_supported",aabb$type)
  aabb_num <- aggregate(aabb$target,list(aabb$source.x,aabb$dap_source_symbol,aabb$type),length)
  aabb_num_total <- rbind(aabb_num_total,aabb_num)
  aabb_total <- rbind(aabb_total,aabb)
}
head(aabb_num_total)
write.csv(aabb_total,"DAP_link.csv",row.names = F)
library(VennDiagram)
#绘制柱形图
library(ggpubr)
pdf("DAP-link_overlap.pdf",family="ArialMT")
ggbarplot(aabb_num_total, "Group.2", "x",
          fill = "Group.3", 
          label = TRUE, lab.pos = "in", lab.col = "white",
          palette = c("#FF7F00", "#4DAF4A"),
          position = position_stack()) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,5200)) +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()
#计算比例
library(dplyr)
proportions <- aabb_num_total %>%
  group_by(Group.1, Group.2) %>%
  mutate(Proportion = x / sum(x)) %>%
  ungroup()
head(proportions)
write.csv(proportions,"DAP_link_overlap.csv",row.names=F)

####调取ERF5 DAP-seq与网络的交集，可视化并做GO
setwd("D:/胚乳/TF印记网络/DAP")
dap <- read.csv("DAP_link.csv")
dap <- dap[which(dap$Source_symbol == "AP2-ERF-6D-1"),]
sou <- dap[which(dap$Type == "DAP_supported"),]
head(sou)
write.csv(sou,"ERF5_DAP.csv",row.names = F)

###TGT GO富集分析
library(data.table)
library(tidyr)
setwd("D:/胚乳/TF印记网络/DAP")
go <- read.csv("TGT_GOTable20240808215944.csv",row.names = 1)
#go <- go[which(go$Sample == "Balance_unblance_unblance" | go$Sample == "Balance_unblance_Balance" | go$Sample == "unalance_unblance_Balance" | go$Sample == "Balance_unblance_Balance") ,]
go <- go[,c(1,2,3,5,7)]
go <- separate(data = go, col = Ratio.in.foreground.list, into = c("value1", "value2"), sep = "/")
go$value <- as.numeric(go$value1)/as.numeric(go$value2)
go$value1 <- as.numeric(go$value1)
go$value2 <- as.numeric(go$value2)
go$Description <- factor(go$Description,levels=unique(go$Description),ordered=TRUE)
head(go)
top<- go %>% group_by(Group) %>% top_n(n = 10, wt = -log10(FDR))
#write.csv(top,"top.csv")
#top <- read.csv("top.csv")
top <- top[order(top$Group),]  ##跑两次这个和下面两个，点的排队趋势就对了
top$Description <- factor(top$Description,levels=unique(top$Description),ordered=TRUE)
#G0 气泡图.

pdf("GO_ERF5_DAP.pdf",width=10,family="ArialMT")
ggplot(top,aes(x=Group,y=Description))+
  geom_point(aes(size=`value1`,
                 color=`FDR`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  scale_color_gradient(low="#61C0BA",high="#EEF8B7") +
  labs(x=NULL,y=NULL)
dev.off()
###DAP和网络的venn图
library(VennDiagram)
setwd("D:/胚乳/TF印记网络/DAP")
dap <- read.table("matched_output.bed",sep="\t")
link <- read.table("dap_tf_link.txt",header=T,sep="\t")
link <- link[,c(1,3)]
head(link)
venn_list <- list(Network = unique(link[which(link$source == "TraesCS6D02G225700"),2]), DAP = unique(dap[which(dap$V2 == "AP2-ERF-6D-1"),3]))
venn.plot<-venn.diagram(venn_list, filename = NULL, 
                        fill = c('#7FC97F', '#BEAED4'), alpha = 0.50, 
                        cat.col = c('#7FC97F', '#BEAED4'), cat.cex = 1.5, cat.fontfamily = 'serif',
                        col = c('#7FC97F', '#BEAED4'), cex = 1.5, fontfamily = 'serif')
pdf(file="dap_link.pdf",family="ArialMT")
grid.draw(venn.plot)
dev.off()

####桑基图显示不同时期不同类型基因不对称 表达
#准备数据
setwd("D:/胚乳/大修/H3K9ac与ATAC-seq的peak的overlap")
tridat01 <- read.csv("DAP0_RNA_tridat.csv",header=T,row.names=1)
tridat02 <- read.csv("DAP2_RNA_tridat.csv",header=T,row.names=1)
tridat03 <- read.csv("DAP4_RNA_tridat.csv",header=T,row.names=1)
tridat04 <- read.csv("DAP6_RNA_tridat.csv",header=T,row.names=1)
tridat05 <- read.csv("DAP8_RNA_tridat.csv",header=T,row.names=1)
tridat06 <- read.csv("DAP12_RNA_tridat.csv",header=T,row.names=1)
tridat07 <- read.csv("DAP16_RNA_tridat.csv",header=T,row.names=1)
tridat08 <- read.csv("DAP22_RNA_tridat.csv",header=T,row.names=1)
#row.names(tridat01) <- tridat01$A_id
#row.names(tridat02) <- tridat02$A_id
#row.names(tridat03) <- tridat03$A_id
#row.names(tridat04) <- tridat04$A_id
head(tridat01)
head(tridat02)
head(tridat03)
head(tridat04)
total1<-merge(tridat01,tridat02,by="row.names",all=T)
head(total1)
rownames(total1) <- total1[,1]
total2<-merge(total1,tridat03,by="row.names",all=T)
rownames(total2) <- total2[,1]
total3<-merge(total2,tridat04,by="row.names",all=T)
#total3 <- total3[,c(1,7,15,23,31)]
total3 <- total3[,c(1,7,12,17,22)]
head(total4)
rownames(total3) <- total3[,1]
total4<-merge(total3,tridat05,by="row.names",all=T)
#total3 <- total3[,c(1,7,15,23,31)]
total4 <- total4[,c(1,3,4,5,6,10)]
rownames(total4) <- total4[,1]
total5<-merge(total4,tridat06,by="row.names",all=T)
#total3 <- total3[,c(1,7,15,23,31)]
total5 <- total5[,c(1,3:7,11)]
rownames(total5) <- total5[,1]
total6<-merge(total5,tridat07,by="row.names",all=T)
#total3 <- total3[,c(1,7,15,23,31)]
total6 <- total6[,c(1,3:8,12)]
rownames(total6) <- total6[,1]
total7<-merge(total6,tridat08,by="row.names",all=T)
total7 <- total7[,c(1,3:9,13)]
write.csv(total7,"gene_triad_variation.csv")
total7$cluster <- paste(total7$group.x,total7$group.y,total7$group.x.1,total7$group.y.1,
                        total7$group.x.2,total7$group.y.2,total7$group.x.3,total7$group.y.3,sep="_")
#统计数量
total7_number <- aggregate(total7$Row.names,list(total7$group.x,total7$group.y,total7$group.x.1,total7$group.y.1,
                                                 total7$group.x.2,total7$group.y.2,total7$group.x.3,total7$group.y.3),FUN=length)
head(total7_number)
write.csv(total7_number,"sankeyR1.csv")
#excel操作后导入数据
#install.packages("networkD3")
#install.packages("d3Network")
library(networkD3)
library("d3Network")
#san <- read.csv("sankey_excel.csv",header=T)
rna_12 <- aggregate(total7$group.x,list(total7$group.x,total7$group.y ),length)
rna_23 <- aggregate(total7$group.y,list(total7$group.y,total7$group.x.1 ),length)
rna_34 <- aggregate(total7$group.x.1,list(total7$group.x.1,total7$group.y.1 ),length)
rna_45 <- aggregate(total7$group.y.1,list(total7$group.y.1,total7$group.x.2 ),length)
rna_56 <- aggregate(total7$group.x.2,list(total7$group.x.2,total7$group.y.2 ),length)
rna_67 <- aggregate(total7$group.x.2,list(total7$group.y.2,total7$group.x.3 ),length)
rna_78 <- aggregate(total7$group.x.2,list(total7$group.x.3,total7$group.y.3 ),length)
rna_12$Group.1 <-  paste("DAP0",rna_12$Group.1,sep="_")
rna_12$Group.2 <- paste("DAP2",rna_12$Group.2,sep="_")
rna_23$Group.1 <-  paste("DAP2",rna_23$Group.1,sep="_")
rna_23$Group.2 <- paste("DAP4",rna_23$Group.2,sep="_")
rna_34$Group.1 <-  paste("DAP4",rna_34$Group.1,sep="_")
rna_34$Group.2 <- paste("DAP6",rna_34$Group.2,sep="_")
rna_45$Group.1 <-  paste("DAP6",rna_45$Group.1,sep="_")
rna_45$Group.2 <- paste("DAP8",rna_45$Group.2,sep="_")
rna_56$Group.1 <-  paste("DAP8",rna_56$Group.1,sep="_")
rna_56$Group.2 <- paste("DAP12",rna_56$Group.2,sep="_")
rna_67$Group.1 <-  paste("DAP12",rna_67$Group.1,sep="_")
rna_67$Group.2 <- paste("DAP16",rna_67$Group.2,sep="_")
rna_78$Group.1 <-  paste("DAP16",rna_78$Group.1,sep="_")
rna_78$Group.2 <- paste("DAP22",rna_78$Group.2,sep="_")
head(rna_12)
head(rna_23)
head(rna_34)
#桑基图
sankey <- rbind(rna_12,rna_23,rna_34,rna_45,rna_56,rna_67,rna_78)
library(networkD3)
library("d3Network")
san <- sankey
colnames(san) <- c("Source","Target","Value")
head(san)
head(san)
Sankeylinks<-san
Sankeynodes<-data.frame(name=unique(c(Sankeylinks$Source,Sankeylinks$Target)),stringsAsFactors=FALSE)  
Sankeynodes$index<-0:(nrow(Sankeynodes) - 1)
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="Source",by.y="name")
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="Target",by.y="name")
Sankeydata<-Sankeylinks[,c(4,5,3)];names(Sankeydata)<-c("Source","Target","Value")
Sankeyname<-Sankeynodes[,1,drop=FALSE]
##绘图
sankeyNetwork(Links=Sankeydata,Nodes=Sankeyname, Source ="Source",
              Target = "Target", Value = "Value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30,
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"))
#link着色
# Colour links
Sankeydata$energy_type <- sub(' .*', '', Sankeyname[Sankeydata$Source + 1, 'name'])
sankeyNetwork(Links=Sankeydata,Nodes=Sankeyname, Source ="Source",
              Target = "Target", Value = "Value", NodeID = "name",
              LinkGroup="energy_type",
              units = "TWh", fontSize = 12, nodeWidth = 30)

######不同表达模式基因的表达热图
setwd("D:/胚乳/时序变化基因")
spe <- read.csv("单时期特异和前期后期高表达基因列表R1.csv",header=T,row.name=1,check.names = F)
#spe <- spe[order(spe[,5],spe[,1],spe[,2],spe[,3],spe[,4]),]
head(spe)
library(tibble)
spe <- spe %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer( cols =  c("4Days":"18Days"),
                names_to = 'stage',
                values_to = 'expr')
head(spe)
library(RColorBrewer)
spe$cluster <- as.factor(spe$cluster)
spe$stage <- factor(spe$stage,levels=c("4Days","7Days","14Days","18Days"),ordered=T)
spe$logtpm <- log2(spe$expr+1)
# Calculate counts per cluster
cluster_counts <- spe %>%
  group_by(cluster, stage) %>%
  summarise(count = n(), .groups = 'drop')
# Create the plot
pdf("dynamic_time_heatmap.pdf",family="ArialMT",height=6,width = 8)
ggplot(spe, aes(stage, logtpm, group = sample)) +
  geom_line(aes(colour = cluster), size = 0.01) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", size = 1.6, color = "#c51b7d") +
  facet_wrap(cluster ~ ., scales = "free_y") +
  scale_colour_manual(values = brewer.pal(7, "Set2")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 8, face = "bold")
  ) +
  geom_text(
    data = cluster_counts,
    aes(x = Inf, y = Inf, label = paste("Count:", count)),
    hjust = 1, vjust = 1,
    size = 3,
    color = "black",
    inherit.aes = FALSE
  )
dev.off()

###Crispr-表型可视化
setwd("D:/胚乳/大修/雨琦/突变体图片")
cr <- read.csv("cripsr突变体的表型.csv")
head(cr)
colnames(cr) <- c("name","value","type")
my_comparisons <- list(c("D1-4", "Fielder"))
cr1 <- cr[which(cr$type == "GL"),]
q1 <- ggboxplot(cr1, x = "name", y = "value", 
          color = "name", palette = "jco",
          add = c("jitter"),
          #bxp.errorbar=T,width = 0.5,
          add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") + 
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab(NULL) + ylab("Grain length") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test") 
cr2 <- cr[which(cr$type == "GW"),]
q2 <- ggboxplot(cr2, x = "name", y = "value", 
          color = "name", palette = "jco",
          add = c("jitter"),
          #bxp.errorbar=T,width = 0.5,
          add.params = list(size = 3,fill="Haplotype",alpha= 0.2)) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#9467bd")) +
  #scale_fill_manual(values=brewer.pal(12,"Paired")) +
  #scale_fill_manual(values=c("#E24E59","#E6A15B","#67C3CD","#9467bd")) +
  #ggtitle("total_triad") + 
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  xlab(NULL) + ylab("Grain length") +
  stat_compare_means(comparisons = my_comparisons,
                     #aes(group = Haplotype),
                     #label = "p.signif",
                     method = "t.test")
##组图
grid.arrange(q1,q2,
             nrow=1,ncol=2)     %>%  ggsave("crispr-cas9统计.pdf",.,width=150,height=100, units="mm")                 
