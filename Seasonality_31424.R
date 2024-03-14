library(edgeR)
library(ggplot2)
library(stringr)
library(pcaMethods)
library(mixOmics)

##Load in human counts from server
human_counts_all.txt <- read.delim("/Volumes/projects-t3/SerreDLab-3/kieran.tebben/mali/hisat2_outputs/count_tables/human_counts_all.txt.gz", comment.char="#")

##Read in counts, remove globin genes
globin_genes <- c("HBB", "HBD", "HBA1", "HBA2", "HBE1", "HBEGF", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
human_counts_all.txt <- human_counts_all.txt[!(human_counts_all.txt$Geneid %in% globin_genes),]

##Subset by project
projects <- read.delim("percent_duplicate_human_edited.txt")
sample_list <- projects[str_detect(projects$Project, "seasonality"), ]
sample_list$Sample <- paste("X", sample_list$Sample, sep="")
sample_list$Sample <- paste0(sample_list$Sample, ".subset.human.rmdup.bam")
sample_list <- as.list(sample_list$Sample)

counts <- human_counts_all.txt[,colnames(human_counts_all.txt) %in% sample_list] 

##Clean up colnames to be only sample ID
for ( col in 1:ncol(counts)){
  colnames(counts)[col] <-  sub(".subset.human.rmdup.bam*", "", colnames(counts)[col])
}

for ( col in 1:ncol(counts)){
  colnames(counts)[col] <-  sub("*X", "", colnames(counts)[col])
}

rownames(counts) <- human_counts_all.txt$Geneid

all <- read.csv("mali_samples_allvariables.csv")
all <- all[all$id %in% colnames(counts),]
all <- all[ order(match(all$id, colnames(counts))), ]
all <- all[-54,]

label <- read.delim("seasonality_variables.txt")
label <- label[ order(match(label$full_id, colnames(counts))), ]
label$season_time <- sub(" ", "_", label$season_time)

### Paired samples early first, late second ###
label_pairedEvL <- subset(label, include_paired_early_vs_late == "Y")
counts_pairedEvL <- counts[,colnames(counts) %in% label_pairedEvL$full_id]

dgList_pairedEvL <- DGEList(counts=counts_pairedEvL, genes=human_counts_all.txt$Geneid)

#Filter out lowly expressed genes 
keep <- rowSums(cpm(dgList_pairedEvL)>10) >= (ncol(counts_pairedEvL)/2)
table(keep)
dgList_pairedEvL_filtered <- dgList_pairedEvL[keep, , keep.lib.sizes=FALSE]

dim(dgList_pairedEvL_filtered)

dgList_pairedEvL_filtered <- calcNormFactors(dgList_pairedEvL_filtered, method="TMM")

label_pairedEvL <- label_pairedEvL[order(match(label_pairedEvL$full_id, colnames(dgList_pairedEvL_filtered$counts))),]
dgList_pairedEvL_filtered$samples$group <- label_pairedEvL$season_time

all_EvL <- all[all$id %in% colnames(dgList_pairedEvL_filtered),]
all_EvL$full_id <- all_EvL$id
allvars_EvL <- merge(all_EvL, label_pairedEvL, by = "full_id")

#Compare parasitemia between groups
ggplot(allvars_EvL, aes(x = season_time, y = log(parasitemiepf), fill = season_time)) + 
  geom_boxplot() + 
  theme_classic() + 
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 13) + 
  geom_line(aes(group=id.y), linetype="11") + 
  scale_x_discrete(name = "Season Timing", labels = c("Early Wet", "Late Wet")) +
  theme(text =element_text(size = 20)) + 
  ylab("Log parasitemia") + 
  scale_fill_manual(name = "Season Timing", labels = c("Early Wet", "Late Wet"), values = c("white", "grey"))

#Differential expression, UNADJUSTED
id <- as.factor(label_pairedEvL$id)
parasitemia <- as.numeric(log(allvars_EvL$parasitemiepf))

design_pairedEvL <- model.matrix(~dgList_pairedEvL_filtered$samples$group + id + parasitemia)
rownames(design_pairedEvL) <- colnames(dgList_pairedEvL_filtered)
design_pairedEvL

y <- estimateDisp(dgList_pairedEvL_filtered, design_pairedEvL)
fit <- glmFit(y, design_pairedEvL)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
summary(decideTests(lrt))

results_pairedEvL <- as.data.frame(lrt$table)
results_pairedEvL$FDR <- p.adjust(results_pairedEvL$PValue, method="fdr")

results_pairedEvL$deexpressed <- ifelse(results_pairedEvL$FDR <= 0.1 & results_pairedEvL$logFC >= 0, "UP", 
                                        ifelse(results_pairedEvL$FDR <= 0.1 & results_pairedEvL$logFC < 0, "DOWN", "NO"))
results_pairedEvL$gene <- rownames(results_pairedEvL)
sigEvL <- subset(results_pairedEvL, FDR <= 0.1)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
my_genes <- c("CCL5", "ADA", "GNLY", "FGFBP2", "IL18", "GBP1", "GBP4", "GBP5", "PARP14", "CLIC4", "LRRK2")

ggplot(results_pairedEvL) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.5, size = 2)) + 
  scale_colour_manual(values = mycolors) + 
  geom_label_repel(data = subset(results_pairedEvL, gene %in% my_genes), 
                   aes(logFC,-log10(PValue),label=gene),
                   size = 8,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + 
  theme_classic() + 
  theme(text = element_text(size = 30), legend.position = "none")

#Differential expression, ADJUSTED FOR IMMUNE CELL PROPORTIONS
cells <- read.delim("seasonality_cellprops.txt")
cells_pairedEvL <- cells[cells$Mixture %in% label_pairedEvL$full_id,]

colnames(cells_pairedEvL)[colnames(cells_pairedEvL) == 'Mixture'] <- 'full_id'
cells_pairedEvL <- merge(cells_pairedEvL, label_pairedEvL, by = "full_id")
cells_pairedEvL$id <- c("150", "150", "158A", "158B", "158B", "158A", "32", "32", "35", "35", "76", "76", "79", "79", "87", "87")

#Compare proportions of immune cell types
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(cells_pairedEvL, aes(x = season_time, y = Dendritic.cells.activated)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 1, fill = "darkgrey") + 
  geom_line(aes(group=id), linetype="11") + 
  scale_x_discrete(name = "Season Timing", labels = c("Early Wet", "Late Wet")) + 
  ylab("Proportion of Activated Dendritic Cells") + 
  theme_classic() + theme(text =element_text(size = 20)) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, color = "red") + 
  stat_compare_means(method = "t.test", paired = TRUE, size = 8, label.y = 0.06)

#Corrected for cell composition
cells_pairedEvL <- cells_pairedEvL[order(match(cells_pairedEvL$full_id, colnames(dgList_pairedEvL_filtered$counts))),]

cells_pairedEvL$bcells <- cells_pairedEvL$B.cells.naive + cells_pairedEvL$B.cells.memory + cells_pairedEvL$Plasma.cells
cells_pairedEvL$tcells <- cells_pairedEvL$T.cells.CD8 + cells_pairedEvL$T.cells.CD4.naive + cells_pairedEvL$T.cells.CD4.memory.resting + cells_pairedEvL$T.cells.CD4.memory.activated + cells_pairedEvL$T.cells.follicular.helper + cells_pairedEvL$T.cells.regulatory..Tregs. + cells_pairedEvL$T.cells.gamma.delta
cells_pairedEvL$NK <- cells_pairedEvL$NK.cells.activated + cells_pairedEvL$NK.cells.resting
cells_pairedEvL$APCs <- cells_pairedEvL$Monocytes + cells_pairedEvL$Macrophages.M0 + cells_pairedEvL$Macrophages.M1 + cells_pairedEvL$Macrophages.M2 + cells_pairedEvL$Dendritic.cells.activated + cells_pairedEvL$Dendritic.cells.resting
cells_pairedEvL$mast <- cells_pairedEvL$Mast.cells.activated + cells_pairedEvL$Mast.cells.resting
cells_pairedEvL$neut <- cells_pairedEvL$Neutrophils

bcells <- as.numeric(cells_pairedEvL$bcells)
tcells <- as.numeric(cells_pairedEvL$tcells)
NK <- as.numeric(cells_pairedEvL$NK)
mast <- as.numeric(cells_pairedEvL$mast)
neut <- as.numeric(cells_pairedEvL$neut)
APC <- as.numeric(cells_pairedEvL$APCs)

id <- as.factor(all_EvL$pid.x)
parasitemia <- as.numeric(log(all_EvL$parasitemiepf))

design_pairedEvL_adj <- model.matrix(~dgList_pairedEvL_filtered$samples$group + id + bcells + tcells + NK + APC + neut + parasitemia)

rownames(design_pairedEvL_adj) <- colnames(dgList_pairedEvL_filtered)
design_pairedEvL_adj

y <- estimateDisp(dgList_pairedEvL_filtered, design_pairedEvL_adj)
fit <- glmFit(y, design_pairedEvL_adj)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
summary(decideTests(lrt))

results_pairedEvL <- as.data.frame(lrt$table)
results_pairedEvL$FDR <- p.adjust(results_pairedEvL$PValue, method="fdr")

results_pairedEvL$deexpressed <- ifelse(results_pairedEvL$FDR <= 0.1 & results_pairedEvL$logFC >= 0, "UP", 
                                        ifelse(results_pairedEvL$FDR <= 0.1 & results_pairedEvL$logFC < 0, "DOWN", "NO"))
results_pairedEvL$gene <- rownames(results_pairedEvL)
sigEvL <- subset(results_pairedEvL, FDR <= 0.1)

#Volcano plot
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
mygenes_adj <- c("MYL9")


ggplot(results_pairedEvL) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.5, size = 2)) + 
  scale_colour_manual(values = mycolors) + 
  geom_label_repel(data = subset(results_pairedEvL, FDR <= 0.1), 
                   aes(logFC,-log10(PValue),label=gene),
                   size = 8,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + 
  theme_classic() + 
  theme(text = element_text(size = 30), legend.position = "none")

### Paired samples late first, early second ###
label_pairedLvE <- subset(label, include_paired_late_vs_early == "Y")
label_pairedLvE <- subset(label_pairedLvE, id != 295 & id != 35)
counts_pairedLvE <- counts[,colnames(counts) %in% label_pairedLvE$full_id]

dgList_pairedLvE <- DGEList(counts=counts_pairedLvE, genes=human_counts_all.txt$Geneid)

#Filter out lowly expressed genes 
keep <- rowSums(cpm(dgList_pairedLvE)>10) >= (ncol(counts_pairedLvE)/2)
table(keep)
dgList_pairedLvE_filtered <- dgList_pairedLvE[keep, , keep.lib.sizes=FALSE]

dim(dgList_pairedLvE_filtered)

dgList_pairedLvE_filtered <- calcNormFactors(dgList_pairedLvE_filtered, method="TMM")

label_pairedLvE <- label_pairedLvE[order(match(label_pairedLvE$full_id, colnames(dgList_pairedLvE_filtered$counts))),]
dgList_pairedLvE_filtered$samples$group <- label_pairedLvE$season_time

#Differential expression, unadjusted
all_LvE <- all[all$id %in% colnames(dgList_pairedLvE_filtered$counts),]
all_LvE <- all_LvE[order(match(all_LvE$id, colnames(dgList_pairedLvE_filtered))),]
all_LvE$full_id <- all_LvE$id
allvars_LvE <- merge(all_LvE, label_pairedLvE, by = "full_id")
allvars_LvE$season_time <- factor(allvars_LvE$season_time, levels = c("late_wet", "early_wet"))

id <- factor(label_pairedLvE$id)
parasitemia <- as.numeric(log(all_LvE$parasitemiepf))

design_pairedLvE <- model.matrix(~dgList_pairedLvE_filtered$samples$group + id + parasitemia)
rownames(design_pairedLvE) <- colnames(dgList_pairedLvE_filtered)
design_pairedLvE

y <- estimateDisp(dgList_pairedLvE_filtered, design_pairedLvE)
fit <- glmFit(y, design_pairedLvE)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

summary(decideTests(lrt))

results_pairedLvE <- as.data.frame(lrt$table)
results_pairedLvE$FDR <- p.adjust(results_pairedLvE$PValue, method="fdr")

results_pairedLvE$deexpressed <- ifelse(results_pairedLvE$FDR <= 0.1 & results_pairedLvE$logFC >= 0, "UP", 
                                        ifelse(results_pairedLvE$FDR <= 0.1 & results_pairedLvE$logFC < 0, "DOWN", "NO"))
results_pairedLvE$gene <- rownames(results_pairedLvE)
sigLvE <- subset(results_pairedLvE, FDR <= 0.1)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results_pairedLvE) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.5, size = 2)) + 
  scale_colour_manual(values = mycolors) + 
  geom_label_repel(data = subset(results_pairedLvE, FDR <= 0.1), 
                   aes(logFC,-log10(PValue),label=gene),
                   size = 8,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + 
  theme_classic() + 
  theme(text = element_text(size = 30), legend.position = "none") + ylim(0,7)

#Differences in immune cell composition
cells <- read.delim("seasonality_cellprops.txt")
cells_pairedLvE <- cells[cells$Mixture %in% label_pairedLvE$full_id,]
cells_pairedLvE <- cells_pairedLvE[,-c(24:26)]

colnames(cells_pairedLvE)[colnames(cells_pairedLvE) == 'Mixture'] <- 'full_id'
cells_pairedLvE <- merge(cells_pairedLvE, label_pairedLvE, by = "full_id")

cells_pairedLvE$season_time <- factor(cells_pairedLvE$season_time, levels = c("late_wet", "early_wet"))

ggplot(cells_pairedLvE, aes(x = season_time, y = B.cells.naive)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 1) + 
  geom_line(aes(group=id), linetype="11") + 
  scale_x_discrete(name = "Season Timing", labels = c("Early Wet", "Late Wet")) + 
  ylab("Proportion of Naive B Cells") + 
  theme_classic() + theme(text =element_text(size = 20)) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, color = "red") + 
  stat_compare_means(method = "t.test", paired = TRUE, size = 8, label.y = 0.105)

#### PLASMODIUM ####
pf_counts_all.txt <- read.delim("/Volumes/projects-t3/SerreDLab-3/kieran.tebben/mali/hisat2_outputs/count_tables/pf_counts_all.txt.gz", comment.char="#")
counts <- pf_counts_all.txt[,c(7:240)]

##Clean up colnames to be only sample ID
for ( col in 1:ncol(counts)){
  colnames(counts)[col] <-  sub(".pf.subset.rmdup.bam*", "", colnames(counts)[col])
}

for ( col in 1:ncol(counts)){
  colnames(counts)[col] <-  sub("*X", "", colnames(counts)[col])
}
rownames(counts) <- pf_counts_all.txt$Geneid

all <- read.csv("mali_samples_allvariables.csv")
all <- all[all$id %in% colnames(counts),]
all <- all[ order(match(all$id, colnames(counts))), ]
all <- all[-54,]

label <- read.delim("seasonality_variables.txt")
label <- label[ order(match(label$full_id, colnames(counts))), ]
label$season_time <- sub(" ", "_", label$season_time)

#Paired early first late second
label_pairedEvL <- subset(label, include_paired_early_vs_late == "Y")
counts_pairedEvL <- counts[,colnames(counts) %in% label_pairedEvL$full_id]

dgList_pairedEvL <- DGEList(counts=counts_pairedEvL, genes=rownames(counts_pairedEvL))

#Filter out lowly expressed genes 
keep <- rowSums(cpm(dgList_pairedEvL)>10) >= (ncol(counts_pairedEvL)/2)
table(keep)
dgList_pairedEvL_filtered <- dgList_pairedEvL[keep, , keep.lib.sizes=FALSE]

dim(dgList_pairedEvL_filtered)

dgList_pairedEvL_filtered <- calcNormFactors(dgList_pairedEvL_filtered, method="TMM")

label_pairedEvL <- label_pairedEvL[order(match(label_pairedEvL$full_id, colnames(dgList_pairedEvL_filtered$counts))),]
dgList_pairedEvL_filtered$samples$group <- label_pairedEvL$season_time
all_EvL <- all[all$id %in% label_pairedEvL$full_id,]

#Differential expression, unadjusted
id <- factor(label_pairedEvL$id)

all_EvL <- all_EvL[order(match(all_EvL$id, colnames(dgList_pairedEvL_filtered))),]
parasitemia <- as.numeric(log(all_EvL$parasitemiepf))

design_pairedEvL <- model.matrix(~dgList_pairedEvL_filtered$samples$group + id + parasitemia)
rownames(design_pairedEvL) <- colnames(dgList_pairedEvL_filtered)
design_pairedEvL

y <- estimateDisp(dgList_pairedEvL_filtered, design_pairedEvL)
fit <- glmFit(y, design_pairedEvL)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
summary(decideTests(lrt))

results_pairedEvL <- as.data.frame(lrt$table)
results_pairedEvL$FDR <- p.adjust(results_pairedEvL$PValue, method="fdr")

results_pairedEvL$deexpressed <- ifelse(results_pairedEvL$FDR <= 0.1 & results_pairedEvL$logFC >= 0, "UP", 
                                        ifelse(results_pairedEvL$FDR <= 0.1 & results_pairedEvL$logFC < 0, "DOWN", "NO"))
results_pairedEvL$gene <- rownames(results_pairedEvL)
sigEvL <- subset(results_pairedEvL, FDR <= 0.1)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
results_pairedEvL$label <- ifelse(results_pairedEvL$gene == "PF3D7_1478600", "PfPTP3", 
                                  ifelse(results_pairedEvL$gene == "PF3D7_0301700", "PfPTP7",
                                         ifelse(results_pairedEvL$gene == "PF3D7_1219000", "PfFRM2", "NA")))


ggplot(results_pairedEvL) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.5, size = 2)) + 
  scale_colour_manual(values = mycolors) + 
  geom_label_repel(data = subset(results_pairedEvL, label != "NA"), 
                   aes(logFC,-log10(PValue),label=label),
                   size = 8,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') + 
  theme_classic() + 
  theme(text = element_text(size = 30), legend.position = "none")

#Paired samples late first, early second
label_pairedLvE <- subset(label, include_paired_late_vs_early == "Y")
label_pairedLvE <- subset(label_pairedLvE, id != 295 & id != 35)
counts_pairedLvE <- counts[,colnames(counts) %in% label_pairedLvE$full_id]
rownames(counts_pairedLvE) <- pf_counts_all.txt$Geneid

dgList_pairedLvE <- DGEList(counts=counts_pairedLvE, genes=rownames(counts_pairedLvE))

#Filter out lowly expressed genes 
keep <- rowSums(cpm(dgList_pairedLvE)>10) >= (ncol(counts_pairedLvE)/2)
table(keep)
dgList_pairedLvE_filtered <- dgList_pairedLvE[keep, , keep.lib.sizes=FALSE]

dim(dgList_pairedLvE_filtered)

dgList_pairedLvE_filtered <- calcNormFactors(dgList_pairedLvE_filtered, method="TMM")

label_pairedLvE <- label_pairedLvE[order(match(label_pairedLvE$full_id, colnames(dgList_pairedLvE_filtered$counts))),]
dgList_pairedLvE_filtered$samples$group <- label_pairedLvE$season_time

#Differential expression, unadjusted
all_LvE <- all[all$id %in% colnames(dgList_pairedLvE_filtered$counts),]
all_LvE <- unique(all_LvE)
all_LvE <- all_LvE[order(match(all_LvE$id, label_pairedLvE$full_id)),]
id <- factor(label_pairedLvE$id)
age <- as.numeric(all_LvE$age)
parasitemia <- as.numeric(log(all_LvE$parasitemiepf))

design_pairedLvE <- model.matrix(~dgList_pairedLvE_filtered$samples$group + id + parasitemia)
rownames(design_pairedLvE) <- colnames(dgList_pairedLvE_filtered)
design_pairedLvE

design_pairedLvE_interaction <- model.matrix(~dgList_pairedLvE_filtered$samples$group + id + dgList_pairedLvE_filtered$samples$group:parasitemia)
rownames(design_pairedLvE_interaction) <- colnames(dgList_pairedLvE_filtered)
design_pairedLvE_interaction

y <- estimateDisp(dgList_pairedLvE_filtered, design_pairedLvE)
fit <- glmFit(y, design_pairedLvE)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

summary(decideTests(lrt))

results_pairedLvE <- as.data.frame(lrt$table)
results_pairedLvE$FDR <- p.adjust(results_pairedLvE$PValue, method="fdr")

results_pairedLvE$deexpressed <- ifelse(results_pairedLvE$FDR <= 0.1 & results_pairedLvE$logFC >= 0, "UP", 
                                        ifelse(results_pairedLvE$FDR <= 0.1 & results_pairedLvE$logFC < 0, "DOWN", "NO"))
results_pairedLvE$gene <- rownames(results_pairedLvE)
sigLvE <- subset(results_pairedLvE, FDR <= 0.1)

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")


ggplot(results_pairedLvE) + 
  geom_point(aes(logFC,-log10(PValue),col= deexpressed, alpha = 0.5, size = 2)) + 
  scale_colour_manual(values = mycolors) + 
  theme_classic() + 
  theme(text = element_text(size = 30), legend.position = "none")


