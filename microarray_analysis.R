###Script to analize microarray data
###Written by Angelina Volkova
### 09/04/2018
library(limma)
library(multtest)
options(scipen=999)
setwd("/Users/angelina_volkova/Desktop/Google_Drive/microarray_project/files/")
filenames <- list.files(pattern="*.txt")
file_list <- list()
###Read all the files
for (i in 1:length(filenames)) {
  file <- read.table(paste0("/Users/angelina_volkova/Desktop/Google_Drive/microarray_project/files/",filenames[i]), 
                     header=T, sep="\t")
  ###Calculate the mean of the duplicate genes in the file
  same_genes_mean <- aggregate(file[,-1], list(file$Gene.name.d), mean)
  file_list[[i]] <- same_genes_mean[-1,]
}
###Group data by microarrays
grouped_data <- list()
for (i in c(1,3,5,7)){
  group <- cbind(file_list[[i]],file_list[[i+1]])
  row.names(group) <- group$Group.1
  group <- group[,c(2:5,7:11)]
  grouped_data[[i]] <- group
}

###Quantile normalization of the arrays
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

micro_norm1 <- quantile_normalisation(grouped_data[[1]])
micro_norm2 <- quantile_normalisation(grouped_data[[3]])
micro_norm3 <- quantile_normalisation(grouped_data[[5]])
micro_norm4 <- quantile_normalisation(grouped_data[[5]])

###Find common genes between 4 arrays
common_genes <- Reduce(intersect, list(row.names(micro_norm1), row.names(micro_norm2),
                                        row.names(micro_norm3), row.names(micro_norm4)))

###Select the common genes from each of the data frames
micro_common1 <- as.data.frame(micro_norm1[row.names(micro_norm1)%in%c(common_genes),])
micro_common2 <- as.data.frame(micro_norm2[row.names(micro_norm2)%in%c(common_genes),])
micro_common3 <- as.data.frame(micro_norm3[row.names(micro_norm3)%in%c(common_genes),])
micro_common4 <- as.data.frame(micro_norm4[row.names(micro_norm4)%in%c(common_genes),])
micro_final <- cbind(micro_common1, micro_common2, micro_common3, micro_common4)
###Rearrange the columns
micro_final <- micro_final[,c(1,10,19,28,2,11,20,29,3,12,21,30,4,13,22,31,
                              5,14,23,32,6,15,24,33,7,16,25,34,8,17,26,35,
                              9,18,27,36)]
###Treat same cell lines as replicates andfind the DE genes
trts=factor(c(rep("CO.SW.620",4),rep("ME.MALME.3M",4), rep("OV.IGROV1",4), rep("OV.OVCAR.4", 4),
              rep("CO.HCT.116",4), rep("CO.HT29",4), rep("LE.SR",4), rep("LC.NCI.H226",4), rep("OV.SK.OV.3",4)))
design.trt=model.matrix(~0+trts)
fit_mean <- lmFit(micro_final, design.trt)
contrast.matrix=makeContrasts(trtsCO.SW.620-trtsCO.HCT.116, trtsCO.SW.620-trtsCO.HT29, trtsCO.SW.620-trtsLE.SR,
                              trtsCO.SW.620-trtsLC.NCI.H226, trtsCO.SW.620-trtsOV.SK.OV.3,
                              trtsME.MALME.3M-trtsCO.HCT.116, trtsME.MALME.3M-trtsCO.HT29, trtsME.MALME.3M-trtsLE.SR,
                              trtsME.MALME.3M-trtsLC.NCI.H226, trtsME.MALME.3M-trtsOV.SK.OV.3,
                              trtsOV.IGROV1-trtsCO.HCT.116, trtsOV.IGROV1-trtsCO.HT29, trtsOV.IGROV1-trtsLE.SR,
                              trtsOV.IGROV1-trtsLC.NCI.H226, trtsOV.IGROV1-trtsOV.SK.OV.3,
                              trtsOV.OVCAR.4-trtsCO.HCT.116, trtsOV.OVCAR.4-trtsCO.HT29, trtsOV.OVCAR.4-trtsLE.SR,
                              trtsOV.OVCAR.4-trtsLC.NCI.H226, trtsOV.OVCAR.4-trtsOV.SK.OV.3,
                              levels=design.trt)
fit.contrast=contrasts.fit(fit_mean,contrast.matrix)
efit.contrast=eBayes(fit.contrast)
comparisons <- list()
for (i in 1:20){
  comp <- topTable(efit.contrast,coef=i,adjust.method="BH",p.value=0.05,genelist=common_genes, number=500)
  comparisons[[i]] <- comp
}
###List of significant genes
results <- as.data.frame(decideTests(efit.contrast))
results_not_0 <- results[!rowSums(results)==0,]
sig_genes <- c(row.names(results_not_0))
final_logfc <- merge(comparisons[[1]][,c(1,2)], comparisons[[2]][,c(1,2)], by="ID", all=TRUE)
for (i in 3:20) {
  final_logfc <- merge(final_logfc, comparisons[[i]][,c(1,2)], by="ID", all=TRUE)
}

names(final_logfc) <- c("Genes", "CO.SW.620_vs_CO.HCT.116", "CO.SW.620_vs_CO.HT29", "CO.SW.620_vs_LE.SR",
                        "CO.SW.620_vs_LC.NCI.H226", "CO.SW.620_vs_OV.SK.OV.3",
                        "ME.MALME.3M_vs_CO.HCT.116", "ME.MALME.3M_vs_CO.HT29","ME.MALME.3M_vs_LE.SR",
                        "ME.MALME.3M_vs_LC.NCI.H226", "ME.MALME.3M_vs_OV.SK.OV.3",
                        "OV.IGROV1_vs_CO.HCT.116", "OV.IGROV1_vs_CO.HT29", "OV.IGROV1_vs_LE.SR",
                        "OV.IGROV1_vs_LC.NCI.H226", "OV.IGROV1_vs_OV.SK.OV.3",
                        "OV.OVCAR.4_vs_CO.HCT.116", "OV.OVCAR.4_vs_CO.HT29", "OV.OVCAR.4_vs_LE.SR",
                        "OV.OVCAR.4_vs_LC.NCI.H226", "OV.OVCAR.4_vs_OV.SK.OV.3")
final_logfc <- final_logfc[final_logfc$Genes%in%c(sig_genes),]
write.table(final_logfc, "/Users/angelina_volkova/Desktop/microarray_results.txt", sep="\t", quote=F, col.names=T, row.names = F)


