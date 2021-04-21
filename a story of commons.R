# Cross-species comparison injured vs uninjured

# Common up genes between Injured_vs_Uninjured human and mouse
human_InjuredVSUninjuredAllDEG_up <- homologene(AllDEG_up$X, inTax = 9606, outTax = 10090) [, 2]
save(human_InjuredVSUninjuredAllDEG_up, file = "human_InjuredVSUninjuredAllDEG_up.RData")
Mouse_InjuredVSUninjuredAllDEG_up <- AllDEG_up$X
save(Mouse_InjuredVSUninjuredAllDEG_up, file = "Mouse_InjuredVSUninjuredAllDEG_up.RData")
common_up <- intersect(human_InjuredVSUninjuredAllDEG_up, Mouse_InjuredVSUninjuredAllDEG_up)
common_up_for_human <- unique(homologene(common_up, inTax = 10090, outTax = 9606) [, 2])
save(common_up, file = "common_upregulated_DEG_in_injured_compared_to_uninjured_between_MandH.RData")

# Common down genes between Injured_vs_Uninjured human and mouse
human_InjuredVSUninjuredAllDEG_down <- homologene(AllDEG_down$X, inTax = 9606, outTax = 10090) [, 2]
save(human_InjuredVSUninjuredAllDEG_down, file = "human_InjuredVSUninjuredAllDEG_down.RData")
Mouse_InjuredVSUninjuredAllDEG_down <- AllDEG_down$X
save(Mouse_InjuredVSUninjuredAllDEG_down, file = "Mouse_InjuredVSUninjuredAllDEG_down.RData")
common_down <- intersect(human_InjuredVSUninjuredAllDEG_down, Mouse_InjuredVSUninjuredAllDEG_down)
common_down_for_human <- unique(homologene(common_down, inTax = 10090, outTax = 9606) [, 2])
save(common_down, file = "common_downregulated_DEG_in_injured_compared_to_uninjured_between_MandH.RData")

# Common up regulons between Injured_vs_Uninjured human and mouse
## human - UninjuredInjuredAllDER_up
human_UninjuredInjuredAllDER_up <- homologene(UninjuredInjuredAllDER_up$X, inTax = 9606, outTax = 10090)[, 2]
save(human_UninjuredInjuredAllDER_up, file = "human_UninjuredInjuredAllDER_up.RData")
Mouse_UninjuredInjuredAllDER_up <- rownames(AllDER_up)
save(Mouse_UninjuredInjuredAllDER_up, file = "Mouse_UninjuredInjuredAllDER_up.RData")
common_up <- intersect(human_UninjuredInjuredAllDER_up, Mouse_UninjuredInjuredAllDER_up)
save(common_up, file = "common_upregulated_DER_in_injured_compared_to_uninjured_between_MandH.RData")

# Common down regulons between Injured_vs_Uninjured human and mouse
## human - UninjuredInjuredAllDER_down
human_UninjuredInjuredAllDER_down <- homologene(UninjuredInjuredAllDER_down$X, inTax = 9606, outTax = 10090)[, 2]
save(human_UninjuredInjuredAllDER_down, file = "human_UninjuredInjuredAllDER_down.RData")
Mouse_UninjuredInjuredAllDER_down <- rownames(AllDER_down)
save(Mouse_UninjuredInjuredAllDER_down, file = "Mouse_UninjuredInjuredAllDER_down.RData")
common_down <- intersect(human_UninjuredInjuredAllDER_down, Mouse_UninjuredInjuredAllDER_down)
save(common_down, file = "common_downregulated_DER_in_injured_compared_to_uninjured_between_MandH.RData")

# rownames(integrated)
new_human_gene_names <- homologene(rownames(integrated), inTax = 9606, outTax = 10090) [, 2]
length(new_human_gene_names)
mouse_gene_names <- rownames(mouse_integrated)
common_gene_names <- intersect(new_human_gene_names, mouse_gene_names)
length(common_gene_names)

###################################
## Find one-to-one match homologene
###################################

mouse_genes <- rownames(mouse_integrated)
human_genes <- rownames(human_integrated)

# define %!in% function
'%!in%' <- function(x, y) {!'%in%'(x, y)}

human_to_mouse_full_list <- human2mouse(human_genes)
mouse_table <- dplyr::count(human_to_mouse_full_list, mouseGene)
mouse_to_exclude <- mouse_table %>% dplyr::filter (n > 1)
final_list <- human_to_mouse_full_list %>% dplyr::filter (mouseGene %!in% mouse_to_exclude$mouseGene)

# reverse check 
mouse_to_human_full_list <- mouse2human((final_list$mouseGene)) 
human_table <- dplyr::count(human_to_mouse_full_list, humanGene)
human_to_exclude <- human_table %>% dplyr::filter (n > 1)
final_list1 <- final_list %>% dplyr::filter (humanGene %!in% human_to_exclude$humanGene)
# new human object only contains the one-to-one genes and to make the integration work, the new human dataset will be using the mouse gene format
gene.use <- final_list1$humanGene

matrix <- human_integrated@assays$RNA@counts[gene.use, ]
new_rownames <- human2mouse(rownames(matrix))$mouseGene 
rownames(matrix) <- new_rownames
human_counts <- as.data.frame(matrix)
human_counts <- human_counts %>% filter(rownames(human_counts) %in% common_gene_names)
human_meta <- human_integrated@meta.data %>% filter (rownames(human_integrated@meta.data) %in% colnames(human_counts))
newhuman_meta <- human_meta[, c(4:8, 11, 43)]
newhuman_meta1 <- newhuman_meta[, 2:3]
newhuman_meta1$detailedage <- newhuman_meta1$age
newhuman_meta1$condition <- newhuman_meta$newcondition
newhuman_meta1$detailedcondition <- newhuman_meta$condition
newhuman_meta1$source <- newhuman_meta$AccessionNo
newhuman_meta1$group <- newhuman_meta$condition

trimmed_human <- CreateSeuratObject(human_counts, project = "human", meta.data = newhuman_meta1)
trimmed_human <- NormalizeData(trimmed_human)

#########################
### Merge human and mouse 
#########################
## only common genes are used
## human and mouse objects are trimmed to include only the common genes
# rownames(integrated)
new_human_gene_names <- homologene(rownames(integrated), inTax = 9606, outTax = 10090) [, 2]
length(new_human_gene_names)
mouse_gene_names <- rownames(mouse_integrated)
common_gene_names <- intersect(new_human_gene_names, mouse_gene_names)
length(common_gene_names)

# building mouse
mouse_counts <- as.data.frame(mouse_integrated@assays$RNA@counts)
mouse_counts <- mouse_counts %>% filter(rownames(mouse_counts) %in% common_gene_names)
mouse_meta <- mouse_integrated@meta.data %>% filter (rownames(mouse_integrated@meta.data) %in% colnames(mouse_counts))
newmouse_meta <- mouse_meta[, c(4:9, 12, 34)]
colnames(newmouse_meta)
#[1] "species"           "age"              
#[3] "detailedage"       "condition"        
#[5] "detailedcondition" "source"           
#[7] "group"  

trimmed_mouse <- CreateSeuratObject(counts = mouse_counts, project = "mouse", meta.data = newmouse_meta)
trimmed_mouse <- NormalizeData(trimmed_mouse)

# put the two objects in a list
list <- c(trimmed_mouse, trimmed_human)

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)

# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:50)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)

for (i in 1:20) {
  combined <- FindClusters(combined, resolution = 0.1*i)
}

png("clustree.png", res = 300, width = 4000, height = 3000)
clustree(combined) # cluster tree visualisation
dev.off()

png("clustree_stability.png", res = 300, width = 4000, height = 3000)
clustree(combined, node_colour = "sc3_stability") # cluster tree visualisation
dev.off()

combined <- FindClusters(combined, resolution = 0.3)
saveRDS(combined, file = "combined.rds")

# It's very important to have an overall view of the dataset and for this purpose, I build a new data frame using the UMAP coordinates and the metadata.
df <- cbind(combined@reductions$umap@cell.embeddings, combined@meta.data)
# Plotly is used to generate beautiful feature plot here.
t <- list(
  family = "Open Sans",
  size = 16,
  color = 'black')

#################
# plot by species
#################
col.touse1 <- c("#052955", "#A7B8F8")
p1 <- plot_ly(data = df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ species, colors = col.touse1, marker = list(size = 12, line = list(color = "black", width = 1)), width = 1000, height = 1000) %>% layout(title = "", yaxis = list(zeroline = FALSE), xaxis = list(zeroline = FALSE), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 40))) %>% hide_colorbar()%>% hide_legend()
## prepre data for pie chart
df1 <- df %>% group_by(species) %>% summarise(count = n()) %>% mutate (percentage = count/nrow(df)*100)
fig1 <- plot_ly(df1, x = ~ species, y = ~ count, type = 'bar', marker = list(color = col.touse1, line = list(color = "black", width = 1)))
fig1_pie <- plot_ly(df1, labels = ~ species, values = ~ percentage, type = 'pie', marker = list(colors = col.touse1, line = list(color = '#FFFFFF', width = 1)), textfont = list(size = 40), insidetextorientation = "radial") %>% layout(width = 500, height = 500, uniformtext = list(minsize = 12, mode ='hide')) %>% hide_legend()

#################
# plot by species
#################

# COLOUR BUILT #
p <- DimPlot(combined, label = TRUE, label.size = 6) + NoLegend()
# isolate colour used in the Seurat package
colourbuilt <- ggplot_build(p)[[1]][[1]][,  c(1, 5)]
# reduce to unique colour and cluster information
further <- distinct(colourbuilt)
further <- further[order(further$group), ]
further.colour <- further$colour

p2 <- plot_ly(data = df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ seurat_clusters, colors = further$colour, marker = list(size = 12, line = list(color = "black", width = 1)), width = 1000, height = 1000) %>% layout(title = "", yaxis = list(zeroline = FALSE), xaxis = list(zeroline = FALSE), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 40))) %>% hide_colorbar()%>% hide_legend()
## prepre data for pie chart
df2 <- df %>% group_by(seurat_clusters) %>% summarise(count = n()) %>% mutate (percentage = count/nrow(df)*100)
fig2 <- plot_ly(df2, x = ~ seurat_clusters, y = ~ count, type = 'bar', marker = list(color = further.colour, line = list(color = "black", width = 1)))
fig2_pie <- plot_ly(df2, labels = ~ seurat_clusters, values = ~ percentage, type = 'pie', marker = list(colors = further.colour, line = list(color = '#FFFFFF', width = 1)), textfont = list(size = 40), insidetextorientation = "radial") %>% layout(width = 500, height = 500, uniformtext = list(minsize = 12, mode ='hide')) %>% hide_legend()

#################
# plot by age
#################
age.colours <- c("orange", "purple", "green")

p3 <- plot_ly(data = df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ age, colors = age.colours, marker = list(size = 12, line = list(color = "black", width = 1)), width = 1000, height = 1000) %>% layout(title = "", yaxis = list(zeroline = FALSE), xaxis = list(zeroline = FALSE), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 40))) %>% hide_colorbar()%>% hide_legend()

df3 <- df %>% group_by(age) %>% summarise(count = n()) %>% mutate (percentage = count/nrow(df)*100)
fig3 <- plot_ly(df3, x = ~ age, y = ~ count, type = 'bar', marker = list(color = age.colours, line = list(color = "black", width = 1)))
fig3_pie <- plot_ly(df3, labels = ~ age, values = ~ percentage, type = 'pie', marker = list(colors = age.colours, line = list(color = '#FFFFFF', width = 1)), textfont = list(size = 40), insidetextorientation = "radial") %>% layout(width = 500, height = 500, uniformtext = list(minsize = 12, mode ='hide')) %>% hide_legend()

################
################
##### plot #####
################
################
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

combined@meta.data <- combined@meta.data %>% mutate(newgroup = case_when(group == "P6 Healthy" ~ "Mouse P6 Uninjured", group == "P10 Healthy" ~ "Mouse P10 Uninjured", group == "Adult Healthy" ~ "Mouse Adult Healthy", group == "P2 MI3W" ~ "Mouse P2 MI3W", group == "Adult MI1D" ~ "Mouse Adult MI1D", group == "Adult MI3D" ~ "Mouse Adult MI3D", group == "Adult MI7D" ~ "Mouse Adult MI7D", group == "Adult MI2W" ~ "Mouse Adult MI2W", group == "Adult MI4W" ~ "Mouse Adult MI4W", group == "Normal" & age == "Fetal" ~ "Human Fetal Uninjured", group == "Normal" & age == "Adult" & species == "Human"  ~ "Human Adult Uninjured", group == "cHF" ~ "Human Adult cHF", group == "dHF" ~ "Human Adult dHF"))

combined@meta.data$newgroup <- factor(combined@meta.data$newgroup, levels = c("Mouse P6 Uninjured", "Mouse P10 Uninjured", "Mouse P2 MI3W", "Mouse Adult Uninjured", "Mouse Adult MI1D", "Mouse Adult MI3D", "Mouse Adult MI7D", "Mouse Adult MI2W", "Mouse Adult MI4W", "Human Fetal Uninjured", "Human Adult Uninjured", "Human Adult cHF", "Human Adult dHF"))

common_up_DEG <- c("Jun", "Sparcl1", "Cd74", "Hspa1a", "Dnajb1", "Itga6", "Fosb", "Egr1", "B2m", "Zfp36", "H2-T23", "Ddx3y",  "Timp3",  "Rps4x",  "Klf4", "Apold1", "H2-K1", "H2-D1", "H2-Q4", "Hsp90aa1", "Hes1", "Adamts1", "Slfn5", "Cdkn1a", "Ier5", "Jund", "Fos", "Cebpd", "Zfp36l2", "Nfkbia", "Dnaja1", "Gimap4", "Nedd9", "Lgals3bp", "Sgk1", "Klf2", "Gimap6", "Tinagl1", "Rhob", "Tsc22d1", "Ets1", "Gas6", "Tsc22d3", "Tmod3", "Ubc", "H3f3b", "H3f3a")

common_down_DEG <- c("Col3a1", "Aes", "Calcrl", "Ckap4", "Fhl2", "Rtn3", "Lims1", "Lama4", "Plscr4", "Lamc1", "Lamb1", "Uqcr11", "Nid1")

common_up_DER <- c("Nfkb1", "Bach1", "Stat1", "Twist1", "Rela", "Ets2", "Klf4", "Foxo4", "Smad4", "Irf1", "Egr1", "Atf2", "Hsf1", "Vdr", "Nfic", "Elk1", "E2f1", "Stat3", "Fli1", "Arntl", "Wt1")

common_down_DER <- c("Tfap2c", "Pax6", "Foxa2", "Klf13", "Gata3")

c1 <- DotPlot(combined, features = common_up_DEG, group.by = "newgroup", cols = c("navy", "plum")) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

c2 <- DotPlot(combined, features = common_down_DEG, group.by = "newgroup", cols = c("navy", "plum")) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

c3 <- DotPlot(combined, features = common_up_DER, group.by = "newgroup", cols = c("navy", "plum"), assay = "dorothea") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

c4 <- DotPlot(combined, features = common_down_DER, group.by = "newgroup", cols = c("navy", "plum"), assay = "dorothea") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## Run dorothea

# load library
library(dorothea)
library(viper)
## read Dorothea Regulons for mouse:
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))

## obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_mouse %>% dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores 
combined <- run_viper(combined, regulon, options = list(method = "scale", minsize = 4, eset.filter = FALSE, cores = 2, verbose = FALSE))

## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = combined) <- "dorothea"
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = rownames(combined), verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:30, verbose = FALSE)

#########################
#########################
# find neighbours
for (i in 1:20) {
  combined <- FindClusters(combined, resolution = 0.1*i)
}

library(clustree)
png("dorothea_integrated1_clustree.png", res = 300, width = 4000, height = 3000)
clustree(combined) # cluster tree visualisation
dev.off()

png("dorothea_integrated1_clustree_stability.png", res = 300, width = 4000, height = 3000)
clustree(combined, node_colour = "sc3_stability") # cluster tree visualisation
dev.off()

# resolution 0.4 is a good choice
combined@meta.data$seurat_clusters <- combined@meta.data$dorothea_snn_res.0.3
combined <- FindClusters(combined, resolution = 0.3) # run it again
saveRDS(combined, file = "combined.rds")
#########################
#########################

set.seed(1024)
combined <- RunUMAP(combined, dims = 1:30, umap.method = "uwot", metric = "cosine")

########################
### get that plot up ###
########################
## DEG
DefaultAssay(combined) <- "RNA"
list1 <- c("Sdpr", "Slc28a2", "Rps28", "Hspa8")
list2 <- c("Cebpd", "Elob", "Gas5", "Rack1", "Selenof", "Selenok", "Selenom", "Selenow", "Sem1", "Atp6v0c", "Cavin2", "Vps28", "Nme2", "Lgals1", "Socs3", "Tmem176b", "Emp1", "Capg", "Ndflp1", "Apoe", "Atp1b3", "Rps28", "Ier3", "Rplp0", "Pabpc1", "S100a6", "Gadd45g", "Vim", "Atp1a1", "Dnajb1", "Ctsz", "Hmgn1", "Junb", "Actn1", "Nucb1", "Zfp36", "Sat1", "Fxyd5", "Itm2c", "Nmt1", "Fos", "Hspa8", "Btg1", "Slc38a2", "Hes1", "Ackr3", "lmna")
DotPlot(combined, features = list1, group.by = "newgroup", cols = c("navy", "plum"), assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## Check through all other figures and plots
## figure 1
list1 <- c("Btg1", "Cebpd", "Elob", "Fos", "Fxyd5", "Gas5", "H2afz", "Junb", "Nme2", "Rbm3")
list2 <- c("Actc1", "Aes", "Fabp3", "Gltscr2", "Mb", "Myh6")
DotPlot(combined, features = list1, group.by = "newgroup", cols = c("navy", "plum"), assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## figure 2
list1 <- c("Xist", "Egfl7", "Tmsb10", "Fabp5", "Pabpc1", "Fabp4", "Rps9", "Rpl41", "Mest")
list2 <- c("Ttn", "Ryr2", "Zbtb20", "Cacna1c", "Atp2a2", "Fhl2", "Slc8a1", "Rora", "Dmd", "Pcdh7")
list3 <- c("Actc1", "Actn2", "Airn", "Ankrd1", "Atp2a2", "Cacna1c", "Celf", "Col1a2", "Col3a1", "Ctnna3")
list4 <- c("AWW112010", "Actg1", "B2m", "Eef1a1", "Fabp4", "Fau", "Ftl1", "H2-D1", "H2-K1")
list5 <- c("Actc1", "Airn", "Ankrd1", "Fabp3", "Mb", "Myh6", "Myl2", "Myl3", "Tnnc1", "Tnnt2")
list6 <- c("Actg1", "B2m", "Fau", "Ftl1", "H2-D1", "H2-K1", "Hspb1", "Ifitm3")
list7 <- c("Egfl7", "Igf2", "H19", "Pabpc1", "Smc4", "Gnas", "Mest", "Sem1", "Actb", "Gas5")
list8 <- c("Mylip", "Anks1b", "Asxl3", "Chl1", "Dlgap1", "Nkain2", "Dgkb", "Rps27rt")
list9 <- c("Aes", "Cxx1a", "Cyr61", "Gltscr2")
list10 <- c("Cct6a", "Cops9", "Elob", "Gas5", "Ndfip1", "Nme2", "Rack1", "Rflnb", "Selenof")
DotPlot(combined, features = list1, group.by = "newgroup", cols = c("navy", "plum"), assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## figure 3
list1 <- c("Cebpd", "Elob", "Gas5", "Rack1", "Selenof", "Selenok", "Selenom", "Selenop", "Selenow", "Sem1", "Atp6v0c", "Cavin2", "Vps28", "Nme2", "Lgals1", "Socs3", "Tmem176b", "Emp1", "Capg", "Ndfip1", "Apoe", "Atp1b3", "Rps28", "Ier3", "Rplp0", "Pabpc1", "S100a6", "Gadd45g", "Vim", "Atp1a1", "Dnajb1", "Ctsz", "Hmgn1", "Junb", "Actn1", "Nucb1", "Zfp36", "Sat1", "Fxyd5", "Itm2c", "Nmt1", "Fos", "Hspa8", "Btg1", "Slc38a2", "Hes1", "Ackr3", "Lmna")
list2 <- c("Sdpr", "Slc28a2", "Rps28", "Hspa8")
DotPlot(combined, features = list1, group.by = "newgroup", cols = c("navy", "plum"), assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## figure 5
list1 <- c("Fth1", "Gapdh", "Eef1a1", "Actb", "Rps26", "Hbg2", "Rps23", "Rpl10", "Tuba1b", "Rpl30")
DotPlot(combined, features = list1, group.by = "newgroup", cols = c("navy", "plum"), assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## figure 6
list1 <- c("Hbg2", "Eef1a1", "Rpl10", "Hba1", "Hba2", "Mdk", "Rpl15", "Rps23", "Prs26", "Gapdh", "Myh6", "Nppa", "Ftl", "Fth1", "Ccl2", "Itln1", "Smarca4", "Pla2g2a", "Apod", "Cxcl8", "B2m", "Myl12A", "Myl3", "Acta1", "Actb")
list2 <- c("Ccl2", "Rpl39", "Rpl11", "Ptma", "Rps11", "Fn1", "Nap1l1", "Efnb2", "Rpl10ap6", "Ltbp4", "Ccdc80", "B2m", "Acta1", "Apod", "Mtrnr2l12", "Smarca4", "Msrb3", "Fth1")
DotPlot(combined, features = list1, group.by = "newgroup", cols = c("navy", "plum"), assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

speciallist <- c("Rps9", "Rps28", "Klf4", "Lgals1", "Pabpc1", "Plvap")
DotPlot(combined, features = speciallist, group.by = "newgroup", cols = c("navy", "plum"), assay = "RNA") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

### regulon
## figure 7
DefaultAssay(combined) <- "dorothea"
list1 <- c("Rfx", "nR2F2", "Pdx1", "Ppara", "Tp73", "Pou5F1", "Nfkb1", "Bach1", "Stat1", "Irf4")
list2 <- c("Myc", "Znf263", "Cebpd", "Tfap2c", "Bcl6", "Foxp1", "Maf", "Zkscan1", "Gata1", "Pax6")
DotPlot(combined, assay = "dorothea", features = list1, group.by = "newgroup", cols = c("navy", "plum")) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## figure 4
list1 <- c("Usf2", "Ssrp1", "Myc", "Usf1", "Rel", "Rela", "Trp53", "Foxj2", "Atf2", "Sp1", "Etv4", "Ets2", "Jun", "Egr1", "Gata2", "Nanog", "Ehf", "Rfx1", "Foxo1", "Tcf7l2")
DotPlot(combined, assay = "dorothea", features = list1, group.by = "newgroup", cols = c("navy", "plum")) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## CME genes Myl2, Mb, Myl3, Tnnt2, Tnni3, Actc1
cmelist <- c("Mb", "Myl2", "Myl3", "Tnni3", "Actc1", "Tnnt2", "Pln", "Myh6", "Tnnc1", "Ttn", "Actn2", "Ryr2", "Myom2", "Myl4", "Myl7", "Sln")
DotPlot(combined, assay = "RNA", features = cmelist, group.by = "newgroup", cols = c("navy", "plum")) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

## subset combined dataset for Mairi's talk
#(1) Mouse ‘regenerative window’ uninjured 
#(2) Mouse ‘regenerative window’ injured [i.e. exclude the mouse P10 timepoint and all adult timepoints for this purpose] 
#(3) human fetal uninjured 
#(4) human adult uninjured 
#(5) human HF [where cHF and dHF are grouped together for this presentation purpose] 
df <- subset(combined, subset = newgroup == c("Human Adult dHF", "Human Adult cHF", "Human Adult Uninjured", "Human Fetal Uninjured", "Mouse P2 MI3W", "Mouse P6 Uninjured"))
df@meta.data <- df@meta.data %>% mutate(artificialgroup = case_when(newgroup == "Human Adult dHF" ~ "Human HF", newgroup == "Human Adult cHF" ~ "Human HF", newgroup == "Human Adult Uninjured" ~ "Human adult uninjured", newgroup == "Human Fetal Uninjured" ~ "Human fetal uninjured", newgroup == "Mouse P2 MI3W" ~ "Mouse 'regenerative' window injured", newgroup == "Mouse P6 Uninjured" ~ "Mouse 'regenerative' window uninjured"))
df@meta.data$artificialgroup <- factor(df@meta.data$artificialgroup, levels = c("Human HF", "Human adult uninjured", "Human fetal uninjured", "Mouse 'regenerative' window injured", "Mouse 'regenerative' window uninjured"))
DotPlot(df, assay = "RNA", features = cmelist, group.by = "artificialgroup", cols = c("navy", "plum")) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

# plotly heatmap
library(plotly)

