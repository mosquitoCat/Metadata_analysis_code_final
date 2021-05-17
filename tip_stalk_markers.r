### Top-50 marker genes of pooled traditional EC subtypes heart from Peter Carmeliet's Cell paper "Single-Cell Transcriptome Atlas of Murine Endothelial cells"
marker.file <- read.csv("PC_Cell_Paper_markers.csv", header = TRUE)
# since these are murine markers, I pulled out the corresponding human genes
artery.marker <- marker.file$artery
capillary.marker <- marker.file$capillary
vein.marker <- marker.file$vein

tip.marker <- c("Adm", "Ankrd37", "C1qtnf6", "Cldn5", "Col4a1", "Cotl1", "Dll4", "Ednrb", "Fscn1", "Gpihbp1", "Hspg2", "Igfbp3", "Inhbb", "Jup", "Kcne3", "Kcnj8", "Lama4", "Lamb1", "Lxn", "Marcksl1", "Mcam", "Mest", "N4bp3", "Nid2", "Notch4", "Plod1", "Plxnd1", "Pmepa1", "Ptn", "Ramp3", "Rbp1", "Rgcc", "Rhoc", "Trp53i11")

stalk.marker <- c("Ackr1", "Aqp1", "C1qtnf9", "Cd36", "Csrp2", "Ehd4", "Fbln5", "Hspb1", "Ligp1", "Il6st", "Jam2", "Lgals3", "Lrg1", "Meox2", "Plscr2", "Sdpr", "Selp", "Spint2", "Tgfbi", "Tgm2", "Tmem176a", "Tmem176b", "Tmem252", "Tspan7", "Vwf")

# Load count data
exprMatrix <- integrated@assays$RNA@data
exprMatrix <- as.matrix(exprMatrix) # convert sparse matrix to matrix
dim(exprMatrix)

# Set gene set
tip.geneSets <- GeneSet(as.character(tip.marker), setName = "Tip_ECs")
tip.geneSets

stalk.geneSets <- GeneSet(as.character(stalk.marker), setName = "Stalk_ECs")
stalk.geneSets

artery.geneSets <- GeneSet(as.character(artery.marker), setName = "Arterial_ECs")
artery.geneSets

capillary.geneSets <- GeneSet(as.character(capillary.marker), setName = "Capillary_ECs")
capillary.geneSets

vein.geneSets <- GeneSet(as.character(vein.marker), setName = "Venous_ECs")
vein.geneSets

# Build gene expression rankings for each cell
detectCores()
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 2, plotStats = TRUE)
cells_rankings
save(cells_rankings, file = "EC_rankings.RData")

# Calculate enrichment for the gene signatures (AUC)
T.cells_AUC <- AUCell_calcAUC(tip.geneSets, cells_rankings)
S.cells_AUC <- AUCell_calcAUC(stalk.geneSets, cells_rankings)
A.cells_AUC <- AUCell_calcAUC(artery.geneSets, cells_rankings)
C.cells_AUC <- AUCell_calcAUC(capillary.geneSets, cells_rankings)
V.cells_AUC <- AUCell_calcAUC(vein.geneSets, cells_rankings)

# add AUC score to the MI object
integrated@meta.data$tip_AUC_score <- t(getAUC(T.cells_AUC))
save(T.cells_AUC, file = "tip_EC_AUC.RData")

integrated@meta.data$stalk_AUC_score <- t(getAUC(S.cells_AUC))
save(S.cells_AUC, file = "stalk_EC_AUC.RData")

integrated@meta.data$artery_AUC_score <- t(getAUC(A.cells_AUC))
save(A.cells_AUC, file = "artery_EC_AUC.RData")

integrated@meta.data$vein_AUC_score <- t(getAUC(V.cells_AUC))
save(V.cells_AUC, file = "vein_EC_AUC.RData")

integrated@meta.data$capillary_AUC_score <- t(getAUC(C.cells_AUC))
save(C.cells_AUC, file = "capillary_EC_AUC.RData")

## integrated
integrated.df <- cbind(integrated@reductions$umap@cell.embeddings, integrated@meta.data)
# [1] "#F8766D" "#E58700" "#C99800" "#A3A500" "#6BB100"
# [6] "#00BA38" "#00BF7D" "#00C0AF" "#00BCD8" "#00B0F6"
# [11] "#619CFF" "#B983FF" "#E76BF3" "#FD61D1" "#FF67A4"
integrated.colour <- c("#F8766D", "#E58700", "#C99800", "#A3A500", "#6BB100", "#00BA38", "#00BF7D", "#00C0AF", "#00BCD8", "#00B0F6", "#619CFF", "#B983FF", "#E76BF3", "#FD61D1", "#FF67A4")

integrated.p <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ seurat_clusters, colors = integrated.colour, marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.tip <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(tip_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.stalk <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(stalk_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.artery <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(artery_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.vein <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(vein_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.capillary <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(capillary_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()
