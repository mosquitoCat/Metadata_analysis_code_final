
# Load count data
exprMatrix <- human_integrated@assays$RNA@data
exprMatrix <- as.matrix(exprMatrix) # convert sparse matrix to matrix
dim(exprMatrix)
# vessel type marker
marker <-read.csv("vessel_type_marker.csv")
new_marker <- marker %>% select(c(X, X.1, X.2, X.3))
new_marker <- new_marker[4:53, ]
colnames(new_marker) <- c("artery", "capillary", "vein", "lymphatic")
head(new_marker)
## mouse marker to human homologene
human_marker_artery <- mouse2human(new_marker[, 1])[2][, 1]
human_marker_capillary <- mouse2human(new_marker[, 2])[2][, 1]
human_marker_vein <- mouse2human(new_marker[, 3])[2][, 1]

# Set gene set
#artery_marker <- new_marker$artery
#capillary_marker <- new_marker$capillary
#vein_marker <- new_marker$vein

arterySets <- GeneSet(human_marker_artery, setName = "artery ECs")
capillarySets <- GeneSet(human_marker_capillary, setName = "capillary ECs")
veinSets <- GeneSet(human_marker_vein, setName = "vein ECs")

arterySets
capillarySets
veinSets

# Build gene expression rankings for each cell
detectCores()
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 2, plotStats = TRUE)
cells_rankings
save(cells_rankings, file = "EC_rankings_human_integrated.RData")

# Calculate enrichment for the gene signatures (AUC)
artery_AUC <- AUCell_calcAUC(arterySets, cells_rankings)
human_integrated@meta.data$artery_AUC_score <- t(getAUC(artery_AUC))

capillary_AUC <- AUCell_calcAUC(capillarySets, cells_rankings)
human_integrated@meta.data$capillary_AUC_score <- t(getAUC(capillary_AUC))

vein_AUC <- AUCell_calcAUC(veinSets, cells_rankings)
human_integrated@meta.data$vein_AUC_score <- t(getAUC(vein_AUC))

saveRDS(human_integrated, file = "human_integrated.rds")

# plot

fig_artery <- human_integrated@meta.data %>%
        plot_ly(
                x = ~newgroup,
                y = ~artery_AUC_score,
                split = ~newgroup,
                type = 'violin',
                box = list(
                        visible = T
                ),
                meanline = list(
                        visible = T
                )
        ) 

fig_capillary <- human_integrated@meta.data %>%
        plot_ly(
                x = ~newgroup,
                y = ~capillary_AUC_score,
                split = ~newgroup,
                type = 'violin',
                box = list(
                        visible = T
                ),
                meanline = list(
                        visible = T
                )
        ) 

fig_vein <- human_integrated@meta.data %>%
        plot_ly(
                x = ~newgroup,
                y = ~vein_AUC_score,
                split = ~newgroup,
                type = 'violin',
                box = list(
                        visible = T
                ),
                meanline = list(
                        visible = T
                )
        ) 

## integrated
human_integrated.df <- cbind(human_integrated@reductions$umap@cell.embeddings, human_integrated@meta.data)
#[1] "#F8766D" "#CD9600" "#7CAE00" "#00BE67" "#00BFC4" "#00A9FF" "#C77CFF"
#[8] "#FF61CC"
human_integrated.colour <- c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF","#FF61CC")

human_integrated.p <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ seurat_clusters, colors = human_integrated.colour, marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

uninjured <- plot_ly(data = human_integrated.df[human_integrated.df$newcondition == "Uninjured", ], x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ seurat_clusters, colors = human_integrated.colour, marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

injured <- plot_ly(data = human_integrated.df[human_integrated.df$newcondition == "Injured", ], x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ seurat_clusters, colors = human_integrated.colour, marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

human_integrated.artery <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(artery_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

human_integrated.vein <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(vein_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

human_integrated.capillary <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(capillary_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

# NOTE OF 14 July 2021
# At the moment the graph doesn't show the results in a best way. A dot plot might do a better job at disseminating the data.
## plot vascular beds signals 
DotPlot(human_integrated, features = c("capillary_AUC_score","vein_AUC_score", "artery_AUC_score"), scale.min = 0, scale.max = 100, col.min = -5, col.max = 5) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
## plot the vascular beds signals across all cells
hist(human_integrated$artery_AUC_score, breaks = 100)
hist(human_integrated$vein_AUC_score, breaks = 100)
hist(human_integrated$capillary_AUC_score, breaks = 100)
## Use AUCell to pull out vein, artery and capillary cells

#######################################################################
# Determine the cells with the given gene signatures or active gene sets
## artery
set.seed(123)
artery_cells_assignment <- AUCell_exploreThresholds(artery_AUC, plotHist = TRUE, assign = FALSE)
artery_cells_assignment$`artery ECs`$aucThr$thresholds # Global_k1 (0.29170126) is good

vein_cells_assignment <- AUCell_exploreThresholds(vein_AUC, plotHist = TRUE, assign = FALSE)
vein_cells_assignment$`vein ECs`$aucThr$thresholds # Global_k1 (0.2960563) is good

capillary_cells_assignment <- AUCell_exploreThresholds(capillary_AUC, plotHist = TRUE, assign = FALSE)
capillary_cells_assignment$`capillary ECs`$aucThr$thresholds # Global_k1 is good but a bit far away from the optimal, which is around 0.29 (0.2479438)

p1 <- DimPlot(human_integrated, cells.highlight = WhichCells(object = human_integrated, expression = artery_AUC_score >= 0.29170126)) + scale_color_manual(labels = c("non-artery EC", "artery EC"), values = c("slategray", "red")) + labs(color = "Cell type")

p2 <- DimPlot(human_integrated, cells.highlight = WhichCells(object = human_integrated, expression = vein_AUC_score >= 0.2960563)) + scale_color_manual(labels = c("non-vein EC", "vein EC"), values = c("slategray", "blue")) + labs(color = "Cell type")

p3 <- DimPlot(human_integrated, cells.highlight = WhichCells(object = human_integrated, expression = capillary_AUC_score >= 0.2479438)) + scale_color_manual(labels = c("non-capillary EC", "capillary EC"), values = c("slategray", "green")) + labs(color = "Cell type")

p1 + p2 + p3

## replot in plotly to improve the visualisation and thresholding can be done directly on the plot function
p4 <- human_integrated.artery <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(artery_AUC_score) >= 0.29170126, colors = c("gray", "red"), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

p5 <- human_integrated.vein <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(vein_AUC_score) >= 0.2960563, colors = c("gray", "blue"), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

p6 <- human_integrated.capillary <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ as.numeric(capillary_AUC_score) >= 0.2479438, colors = c("gray", "green"), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

p7 <- human_integrated.artery_vein <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ (as.numeric(artery_AUC_score) >= 0.29170126 & as.numeric(vein_AUC_score) >= 0.2960563), colors = c("gray", "magenta"), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

p8 <- human_integrated.capillary_vein <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ (as.numeric(capillary_AUC_score) >= 0.2479438 & as.numeric(vein_AUC_score) >= 0.2960563), colors = c("gray", "cyan"), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

p9 <- human_integrated.capillary_artery <- plot_ly(data = human_integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, type = "scatter", mode = "markers", color = ~ (as.numeric(capillary_AUC_score) >= 0.2479438 & as.numeric(artery_AUC_score) >= 0.29170126), colors = c("gray", "yellow"), marker = list(size = 8, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

