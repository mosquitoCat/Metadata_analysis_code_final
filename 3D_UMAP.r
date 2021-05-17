library(plotly)

set.seed(1024)
integrated <- RunUMAP(integrated,
                      dims = 1:50,
                      n.components = 3L)

# Extract tSNE information from Seurat Object
umap_1 <- integrated[["umap"]]@cell.embeddings[,1]
umap_2 <- integrated[["umap"]]@cell.embeddings[,2]
umap_3 <- integrated[["umap"]]@cell.embeddings[,3]

# Create a new dataframe
integrated.df <- FetchData(integrated, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "newgroup","seurat_clusters", "tip_AUC_score", "stalk_AUC_score", "artery_AUC_score", "vein_AUC_score", "capillary_AUC_score"))

integrated.colour <- c("#F8766D", "#E58700", "#C99800", "#A3A500", "#6BB100", "#00BA38", "#00BF7D", "#00C0AF", "#00BCD8", "#00B0F6", "#619CFF", "#B983FF", "#E76BF3", "#FD61D1", "#FF67A4")

integrated.p <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, z = ~ UMAP_3, type = "scatter3d", mode = "markers", color = ~ seurat_clusters, colors = integrated.colour, marker = list(size = 2, line = list(color = "black", width = 0.5)), width = 1000, height = 1000) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.tip <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, z = ~ UMAP_3, type = "scatter3d", mode = "markers", color = ~ as.numeric(tip_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 1000, height = 1000) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.stalk <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, z = ~ UMAP_3, type = "scatter3d", mode = "markers", color = ~ as.numeric(stalk_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.artery <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, z = ~ UMAP_3, type = "scatter3d", mode = "markers", color = ~ as.numeric(artery_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.vein <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, z = ~ UMAP_3, type = "scatter3d", mode = "markers", color = ~ as.numeric(vein_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()

integrated.capillary <- plot_ly(data = integrated.df, x = ~ UMAP_1, y = ~ UMAP_2, z = ~ UMAP_3, type = "scatter3d", mode = "markers", color = ~ as.numeric(capillary_AUC_score), colors = viridisLite::viridis(20), marker = list(size = 6, line = list(color = "black", width = 0.5)), width = 500, height = 500) %>% layout(title = "", yaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), xaxis = list(zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, title = ""), font = t, legend = list(x = 0.9, y = 0.1, font = list(size = 30))) %>% hide_colorbar() %>% hide_legend() %>% hide_guides()
