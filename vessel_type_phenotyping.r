# Load count data
exprMatrix <- immuno@assays$RNA@data
exprMatrix <- as.matrix(exprMatrix) # convert sparse matrix to matrix
dim(exprMatrix)
# vessel type marker
marker <-read.csv("vessel_type_marker.csv")
new_marker <- marker %>% select(c(X, X.1, X.2, X.3))
new_marker <- new_marker[4:53, ]
colnames(new_marker) <- c("artery", "capillary", "vein", "lymphatic")
head(new_marker)

# Set gene set
artery_marker <- new_marker$artery
capillery_marker <- new_marker$capillary
vein_marker <- new_marker$vein

arterySets <- GeneSet(artery_marker, setName = "artery ECs")
capillerySets <- GeneSet(capillery_marker, setName = "capillery ECs")
veinSets <- GeneSet(vein_marker, setName = "vein ECs")

arterySets
capillerySets
veinSets

# Build gene expression rankings for each cell
detectCores()
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores = 2, plotStats = TRUE)
cells_rankings
save(cells_rankings, file = "EC_rankings_immuno.RData")

# Calculate enrichment for the gene signatures (AUC)
artery_AUC <- AUCell_calcAUC(arterySets, cells_rankings)
immuno@meta.data$artery_AUC_score <- t(getAUC(artery_AUC))

capillery_AUC <- AUCell_calcAUC(capillerySets, cells_rankings)
immuno@meta.data$capillery_AUC_score <- t(getAUC(capillery_AUC))

vein_AUC <- AUCell_calcAUC(veinSets, cells_rankings)
immuno@meta.data$vein_AUC_score <- t(getAUC(vein_AUC))

saveRDS(immuno, file = "immuno.rds")

# plot

fig <- immuno@meta.data %>%
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

fig <- immuno@meta.data %>%
        plot_ly(
                x = ~newgroup,
                y = ~capillery_AUC_score,
                split = ~newgroup,
                type = 'violin',
                box = list(
                        visible = T
                ),
                meanline = list(
                        visible = T
                )
        ) 

fig <- immuno@meta.data %>%
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