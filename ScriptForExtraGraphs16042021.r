# plot specific subpopulations across time
## dataframe1
P6_uninjured <- c(9.7, 12.6, 7.1, 13.7, 2.0, 7.7, 21.1, 6.0, 11.6, 18.8, 24.5, 11.3, 8.5, 4.8, 45.7)
P10_uninjured <- c(9.4, 5.8, 4.5, 12.0, 2.5, 6.7, 10.5, 1.4, 18.0, 17.1, 43.9, 24.5, 4.2, 7.5, 20.1)
P2_MI3W <- c(19.1, 2.3, NA, 7.0, NA, 11.3, 2.2, NA, 24.3, 10.2, 14.4, 10.7, 2.4, 14.5, 15.0)
Adult_uninjured <- c(9.4, 14.3, 27.5, 14.7, 16.9, 11.8, 5.6, 12.4, 12.8, 5.3, 5.9, 1.1, 22.3, 2.6, 14.0)
Adult_MI1D <- c(19.1, 8.9, NA, 4.6, NA, 3.0, 4.9, NA, NA, NA, 10.1, 20.2, NA, 39.8, NA)
Adult_MI3D <- c(14.5, 9.6, 11.2, 13.3, 13.0, 12.4, 8.7, 8.8, 10.7, 7.2, NA, 4.3, 12.1, 11.5, NA)
Adult_MI7D <- c(0.3, 19.0, 39.9, 10.4, 28.2, 10.3, 24.8, 43.5, 7.5, 16.3, 1.2, 0.1, 23.5, 3.6, 5.1)
Adult_MI2W <- c(13.0, 13.2, 5.8, 13.9, 11.1, 13.4, 10.8, 8.5, 7.9, 11.4, NA, 7.5, 12.7, 9.8, NA)
Adult_MI4W <- c(5.5, 14.4, 4.0, 10.6, 26.2, 23.3, 11.5, 19.4, 7.2, 13.8, 20.2, 14.4, 20.2, 14.4, NA)

percentage <- c(P6_uninjured, P10_uninjured, P2_MI3W, Adult_uninjured, Adult_MI1D, Adult_MI3D, Adult_MI7D, Adult_MI2W, Adult_MI4W)

group <- c(rep("P6_uninjured", 15), rep("P10_uninjured", 15), rep("P2_MI3W", 15), rep("Adult_uninjured", 15), rep("Adult_MI1D", 15), rep("Adult_MI3D", 15), rep("Adult_MI7D", 15), rep("Adult_MI2W", 15), rep("Adult_MI4W", 15))

cluster <- rep(as.character(0:14), 9)

data <- data.frame(group, cluster, percentage)
data$group <- factor(data$group, levels = c("P6_uninjured", "P10_uninjured", "P2_MI3W", "Adult_uninjured", "Adult_MI1D", "Adult_MI3D", "Adult_MI7D", "Adult_MI2W", "Adult_MI4W"))

fig <- plot_ly(data, x = ~ group, y = ~ percentage, color = ~ cluster, type = 'scatter', mode = 'markers') 

##################################################
p <- DimPlot(integrated, label = TRUE, label.size = 6) + NoLegend()
# isolate colour used in the Seurat package
colourbuilt <- ggplot_build(p)[[1]][[1]][,  c(1, 5)]
# reduce to unique colour and cluster information
further <- distinct(colourbuilt)
further <- further[order(further$group), ]
further.colour <- further$colour
##################################################
# styling plotly
ax <- list(
  zeroline = TRUE,
  showline = FALSE,
  showgrid = TRUE,
  tickangle = -45,
)

ay <- list(
  zeroline = TRUE,
  showline = FALSE,
  showgrid = TRUE
)

use.font <- list(
  family = "Open Sans",
  size = 16
  )
##################################################

## angiogenic clusters
use.colour1 <- further.colour[c(3, 6, 8, 10)] 
smalldata1 <- data %>% filter(cluster %in% c("2", "5", "7", "9"))
smallfig1 <- plot_ly(na.omit(smalldata1), x = ~ group, y = ~ percentage, color = ~ cluster, type = 'scatter', mode = 'lines+markers', marker = list(size = 10), colors = use.colour1) %>% layout(xaxis = ax, yaxis = ay, font = use.font) 

## inflammatory clusters
use.colour2 <- further.colour[c(1, 11, 13, 14)]
smalldata2 <- data %>% filter(cluster %in% c("0", "10", "12", "13"))
smallfig2 <- plot_ly(na.omit(smalldata2), x = ~ group, y = ~ percentage, color = ~ cluster, type = 'scatter', mode = 'lines+markers', marker = list(size = 10), colors = use.colour2) %>% layout(xaxis = ax, yaxis = ay, font = use.font) 

## script for sorting out the common up and down genes in injured vs healthy adult and P2 3WMI vs P6 healthy
'%!in%' <- function(x, y) {!'%in%'(x, y)}
## part of the script is from Module1.Rmd
# Figure 3 A
### first adult
AdultHealthyInjuredAllDEG <- read.csv("AdultHealthyInjuredAllDEG.csv")
#### up in injured
AdultHealthyInjuredAllDEG_up <- AdultHealthyInjuredAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
write.csv(AdultHealthyInjuredAllDEG_up, file = "AdultHealthyInjuredAllDEG_upInInjured.csv")
#### down in injured
AdultHealthyInjuredAllDEG_down <- AdultHealthyInjuredAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
write.csv(AdultHealthyInjuredAllDEG_down, file = "AdultHealthyInjuredAllDEG_downInInjured.csv")
### second neonatal
P6P2MI3WAllDEG <- read.csv("P6P2MI3WAllDEG.csv")
#### up in P2MI
P6P2MI3WAllDEG_up <- P6P2MI3WAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
write.csv(P6P2MI3WAllDEG_up, file = "P6P2MI3WAllDEG_upInP2MI3W.csv")
#### down in P2MI
P6P2MI3WAllDEG_down <- P6P2MI3WAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
write.csv(P6P2MI3WAllDEG_down, file = "P6P2MI3WAllDEG_downInP2MI3W.csv")
#### common up in injured
commonup <- intersect(AdultHealthyInjuredAllDEG_up$X, P6P2MI3WAllDEG_up$X)
#### common down in injured
commondown <- intersect(AdultHealthyInjuredAllDEG_down$X, P6P2MI3WAllDEG_down$X)
#### uniquely up in P2MI
uniqueup_P2MI <- P6P2MI3WAllDEG_up$X[P6P2MI3WAllDEG_up$X %!in% commonup]
#### uniquely down in P2MI
uniquedown_P2MI <- P6P2MI3WAllDEG_down$X[P6P2MI3WAllDEG_down$X %!in% commondown]
##### uniquely up in P2MI not in (Adult vs P6)
uniqueup_P2MI_notage <- uniqueup_P2MI[uniqueup_P2MI %!in% P6AdultAllDEG_up$X]
##### uniquely down in P2MI not in (Adult vs P6)
uniquedown_P2MI_notage <- uniquedown_P2MI[uniquedown_P2MI %!in% P6AdultAllDEG_down$X]
##### uniquely up in P2MI not in Adult and P10 together up
new_uniqueup_P2MI_notage <- uniqueup_P2MI[uniqueup_P2MI %!in% togetherup]
##### uniquely down in P2MI not in Adult and P10 together down
new_uniquedown_P2MI_notage <- uniquedown_P2MI[uniquedown_P2MI %!in% togetherdown]
## new branch
### P6 vs P10
P6P10AllDEG <- read.csv("P6P10AllDEG.csv")
#### up in P10
P6P10AllDEG_up <- P6P10AllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
write.csv(P6P10AllDEG_up, file = "P6P10AllDEG_upInP10.csv")
#### down in P10
P6P10AllDEG_down <- P6P10AllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
write.csv(P6P10AllDEG_down, file = "P6P10AllDEG_downInP10.csv")
### P6 vs Adult
P6AdultAllDEG <- read.csv("P6AdultAllDEG.csv")
#### up in Adult
P6AdultAllDEG_up <- P6AdultAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
write.csv(P6AdultAllDEG_up, file = "P6AdultAllDEG_upINadult.csv")
#### down in Adult
P6AdultAllDEG_down <- P6AdultAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
write.csv(P6AdultAllDEG_down, file = "P6AdultAllDEG_downINadult.csv")
### commonup
commonup <- intersect(P6P10AllDEG_up$X, P6AdultAllDEG_up$X)
### commondown
commondown <- intersect(P6P10AllDEG_down$X, P6AdultAllDEG_down$X)
#################
## need pie chart for the following two
#################
### together up
togetherup <- unique(c(P6P10AllDEG_up$X, P6AdultAllDEG_up$X))
### together down
togetherdown <- unique(c(P6P10AllDEG_down$X, P6AdultAllDEG_down$X))

## unique DEG for each time point in adult MI
### All DEG
#### MI1D
MI1DHealthyAllDEG <- read.csv("AdultMI1D_AdultHealthy_AllDEG.csv")
#### MI3D
MI3DHealthyAllDEG <- read.csv("AdultMI3D_AdultHealthy_AllDEG.csv")
#### MI7D
MI7DHealthyAllDEG <- read.csv("AdultMI7D_AdultHealthy_AllDEG.csv")
#### MI2W
MI2WHealthyAllDEG <- read.csv("AdultMI2W_AdultHealthy_AllDEG.csv")
#### MI4W
MI4WHealthyAllDEG <- read.csv("AdultMI4W_AdultHealthy_AllDEG.csv")
### All ups
MI1DHealthyAllDEG_MIup <- MI1DHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
MI3DHealthyAllDEG_MIup <- MI3DHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
MI7DHealthyAllDEG_MIup <- MI7DHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
MI2WHealthyAllDEG_MIup <- MI2WHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
MI4WHealthyAllDEG_MIup <- MI4WHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC < 0)
### All downs
MI1DHealthyAllDEG_MIdown <- MI1DHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
MI3DHealthyAllDEG_MIdown <- MI3DHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
MI7DHealthyAllDEG_MIdown <- MI7DHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
MI2WHealthyAllDEG_MIdown <- MI2WHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
MI4WHealthyAllDEG_MIdown <- MI4WHealthyAllDEG %>% dplyr::filter (p_val_adj < 0.05, avg_log2FC > 0)
### Common ups 48 in total
common_up <- Reduce(intersect, list(MI1DHealthyAllDEG_MIup$X, MI3DHealthyAllDEG_MIup$X, MI7DHealthyAllDEG_MIup$X, MI2WHealthyAllDEG_MIup$X, MI4WHealthyAllDEG_MIup$X))
### Common downs 76 in total
common_down <- Reduce(intersect, list(MI1DHealthyAllDEG_MIdown$X, MI3DHealthyAllDEG_MIdown$X, MI7DHealthyAllDEG_MIdown$X, MI2WHealthyAllDEG_MIdown$X, MI4WHealthyAllDEG_MIdown$X))
### unique ups
### unique up in MI1D
uniqueUp_MI1D <- MI1DHealthyAllDEG_MIup$X[MI1DHealthyAllDEG_MIup$X %!in% c(MI3DHealthyAllDEG_MIup$X, MI7DHealthyAllDEG_MIup$X, MI2WHealthyAllDEG_MIup$X, MI4WHealthyAllDEG_MIup$X, common_up)]
top20_uniqueUp_MI1D <- MI1DHealthyAllDEG_MIup %>% filter(X %in% uniqueUp_MI1D) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniqueUp_MI1D$X
### unique up in MI3D
uniqueUp_MI3D <- MI3DHealthyAllDEG_MIup$X[MI3DHealthyAllDEG_MIup$X %!in% c(MI1DHealthyAllDEG_MIup$X, MI7DHealthyAllDEG_MIup$X, MI2WHealthyAllDEG_MIup$X, MI4WHealthyAllDEG_MIup$X, common_up)]
top20_uniqueUp_MI3D <- MI3DHealthyAllDEG_MIup %>% filter(X %in% uniqueUp_MI3D) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniqueUp_MI3D$X
### unique up in MI7D
uniqueUp_MI7D <- MI7DHealthyAllDEG_MIup$X[MI7DHealthyAllDEG_MIup$X %!in% c(MI1DHealthyAllDEG_MIup$X, MI3DHealthyAllDEG_MIup$X, MI2WHealthyAllDEG_MIup$X, MI4WHealthyAllDEG_MIup$X, common_up)]
top20_uniqueUp_MI7D <- MI7DHealthyAllDEG_MIup %>% filter(X %in% uniqueUp_MI7D) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniqueUp_MI7D$X
### unique up in MI2W
uniqueUp_MI2W <- MI2WHealthyAllDEG_MIup$X[MI2WHealthyAllDEG_MIup$X %!in% c(MI1DHealthyAllDEG_MIup$X, MI3DHealthyAllDEG_MIup$X, MI7DHealthyAllDEG_MIup$X, MI4WHealthyAllDEG_MIup$X, common_up)]
top20_uniqueUp_MI2W <- MI2WHealthyAllDEG_MIup %>% filter(X %in% uniqueUp_MI2W) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniqueUp_MI2W$X
### unique up in MI4W
uniqueUp_MI4W <- MI4WHealthyAllDEG_MIup$X[MI4WHealthyAllDEG_MIup$X %!in% c(MI1DHealthyAllDEG_MIup$X, MI3DHealthyAllDEG_MIup$X, MI7DHealthyAllDEG_MIup$X, MI2WHealthyAllDEG_MIup$X, common_up)]
top20_uniqueUp_MI4W <- MI4WHealthyAllDEG_MIup %>% filter(X %in% uniqueUp_MI4W) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniqueUp_MI4W$X
### unique downs
### unique down in MI1D
uniquedown_MI1D <- MI1DHealthyAllDEG_MIdown$X[MI1DHealthyAllDEG_MIdown$X %!in% c(MI3DHealthyAllDEG_MIdown$X, MI7DHealthyAllDEG_MIdown$X, MI2WHealthyAllDEG_MIdown$X, MI4WHealthyAllDEG_MIdown$X, common_down)]
top20_uniquedown_MI1D <- MI1DHealthyAllDEG_MIdown %>% filter(X %in% uniquedown_MI1D) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniquedown_MI1D$X
### unique down in MI3D
uniquedown_MI3D <- MI3DHealthyAllDEG_MIdown$X[MI3DHealthyAllDEG_MIdown$X %!in% c(MI1DHealthyAllDEG_MIdown$X, MI7DHealthyAllDEG_MIdown$X, MI2WHealthyAllDEG_MIdown$X, MI4WHealthyAllDEG_MIdown$X, common_down)]
top20_uniquedown_MI3D <- MI3DHealthyAllDEG_MIdown %>% filter(X %in% uniquedown_MI3D) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniquedown_MI3D$X
### unique down in MI7D
uniquedown_MI7D <- MI7DHealthyAllDEG_MIdown$X[MI7DHealthyAllDEG_MIdown$X %!in% c(MI1DHealthyAllDEG_MIdown$X, MI3DHealthyAllDEG_MIdown$X, MI2WHealthyAllDEG_MIdown$X, MI4WHealthyAllDEG_MIdown$X, common_down)]
top20_uniquedown_MI7D <- MI7DHealthyAllDEG_MIdown %>% filter(X %in% uniquedown_MI7D) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniquedown_MI7D$X
### unique down in MI2W
uniquedown_MI2W <- MI2WHealthyAllDEG_MIdown$X[MI2WHealthyAllDEG_MIdown$X %!in% c(MI1DHealthyAllDEG_MIdown$X, MI3DHealthyAllDEG_MIdown$X, MI7DHealthyAllDEG_MIdown$X, MI4WHealthyAllDEG_MIdown$X, common_down)]
top20_uniquedown_MI2W <- MI2WHealthyAllDEG_MIdown %>% filter(X %in% uniquedown_MI2W) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniquedown_MI2W$X
### unique down in MI4W
uniquedown_MI4W <- MI4WHealthyAllDEG_MIdown$X[MI4WHealthyAllDEG_MIdown$X %!in% c(MI1DHealthyAllDEG_MIdown$X, MI3DHealthyAllDEG_MIdown$X, MI7DHealthyAllDEG_MIdown$X, MI2WHealthyAllDEG_MIdown$X, common_down)]
top20_uniquedown_MI4W <- MI4WHealthyAllDEG_MIdown %>% filter(X %in% uniquedown_MI4W) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20_uniquedown_MI4W$X
### saves
#### saves ups
write.csv(uniqueUp_MI1D, file = "uniqueUp_MI1D.csv")
write.csv(uniqueUp_MI3D, file = "uniqueUp_MI3D.csv")
write.csv(uniqueUp_MI7D, file = "uniqueUp_MI7D.csv")
write.csv(uniqueUp_MI2W, file = "uniqueUp_MI2W.csv")
write.csv(uniqueUp_MI4W, file = "uniqueUp_MI4W.csv")
#### saves downs
write.csv(uniquedown_MI1D, file = "uniqueDown_MI1D.csv")
write.csv(uniquedown_MI3D, file = "uniqueDown_MI3D.csv")
write.csv(uniquedown_MI7D, file = "uniqueDown_MI7D.csv")
write.csv(uniquedown_MI2W, file = "uniqueDown_MI2W.csv")
write.csv(uniquedown_MI4W, file = "uniqueDown_MI4W.csv")

## isolate cluster 2, 5, 7, 9 as the angiogenic cells to plot the genes 
Idents(integrated) <- "seurat_clusters"
levels(integrated)
angio <- subset(integrated, idents = c("2", "5", "7", "9"))
DefaultAssay(angio)
angio <- NormalizeData(angio)
DefaultAssay(angio) <- "RNA"
Idents(angio) <- "group"
levels(angio)
saveRDS(angio, file = "angio.rds")
