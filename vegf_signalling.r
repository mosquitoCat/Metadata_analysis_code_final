# Check Vegf signalling pathways

# for combined datasets
features_to_plot <- c("Vegfa", "Vegfb", "Vegfc", "Figf", "Kdr", "Flt1", "Flt4", "Nrp1", "Nrp2") # Vegfd is Figf

features_to_plot <- c("Vegfa", "Vegfb", "Vegfc", "Figf")

features_to_plot <- c("Kdr", "Flt1", "Flt4", "Nrp1", "Nrp2") 

Idents(combined) <- "newgroup"

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

DotPlot(combined, features = features_to_plot, group.by = "newgroup", cols = c("navy", "plum1")) + theme_bw()

FeaturePlot(combined, features = features_to_plot, cols = c("navy", "plum1"))

# for mouse datasets only

features_to_plot <- c("Vegfa", "Vegfb", "Vegfc", "Figf", "Kdr", "Flt1", "Flt4", "Nrp1", "Nrp2") # Vegfd is Figf

features_to_plot <- c("Vegfa", "Vegfb", "Vegfc", "Figf")

features_to_plot <- c("Kdr", "Flt1", "Flt4", "Nrp1", "Nrp2") 

Idents(mouse_integrated) <- "paperID"

DefaultAssay(mouse_integrated) <- "RNA"
mouse_integrated <- NormalizeData(mouse_integrated)
mouse_integrated$paperID <- factor(mouse_integrated$paperID, levels = c("uninjured P6", "uninjured P10", "P2 MI21D", "uninjured adult", "Adult MI1D", "Adult MI3D", "Adult MI7D", "Adult MI14D", "Adult MI28D"))

DotPlot(mouse_integrated, features = features_to_plot, group.by = "paperID", cols = c("navy", "plum1")) + theme_bw()

FeaturePlot(mouse_integrated, features = features_to_plot, cols = c("navy", "plum1"))

# for human datasets only

features_to_plot <- c("VEGFA", "VEGFB", "VEGFC", "FIGF", "KDR", "FLT1", "FLT4", "NRP1", "NRP2") # Vegfd is Figf

features_to_plot <- c("VEGFA", "VEGFB", "VEGFC", "FIGF")

features_to_plot <- c("KDR", "FLT1", "FLT4", "NRP1", "NRP2") 

Idents(human_integrated) <- "newgroup"

DefaultAssay(human_integrated) <- "RNA"
human_integrated <- NormalizeData(human_integrated)


DotPlot(human_integrated, features = features_to_plot, group.by = "newgroup", cols = c("navy", "plum1")) + theme_bw()

FeaturePlot(human_integrated, features = features_to_plot, cols = c("navy", "plum1"))
