library(data.table)
df_1kb <- fread("~/Downloads/KLF4.1.tsv")
head(df_1kb)
dim(df_1kb)
list_targets <- df_1kb$Target_genes
list_check1 <- c("VEGFC", "FLT4", "NRP2", "GAS6", "TSC22D3") # "NRP2" (22.862069)    "GAS6" (9.206897)    "TSC22D3" (4.275862)
intersect1 <- intersect(list_check1, list_targets)
df_1kb %>% filter(Target_genes %in% intersect1) %>% select (`KLF4|Average`)

# top 10 upregulated in cHF EC vs uninjured
list_check2 <- c("MSRB", "MYH6", "MTATP8P1", "AP000251.3", "rP11-777B9.5", "CH507-513H4.5", "MTRNR2L12", "MTND6P3", "MTRNR2L8", "RP11-100K18.1")
intersect2 <- intersect(list_check2, list_targets) # "MTRNR2L8" (176.2414)
df_1kb %>% filter(Target_genes %in% intersect2) %>% select (`KLF4|Average`)

# top 10 upregulated in injured mouse vs uninjured
## within regenerative window
list_check3 <- c("Mylip", "Gm17660", "Anks1b", "Asxl3", "Chl1", "Dlgap1", "Nkain2", "Dgkb", "Rps27rt", "Gm26694")
new_list_check3 <- mouse2human(list_check3)$humanGene
intersect3 <- intersect(new_list_check3, list_targets) # "Mylip" (10.14), "Anks1b" (7.55), "Nkain2" (6.10)
df_1kb %>% filter(Target_genes %in% intersect3) %>% select (`KLF4|Average`)

## without the window
list_check4 <- c("Cct6a", "Cops9", "Elob", "Eloc", "Gas5", "Ndfip1", "Nme2", "Rack1", "Rflnb", "Selenof")
new_list_check4 <- mouse2human(list_check4)$humanGene
intersect4 <- intersect(new_list_check4, list_targets) # "Cct6a" (37.52), "Ndfip1" (18.48), "Nme2" (10.10)
df_1kb %>% filter(Target_genes %in% intersect4) %>% select (`KLF4|Average`)

## commonly upregulated between mouse and human
list_check5 <- c("Sox18", "Rbpms", "Hspa8", "Rps28", "Gm26917")
new_list_check5 <- mouse2human(list_check5)$humanGene
intersect5 <- intersect(new_list_check5, list_targets) # "Hspa8" (21.93), "Rps28" (4.34)
df_1kb %>% filter(Target_genes %in% intersect5) %>% select (`KLF4|Average`)




