rm(list=ls())
library(dplyr)
mouse_cells_meta <- read.csv("experiment/tabula_muris/annotations_facs.csv")
mouse_facs_meta <- read.csv("experiment/tabula_muris/metadata_FACS.csv")

celltype_counts <- mouse_cells_meta %>% group_by(tissue, plate.barcode, cell_ontology_class) %>%
  summarise(count=n())

spleen <- celltype_counts %>% filter(tissue == "Spleen")
fat <- celltype_counts %>% filter(tissue == "Fat")


fat_cells <- mouse_cells_meta %>% filter(tissue == "Fat") %>%
  select(cell, cell_ontology_class, mouse.id, mouse.sex, plate.barcode, tissue)

fat_cells_majority <- fat_cells %>% filter(cell_ontology_class %in% c("B cell", "endothelial cell", "mesenchymal stem cell of adipose",
                                                                      "myeloid cell", "T cell"))

set.seed(2024)
fat_cells_sample <- fat_cells_majority %>% group_by(cell_ontology_class) %>%
  slice_sample(n=100)
rownames(fat_cells_sample) <- gsub("_", "-", fat_cells_sample$cell)


fat_counts <- read.csv("experiment/tabula_muris/FACS/Fat-counts.csv")
rownames(fat_counts) <- gsub("_", "-", fat_counts$X)
fat_counts$X <- NULL
colnames(fat_counts) <- gsub("_", "-", colnames(fat_counts))

library(Matrix)
# change to sparse matrix
fat_counts <- as(as.matrix(fat_counts), "sparseMatrix")
fat_cells_sample$cell <- gsub("_", "-", fat_cells_sample$cell)
fat_counts_sample <- fat_counts[, fat_cells_sample$cell]

# go through Seurat pipeline to find the most variable genes
library(Seurat)
fat_seurat <- CreateSeuratObject(counts = fat_counts_sample,
                                 assay = "RNA", project = "mouse_fat",
                                 meta.data=fat_cells_sample)

VlnPlot(fat_seurat, features = c("nFeature_RNA", "nCount_RNA"))

# find the most variable genes
fat_seurat <- NormalizeData(fat_seurat, normalization.method = "LogNormalize", scale.factor = 100000)
fat_seurat <- FindVariableFeatures(fat_seurat, selection.method = "vst", nfeatures = 1000)
head(VariableFeatures(fat_seurat), 10)
variable_counts <- fat_seurat[["RNA"]]$counts[VariableFeatures(fat_seurat), ]

writeMM(variable_counts, "experiment/tabula_muris/ZINBwave/mouse_fat_template_mat.mtx")
writeLines(rownames(variable_counts), "experiment/tabula_muris/ZINBwave/genes_template.txt")
writeLines(colnames(variable_counts), "experiment/tabula_muris/ZINBwave/cells_template.txt")
write.table(fat_seurat[[]], "experiment/tabula_muris/ZINBwave/cells_metadata.tsv", sep='\t',
          row.names=F, quote=F)


counts_dense <- as.matrix(variable_counts)
log_counts_dense <- log10(1+counts_dense)  



# visualize the samples, confirm that the expression profiles are distinct

library(ggplot2)
library(reshape2)

plot_heatmap <- function(input_matrix, legend_title, has.legend=T, 
                         colormap="viridis", min=0, max=3){
  
  rownames(input_matrix) <- NULL
  colnames(input_matrix) <- NULL
  
  long_matrix <- melt(input_matrix)
  base_plot <- ggplot(data = long_matrix, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() + labs(fill=legend_title) + 
    scale_fill_viridis_c(option = "viridis", limits=c(min, max))+
    theme_minimal()
  
  if (has.legend){
    final_plot <- base_plot + theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.title.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks.y=element_blank(),
                                    axis.title.y=element_blank())
  } else{
    final_plot <- base_plot + theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.title.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks.y=element_blank(),
                                    axis.title.y=element_blank(),
                                    legend.position = "none")
  }
  
  return(final_plot)
}

my_heatmap <- plot_heatmap(log_counts_dense, legend_title="Log10Count",
                           min=0, max=3)


