library(Seurat)
library(harmony)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)

### Pipeline for every clustering iteration of PBMC objects, sorted CD8 T cells and overall CSF object ###

data <- FindVariableFeatures(object = data, nfeatures = 1500)
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB", data@assays$RNA@var.features)] # ^IGH|^IGK|^IGL for B cells
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch", "Donor"), assay.use = "RNA", max.iter.harmony = 30)
data <- RunUMAP(data, dims = 1:15, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:15, reduction = "harmony")
data <- FindClusters(data, resolution = 0.5)

### Pipeline for CSF CD8 T, CD4 T and myeloid cells ###

data.list <- SplitObject(data, split.by = "Donor")
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  x@assays$RNA@var.features <- x@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", x@assays$RNA@var.features)]
  x
})
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 2000)
for(i in names(data.list)){
  data.list[[i]] <- ScaleData(data.list[[i]], features = features, verbose = FALSE)
  data.list[[i]] <- RunPCA(data.list[[i]], features = features, npcs = 30)
}
immune.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features, dims = 1:30, reduction = 'rpca')
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, vars.to.regress = c("nCount_RNA", "percent.mt"))
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
