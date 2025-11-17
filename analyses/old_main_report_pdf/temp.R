

counts <- read.csv("../../../../data/crispr_counts.csv", row.names = 1)
crispr count data

metadata <- read.csv("../../../../data/crispr_metadata.csv", row.names = 1)

OsC <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata)
OsC@assays$RNA$data <- counts


library(Seurat)
library(dplyr)
library(biologicViewerSC)

all.genes <- rownames(OsC)


# Clean the data first
# Remove genes with all zeros
OsC <- subset(OsC, features = rownames(OsC)[Matrix::rowSums(OsC@assays$RNA$counts) > 0])

# Remove cells with very few features
#OsC <- subset(OsC, cells = colnames(OsC)[Matrix::colSums(OsC@assays$RNA$counts) > 10])


OsC <- OsC %>% 
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>% 
    RunPCA(npcs = 20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.2) %>%
    RunUMAP(dims = 1:20)


names(OsC@meta.data) <- gsub("\\.", "_", names(OsC@meta.data))



project_id <- "depmap_crispr"
projectPath = "/nemo/stp/babs/www/shiny/external/users/boeings/"

dbHostURL <- "clvd1-db-u-p-31.thecrick.org"
dbAdminUser <- "babs"

OsC <- JoinLayers(OsC)
params <- biologicViewerSC::scanObjParams(OsC)
#OsC@assays$RNA$data <- GetAssayData(OsC, layer = "data")
#OsC[['RNA']] <- as(OsC[['RNA']], "Assay")



###############################################################################
## Determine default gene                                                    ##

DefaultAssay(OsC) <- "RNA"
my_genes <- rownames(x = OsC@assays$RNA)


## Based on https://github.com/satijalab/seurat/issues/3560 the next two lines were 
# added/altered:
# in large datasets fetch data produces errors. After some trial and error, 
# breaking down the problem into chunks seems to solve the problem. 
cells <- row.names(OsC@meta.data)
cellList <- split(cells, ceiling(seq_along(cells)/10000))

## Define helper function ##
subset_fun <- function(obj, cellIDs, assay = "RNA") {
  z <- subset(obj, subset = cellID %in% cellIDs)@assays$RNA$data %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "gene") %>%
    tidyr::pivot_longer(
      !gene,
      names_to = "cellID",
      values_to = "lg10Expr"
    ) %>%
    filter(lg10Expr > 0)
    return(z)
}
## End of helper function

## Add cellID column ##
pos <- grep("cellID", names(OsC@meta.data))
if (length(pos) == 0){
    OsC@meta.data[["cellID"]] <- row.names(OsC@meta.data)
}

####################################################
## Create long-format list with expression values ##
dfExpr <- purrr::map(cellList, function(x) subset_fun(obj=OsC, cellIDs = x)) %>%
    dplyr::bind_rows() %>%
    data.frame()
## Done                                           ##
####################################################



biologicViewerSC::seuratObjectToViewer(
    params = params,
    project_id = project_id,
    projectPath = projectPath,
    OsC = OsC,
    dataMode = "MySQL",
    host = dbHostURL,
    dbname = "test_data",
    db.pwd = dbAdminPassword,
    db.user = "babs",
    appDomains = c("bioinformatics.crick.ac.uk","10.%"),
    geneDefault = "EZH2", # This will set the default gene in the app
    dfExpr = dfExpr
)