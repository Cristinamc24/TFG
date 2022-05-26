# Este script descarga los datos de expresion del dataset GSE106977 y
# obtiene:
# - matriz de expresion
# - tabla con la variable clinica (pCR) Paso previo a la construccion y
#   analisis de la red de coexpresion

########################################################################

# 0. Libraries and setup
library(GEOquery) library(limma) l ibrary(umap)

probeid2genename = function(probeid){ require(dplyr) df =
gset@featureData@data fila = (dplyr::filter(df , ID ==
probeid))$gene_assignment if (length(fila) == 0){return(NA)} else { genename
= strsplit(fila, split = "//")[[1]][2] return(gsub(" ","",genename)) } }

Sys.setenv("VROOM_CONNECTION_SIZE" = 50000000) # necesario para descargar
datos con getGEO

# 1. Data download
gset <- getGEO("GSE106977", GSEMatrix =TRUE, AnnotGPL=FALSE)

if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <-
1 gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- paste0("00010110000001000000100001010110001100010110010000",
               "01111101110000001100000001010010010100100001011010",
               "1100000110011110111") sml <- strsplit(gsms, split="")[[1]]

# 2. Get the expression matrix
ex <- exprs(gset)
# Its row names are the probe ID, which we need to transforn to Gene names or
# Gene IDs gset@featureData@data is a dataframe with probes annotation. The
# column gene_assignment has the relevant info

rownames(ex) <- sapply(rownames(ex), probeid2genename) # annotate rownames with gene names

ex <- na.omit(ex) # eliminate rows with NAs 
ex <- ex[!duplicated(ex), ]  #remove duplicates 
ex <- ex[2722:68540,] # remove rows with NA in the row name

# 3. Build phenotype dataframe
gs <- factor(sml)
groups <- make.names(c("no pcr","pcr"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

pcrTable <- data.frame(muestras = gset$geo_accession,ValorPCR = as.numeric(sml) , check.names = TRUE)
pcrTable[pcrTable$ValorPCR == 1,]$ValorPCR <- 2
pcrTable[pcrTable$ValorPCR == 0,]$ValorPCR <- 1
pcrTable[pcrTable$ValorPCR == 2,]$ValorPCR <- 0


# 4. save output
write.csv(pcrTable, file = "tabla_pCR.csv")

save(ex, file = "expression_matrix.Rdata")
save(pcrTable, file = "tabla_pCR.Rdata")

