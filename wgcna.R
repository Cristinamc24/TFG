# Este script modela la red de coexpresión desde la matriz de expresión del dataset GSE106977 y obtiene los módulos de genes asociados con pCR

########################################################################

# 0. Libraries and setup
library(WGCNA) 
library(dplyr) 


### 1. Load the expression matrix
load("expression_matrix.Rdata")

# comprueba que la matriz se ha cargado bien
colnames(ex) # deben salir las muestras, 119 muestras
rownames(ex) # deben salir los genes
dim(ex) #65819 genes y 119 muestras

ex2 = ex[!is.na(rownames(ex)),] # eliminamos los genes que no tienen nombre asignado (el nombre es NA)
ex3 = ex2[!(rownames(ex2) == "---"),] # eliminamos los genes que no estan bien identificados

## transformamos la matriz de expresion para WGCNA

datExpr0 = as.data.frame(t(ex3)) #hacemos la traspuesta los genes en columnas y las muestras en filas. 


gsg = goodSamplesGenes(datExpr0, verbose = 3) # miramos si hay genes que tienen pocos datos de expresion a traves de las muestras
gsg$allOK #si sale TRUE es que todos los genes tienen suficientes valores de expresion.

# En caso de que haya genes o muestras que le faltan valores de expresion, se identifican y se eliminan
if (!gsg$allOK)
{
        # Optionally, print the gene and sample names that were removed:
        if (sum(!gsg$goodGenes)>0) 
                printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0) 
                printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
        # Remove the offending genes and samples from the data:
        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


## Clustering de las muestras

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 95, col = "red");
dev.off()

# Determine cluster under the line
# eliminamos la muestra outlier
clust = cutreeStatic(sampleTree, cutHeight = 95, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)



### 2. Load clinical data

traitData = read.csv("tabla_pcr.csv");
allTraits = traitData[, -c(1)]

dim(allTraits)
names(allTraits)


# Form a data frame analogous to expression data that will hold the clinical traits.

Samples = rownames(datExpr);
traitRows = match(Samples, allTraits$muestras);
datTraits = as.data.frame(allTraits[traitRows, 2])
# names(datTraits) <- names(allTraits)[2]
names(datTraits) <- "pCR"
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();



# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(file="sample_clustering-pCR.pdf", width=20, height =9)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()


## guardamos los datos procesados

save(datExpr, datTraits, file = "WCGNA_dataInput.RData")




### 3. Run WGCNA

# Si empezamos por aqui no hace falta volver a procesar los datos
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
library(dplyr)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the data saved in the first part
lnames = load(file = "WCGNA_dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames


# Allow multi-threading within WGCNA. This helps speed up certain calculations.
allowWGCNAThreads()

# Ojo, si se usa RStudio, allowWGCNAThreads da error, hay que usar
disableWGCNAThreads()




## Choose a set of soft-thresholding powers


powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf(file = "soft_threshold.pdf", width=9, height=5)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()



# vemos que podemos elegir 10 o 12. A continuacion te muestro el código para obtener los resultados con 12. Puedes obtenerlos tambien con 10


soft_power = 12

## Build the co-expression network

net = blockwiseModules(datExpr, power = soft_power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)


pdf(file = "coexp_modules.pdf", width=20, height=9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "networkConstruction-auto.RData") # guardamos la red y la info de modulos

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);




pdf(file = "module_pCR-correlation.pdf", width = 4, height = 6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


# los modulos green (correlacion -0.24 p-val 0.008) y black (correlacion 0.18 pval 0.05) muestran correlaciones significativas. Cuando uses otro valor de soft_power te pueden salir otros modulos con correlaciones significativas.
# La asignacion de colores a los modulos es arbitraria, el color no significa nada en si, solo sirve para identificar el modulo. Es posible que cuando corras el codigo, salgan otros colores. Lo importante es que con el soft_power de 12 hay dos modulos interesantes en los que la expresion de los genes esta corelacionada con pCR (-0.24 y 0.18) da igual que color tengan los modulos. 




# Define variable pCR containing the pCR column of datTrait
pCR = as.data.frame(datTraits$pCR);
names(pCR) = "pCR"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, pCR, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(pCR), sep="");
names(GSPvalue) = paste("p.GS.", names(pCR), sep="");

# En este bloque hemos definido dos tablas de una columna cada una, geneTraitSignificance tiene la correlación entre el perfil de expresion de cada gen de forma individual (no por el modulo al que pertenece) y pCR. GSPvalue tiene el p-valor de esa correlacion.


## exploracion de cada modulo
# uso el modulo green, debes sacar un plot como este para cada modulo que te salga interesante en module_pCR-correlation.pdf

module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf(file=paste("modulo_", module, ".pdf" sep=""), width = 7, height =7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body pCR",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()


# genes que pertenecen al modulo que estamos estudiando
names(datExpr)[moduleColors==module] #154 genes en el modulo green. 62 genes en el modulo black


## dataframe con la informacion de modulos y correlacion

# Create the starting data frame
genenames =  names(datExpr)
geneInfo0 = data.frame(gene.name = genenames,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for pCR
modOrder = order(-abs(cor(MEs, pCR, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pCR));
geneInfo = geneInfo0[geneOrder, ]


write.csv(geneInfo, file = "geneInfo.csv") #guardamos la informacion general de genes y todos los modulos



## Sacamos listas de genes en los modulos interesantes y los datos de correlacion de estos genes. Ojo, la variable intModules tiene que ajustarse a tus resultados

# $ Choose interesting modules
intModules = c("green", "black")
for (module in intModules)
{
        # Select module genes
        modGenes = (moduleColors==module)
        # Get their gene name
        modLLIDs = genenames[modGenes];
        # Write them into a file
        fileName = paste("genes-", module, ".txt", sep="")
        fileName2 = paste("data-", module, ".csv", sep="");
        write.table(as.data.frame(modLLIDs), file = fileName,
                    row.names = FALSE, col.names = FALSE, quote=F)
        # get their correlation with pCR and membership to thier module
        mod_data = dplyr::filter(geneInfo, gene.name %in% modLLIDs) %>% dplyr::select(gene.name, GS.pCR, p.GS.pCR, paste("MM.",module, sep=""), paste("p.MM.",module, sep=""))
        write.csv(mod_data, file = fileName2, row.names = FALSE, quote = F)
}




