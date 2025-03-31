###### PEC 1. Analisis de datos omicos ###
### Alejandro Gombau García ###
### IMPORTACiÓN DE LOS DATOS ###
dataset_pec1 <- read.csv("https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv")


### CREACIÓN DEL OBJETO SUMMARIZED EXPERIMENTS ###
# Extraemos metadata de las muestras (colData)
col_data <- dataset_pec1[,c("Patient.ID", "Muscle.loss")]
rownames(col_data) <- dataset_pec1$Patient.ID # Usamos la id del paciente como indicador

# Extraemos la matriz de los metabolitos y la transponemos
metabo_matrix <- as.matrix((dataset_pec1[,-c(1,2)])) # columnas de la 3 a la 65
metabo_matrix <- t(metabo_matrix)

# Creamos metadatos para las filas usando los nombres de los metabolitos
row_data <- data.frame(Metabolito = rownames(metabo_matrix))
rownames(row_data) <- rownames(metabo_matrix)


# Creamos el objeto SummarizedExperiment
se_object <- SummarizedExperiment(
  assays = list(counts = metabo_matrix),
  rowData = row_data,
  colData = col_data
)


### BOXPLOT ###

boxplot(t(assay(se_object)), 
        main = "Distribución de valores por metabolito",
        ylab = "Concentración",
        las = 2,
        cex.axis = 0.7)

### HEATMAP ###

heatmap(cor(t(assay(se_object, "counts"))), 
        main = "Heatmap de correlación entre metabolitos")

### BOXPLOTS ###

# Seleccionamos 10 metabolitos a comparar
metabolitos <- c("X1.6.Anhydro.beta.D.glucose", "X1.Methylnicotinamide", 
                 "X2.Aminobutyrate", "X2.Hydroxyisobutyrate", "X2.Oxoglutarate", 
                 "X3.Aminoisobutyrate", "X3.Hydroxybutyrate", "X3.Hydroxyisovalerate", 
                 "X3.Indoxylsulfate", "X4.Hydroxyphenylacetate")

# Configuramos la disposición de la figura en 2 filas y 5 columnas
par(mfrow = c(2, 5))

# Iteramos sobre cada metabolito para crear un boxplot comparando los grupos
for(m in metabolitos) {
  boxplot(assay(se_object, "counts")[m, ] ~ colData(se_object)$Muscle.loss,
          main = paste("Comparación de", m),
          xlab = "Tipo de paciente",
          ylab = "Concentración",
          col = c("blue", "red"))
}

par(mfrow = c(1, 1))

### PCA ##
# Realizamos un PCA con la matriz de datos transpuesta
pca <- prcomp(t(assay(se_object, "counts")), scale. = TRUE)

# Definimos colores según el grupo (asumiendo que "cachexic" y "control" son los valores en Muscle.loss)
grupos <- colData(se_object)$Muscle.loss
colores <- ifelse(grupos == "cachexic", "red", "blue")

# Configuramos el layout para varios gráficos
par(mfrow = c(1, 3))  # 1 fila, 3 columnas

# PC1 vs PC2 (ya lo vimos, pero lo incluimos para referencia)
plot(pca$x[,1], pca$x[,2],
     main = "PC1 vs PC2",
     xlab = "PC1",
     ylab = "PC2",
     pch = 19, col = colores)
legend("topright", legend = unique(grupos), col = unique(colores), pch = 19)

# PC1 vs PC3
plot(pca$x[,1], pca$x[,3],
     main = "PC1 vs PC3",
     xlab = "PC1",
     ylab = "PC3",
     pch = 19, col = colores)
legend("topright", legend = unique(grupos), col = unique(colores), pch = 19)

# PC2 vs PC3
plot(pca$x[,2], pca$x[,3],
     main = "PC2 vs PC3",
     xlab = "PC2",
     ylab = "PC3",
     pch = 19, col = colores)
legend("topright", legend = unique(grupos), col = unique(colores), pch = 19)

par(mfrow = c(1,1))


### SCREEPLOT ###

screeplot(pca, type = "lines", main = "Scree Plot de PCA (metabolitos)")


### DENDROGRAMA ###

# Generar el dendrograma directamente a partir de la matriz transpuesta
dend <- as.dendrogram(hclust(dist(t(assay(se_object, "counts"))), method = "complete"))

# Asignar colores a las etiquetas según Muscle.loss: rojo para "cachexic", azul para "control"
labels_colors(dend) <- ifelse(
  colData(se_object)$Muscle.loss[match(labels(dend), colnames(assay(se_object, "counts")))] == "cachexic",
  "red", "blue"
)

# Graficar el dendrograma con las etiquetas coloreadas
plot(dend, main = "Dendrograma de Pacientes (color por grupo)", ylab = "Distancia")
legend("topright", legend = c("cachexic", "control"), fill = c("red", "blue"))

### VOLCANO PLOT ###
# Extraemos la información del grupo de cada paciente
grupo <- colData(se_object)$Muscle.loss

# Calculamos la media para cada metabolito en cada grupo
media_cachexic <- apply(assay(se_object, "counts")[, grupo == "cachexic"], 1, mean, na.rm = TRUE)
media_control  <- apply(assay(se_object, "counts")[, grupo == "control"], 1, mean, na.rm = TRUE)

# Calculamos el fold change: media en cachexic / media en control
fold_change <- media_cachexic / media_control

# Calculamos el log2 del fold change, lo cual facilita la interpretación (0 indica igualdad, valores positivos indican mayor concentración en cachexic)
log2_fold_change <- log2(fold_change)

# Creamos una tabla resumen con los resultados
resultados_fc <- data.frame(
  Metabolito    = rownames(assay(se_object, "counts")),
  Media_Cachexic = media_cachexic,
  Media_Control  = media_control,
  Fold_Change    = fold_change,
  Log2_FC        = log2_fold_change
)

# Calcular p-values para cada metabolito mediante un t-test
p_values <- apply(assay(se_object, "counts"), 1, function(x) {
  t.test(x ~ colData(se_object)$Muscle.loss)$p.value
})

# Crear un data frame para el volcano plot
volcano_data <- data.frame(
  Metabolito = rownames(assay(se_object, "counts")),
  Log2_FC = log2_fold_change,
  P_Value = p_values,
  NegLog10_P = -log10(p_values)
)

# Graficar el volcano plot
plot(volcano_data$Log2_FC, volcano_data$NegLog10_P,
     xlab = "Log2 Fold Change",
     ylab = "-Log10(p-value)",
     main = "Volcano Plot",
     pch = 20)
# Agregar una línea horizontal para p=0.05 (significación estadística)
abline(h = -log10(0.05), col = "red", lty = 2)

