####################
#                  #
# Copiar todo esto #
#                  #
####################

# Hecho con gusto por Carla Carolina P�rez Hern�ndez (UAEH)

# Modificado en pr�ctica de Laboratorio 36 por Araceli Guerrero Herrera (ZuRy) :3

# Laboratorio 36 - MAPA DE CALOR -T�RMICO- with pheatmap
# DATOS GEN�TICOS TOMADOS DE Sahir Bhatnagar
# PR�CTICA DE CODERS

# Objetivo: Realizar un heatmap con datos gen�ticos 
# ------------------------------------------------------------------------------------------------------
# En este ejercicio vamos a:
# 1. Cargar nuestra matriz hipot�tica de datos y dataframes adicionales
# 2. Realizar varios heatmaps

#Un mapa de calor es una representaci�n gr�fica de datos que utiliza un sistema de codificaci�n de colores para representar diferentes valores


#Heatmaps with pheatmap 
#Simulated data created by Sahir Bhatnagar.

#possible data pre-processing - normalization - quantile, median, etc., log transform
#not necessary here - we have log fold change data that has already been normalized

#Calculating your distance matrix (see dist objects):
#compute how similar or different you values are
#parametric - distance measures based on Pearson correlation 
#non parametric - spearman rank - replace by ranks and calculate correlation, Kendall's - relative ordering
#euclidean - shortest distance between values (has to be normalized), takes magnitude into account
#city block/Manhattan - sum of distances along each dimension
#distance 1-correlation - of all pairs of items to be clustered

#Cluster your samples (see hclust objects):
#hierarchical, organizes into a tree structure based on similarity - short branches if similar and longer branches as similarity decreases
#repeated cycles where the 2 closest remaining items (smallest distance) get joined by a branch with the length of the branch reflecting the distance between them, the distance between this item and all other remaining items are computed until only one object remains
#single linkage clustering - distance between 2 items is the minimum of all pairwise distances between items contained in x and y - fast b/c no other calculations need to be performed once you have your distance matrix
#complete linkage is the maximum of all paiwise distances between x and y 
#average linkage - mean of all pairwise distances between items contained in x and y
#k-means organize into clusters (self-chosen number) - items are randomly assigned to a cluster - the mean vector fo rall items in each hcluster is computed, items are reassigned to the cluster whose center is closest to them - random starting points so will not always get the same answer, number of trial done to deal with the randomness
#self organizing maps


# Paso 1. Instalar y ejecutar la paqueter�a de pheatmap.
install.packages("pheatmap")
library(pheatmap)

# Paso 2. Importar datos (ubicar la ruta de los 3 archivos que se utilizar�n).
file.choose()

# Paso 3. Copiar y pegar la ubicaci�n del primer archivo para abrirlo como matriz.
genes <- as.matrix(
  read.csv("C:\\Users\\Araceli Guerrero\\Documents\\Laboratorios de R\\LAB36\\heatmap_data.csv",
           sep = ",",
           header = T,
           row.names = 1))

# Paso 4. Visualizar los primeros datos de la matriz genes.
head(genes[,1:5])

# Paso 5. Copiar y pegar la ubicaci�n del archivo de las anotaciones de columnas para leerlo como dataframe.
annotation_col <- read.csv("C:\\Users\\Araceli Guerrero\\Documents\\Laboratorios de R\\LAB36\\annotation_col.csv",
           header = T,
           row.names = 1)

# Paso 6. Visualizar los primeros datos del dataframe annotation_col.
head(annotation_col[,1:2])

# Paso 7. Copiar y pegar la ubicaci�n del archivo de las anotaciones de renglones para leerlo como dataframe.
annotation_row <- read.csv("C:\\Users\\Araceli Guerrero\\Documents\\Laboratorios de R\\LAB36\\annotation_row.csv",
                           header = T,
                           row.names = 1)

# Paso 8. Visualizar los primeros datos del dataframe annotation_row.
head(annotation_row[,1:1])

# Paso 9. Dibujar el heatmap.
pheatmap(genes)

# Paso 10. Cambiar el tama�o de la fuente (font).
pheatmap(genes, fontsize = 6)

# Por default se dibujan los dendogramas (clustering) tanto en los renglones como en las columnas.

# Paso 11. Cl�ster por genes (grupos de genes similares) se encuentra en los renglones. 
# Eliminar dendograma de las columnas (pacientes).
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = F)

# Paso 12. Cl�ster por pacientes (grupos de pacientes similares) se encuentra en las columnas.
# Eliminar dendograma de los renglones (genes).
pheatmap(genes, fontsize = 6, cluster_rows = F, cluster_cols = T)

# Paso 13. Visualizar ambos dendogramas.
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T)

# Paso 14. Identificar patrones subyacentes a las anotaciones de los renglones.
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row)

# Paso 15. Agregar las anotaciones de las columnas.
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col)

# Los genes est�n en cl�ster de acuerdo a los patrones subyacentes a los que pertenecen.
# Ahora se tiene informaci�n sobre el medicamento y la condici�n.
 
# Paso 16. Generar el gr�fico completo quitando cl�sters (�rboles de agrupaci�n o dendogramas).
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col,
         treeheight_row = 0, treeheight_col = 0, main = "Expresi�n Gen�tica")

# Paso 17. Tomar una muestra de la matriz original para crear un subset.
sub <- genes [c(1:5, 55:60), c(1:5, 20:35, 55:60)] 

# Paso 18. Visualizar los primeros datos del subset.
head(sub[,1:5])

# Paso 19. Generar el gr�fico del subset.
pheatmap(sub, fontsize = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col,
         treeheight_row = 0, treeheight_col = 0, main = "Expresi�n Gen�tica")

# Paso 20. Desplegar los valores del gr�fico.
pheatmap(sub, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col,
         treeheight_row = 0, treeheight_col = 0, main = "Expresi�n Gen�tica", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE,
         fontsize_number = 6)

# Paso 21. Instalar y ejecutar la paqueter�a de viridis para tener disponibles las paletas de colores magma, plasma, cividis e inferno.
install.packages("viridis")
library(viridis)

# Paso 22. A�adir colores con las paletas cargadas.
pheatmap(sub, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, 
         annotation_col = annotation_col, treeheight_row = 0, treeheight_col = 0, 
         main = "Expresi�n Gen�tica", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE,
         fontsize_number = 6, col = viridis_pal(option = "plasma") (6))

# Paso 23. Cambiar colores de la paleta.
pheatmap(sub, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, 
         annotation_col = annotation_col, treeheight_row = 0, treeheight_col = 0, 
         main = "Expresi�n Gen�tica", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE,
         fontsize_number = 6, col = viridis_pal(option = "magma") (6))

pheatmap(sub, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, 
         annotation_col = annotation_col, treeheight_row = 0, treeheight_col = 0, 
         main = "Expresi�n Gen�tica", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE,
         fontsize_number = 6, col = viridis_pal(option = "viridis") (6))

# Paso 24.Identificar las distancias entre los genes.
dist(sub)

# Paso 25. Identificar el heatmap de la correlaci�n entre pacientes.
pheatmap(cor(sub))

# Paso 26. Identificar el heatmap de la correlaci�n entre genes con la matriz traspuesta.
trasp <- t(sub)
pheatmap(cor(trasp))

#####################################################################################################################
