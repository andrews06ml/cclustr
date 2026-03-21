# install.packages("devtools")   # si no lo tienes
# library(devtools)

# Primera parte
# Dirección donde quedara la libreria
#setwd("C:/Users/maria/OneDrive/Documentos/1. Andrés/MAESTRIA UBOSQUE/Proyecto de grado/Libreria")

#Crear la libreria local
#usethis::create_package("cclustr")

# Segunda parte
#devtools::document() #para actualizar comentarios
#devtools::load_all() #para actualizar cambios en funciones

?as_mild_list
?cluster_imputations
?consensus_clustering
?validate_clustering
?choose_best_clustering
?plot_consensus_dendrogram
?plot_consensus_matrix
?plot_validation_metrics

#-----------------------------------------#
# Pruebas
#-----------------------------------------#

#-----------------------------------------#
# Prueba para la primera función (imputación)
#-----------------------------------------#
basetotal<- readxl::read_xlsx("C:/Users/maria/OneDrive/Documentos/1. Andrés/MAESTRIA UBOSQUE/Proyecto de grado/Funciones con Documentación Oxygen/WDI.xlsx")
base <- basetotal[,c(5:12)]

library(mice)

# Imputación múltiple con mice (caso 1)
mice_fit <- mice(base, m = 5, method = 'pmm', seed = 123, maxit = 20)

completed <- mice::complete(mice_fit, action = "all")

sapply(completed, function(df) sum(is.na(df)))

# Utilización de la función
test_mice <-as_mild_list(mice_fit)

#-----------------------------------------#
# Prueba para la segunda función (agrupamiento)
#-----------------------------------------#
clustering <- cluster_imputations(test_mice, method = "ward.D2", k = 2:5)

cl <- clustering[["k3"]][["1"]]
plot(df1[,1:2], col = cl, pch = 19)

#-----------------------------------------#
# Prueba para la tercera función (coasignación)
#-----------------------------------------#


#-----------------------------------------#
# Prueba para la cuarta función (validación)
#-----------------------------------------#


#-----------------------------------------#
# Prueba para la quinta función (selección)
#-----------------------------------------#

#-----------------------------------------#
# Prueba para la sexta función (dendograma)
#-----------------------------------------#

#-----------------------------------------#
# Prueba para la septima función (mapa de calor)
#-----------------------------------------#

#-----------------------------------------#
# Prueba para la octava función (validacion por k)
#-----------------------------------------#
