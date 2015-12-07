library("devtools")
library("roxygen2")
library("testthat")
library(Rcpp)
#rm(list=ls())

# Instalacion y prueba del paquete
#actualizamos la documentacion


#remove.packages("predeolica",lib = .libPaths()[2])
#system("make")
system("rm man/*.Rd*")

devtools::clean_dll()
devtools::clean_source("src")
document(".")

build(".", binary=T, args="--no-multiarch")

# Para instalarlo basta esto porque el paquete ya esta descargado
install.packages("../GradDesc_0.1_R_x86_64-redhat-linux-gnu.tar.gz",
                 lib = .libPaths()[2],# Para que lo instale en la libreria local
                 repos = NULL)
