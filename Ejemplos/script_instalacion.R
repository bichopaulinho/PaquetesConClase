#' ---
#' title: "Paquetes con C++lase. Un Crash Course de Interfaces de C++ en R con Rcpp"
#' author: Gonzalo Mateo (gonzmg88@gmail.com) & Paulino Tardáguila (tardaguila@hotmail.com)
#' date: "Madrid, 9 de diciembre de 2015"
#' output:
#'  html_document:
#'    theme: cerulean
#'    highlight: tango
#'    css: ../images/base.css
#' ---

#+ include=FALSE
knitr::opts_chunk$set(cache=FALSE, fig.align='center', message=FALSE, "warning"=FALSE)
rm(list = ls());gc()

#' ## Script de Instalación
library("devtools")
library("roxygen2")
library("testthat")
library(Rcpp)
#rm(list=ls())

#remove.packages("GradDesc",lib = .libPaths()[2])

#' Limpiamos todo para forzar a reinstalar
#+ eval=FALSE
system("rm man/*.Rd*")

devtools::clean_dll()
devtools::clean_source("src")
document(".")

#' Generamos el objeto como binario (util para compartir el paquete. Evitamos que haya que recompilarlo por nuestros usuarios)
#+ eval=FALSE
build(".", binary=T, args="--no-multiarch")

#' Para instalarlo basta esto porque el paquete ya esta descargado
#+ eval=FALSE
install.packages("../GradDesc_0.1_R_x86_64-redhat-linux-gnu.tar.gz",
                 #lib = .libPaths()[2],# Para que lo instale en la libreria local
                 repos = NULL)
