# PaquetesConClase

## *Un Crash-course sobre interfaces de C++ en R* 

### Grupo de Usuarios R de Madrid 

### 10/12/2015 

### Presentación

El presente repositorio contiene la presentación [presentacion.html](presentacion.html) y los 2 paquetes de ejemplo mencionados en ella. 

### Instalación y Dependencias

Hay que tener instalado:

```r
dependencias = c("Rcpp","roxygen2","devtools","rmarkdown","knitr")
install.packages(dependencias)

# Para instalar la versión autocontenida del paquete:
devtools::install_github('bichopaulinho/PaquetesConClase',subdir='GradDesc')

# Para generar los archivos de ejemplo:
library(rmarkdown)
rmarkdown::render("Ejemplos/EjemploFuncionamiento.R")

```

