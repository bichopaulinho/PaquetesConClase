# PaquetesConClase

## *Un Crash-course sobre interfaces de C++ en R* 

### Grupo de Usuarios R de Madrid 

### 10/12/2015 

### Instalación y Dependencias

Hay que tener instalado:

```r
dependencias = c("Rcpp","roxygen2","devtools","rmarkdown","knitr")
install.packages(dependencias)

# Para generar los archivos de ejemplo:
library(rmarkdown)
rmarkdown::render("EjemploFuncionamiento.R")
```

