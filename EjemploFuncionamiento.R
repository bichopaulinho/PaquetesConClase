#' ---
#' title: "Paquetes con C++lase. Un Crash Course de Interfaces de C++ en R con Rcpp"
#' author: Gonzalo Mateo (gonzmg88@gmail.com) & Paulino Tardáguila (tardaguila@hotmail.com)
#' date: "Madrid, 9 de diciembre de 2015"
#' output:
#'  html_document:
#'    toc: true
#'    theme: cerulean
#'    highlight: tango
#'    css: images/base.css
#' ---

#+ include=FALSE
knitr::opts_chunk$set(cache=FALSE, fig.align='center', message=FALSE, "warning"=FALSE)
rm(list = ls());gc()
#' ## Conjunto de Ejemplos de funcionamiento del paquete
#' El sitio ideal para estos ejemplos son las `vignettes` dentro del paquete. ¿[Como añadir vignettes a tu paquete](http://r-pkgs.had.co.nz/vignettes.html)?.

#+ eval=FALSE
browseVignettes("GradDesc")

#' ### Ejemplo 1
#' Aplicamos la formula de minimos cuadrados. Comprobamos que obtenemos los mismos resultados
library(GradDesc)
m=5
n=2
A=matrix(rnorm(m*n),nrow=m)
y =c(2,3,2,1,1)

xsol1 = GradientDescentSumaCuadrados(A,y)
xsol1

y_matrix =matrix(y,ncol=1)
xsol2 = solve(t(A)%*%A,t(A)%*%y_matrix)
t(xsol2)

sum(abs(xsol1 - as.numeric(xsol2)))<1e-6

#' ### Ejemplo 2
#' Usamos el `RCPP_MODULE`. Tenemos una clase RC de R [(¿qué es eso?)](http://adv-r.had.co.nz/OO-essentials.html#rc) llamada `GradientDescentSumaCuadradosClase`.
#' `Rcpp` nos ha generado esta clase de forma automática con la macro `RCPP_MODULE`.
#'
#' Con un objeto de esta clase vamos a cambiar el punto inicial y al llamar a `solve` (aplicar la optimización).
#'

#' Para mostrar el contenido de la clase:
show( GradientDescentSumaCuadradosClase )

#' Creamos el objeto
obj = new(GradientDescentSumaCuadradosClase,A,y)

punto_inicial_alternativo = c(2,2)
obj$starting_point(punto_inicial_alternativo)

#' Podemos comprobar que realmente hemos cambiado el punto inicial en la evaluación inicial. (Esto solo sale en la consola de R)
evaluacion_inicial = sum((A%*%matrix(punto_inicial_alternativo,nrow=2)-matrix(y,ncol=1))^2)
obj$solve()

#' ### Ejemplo 3: Regresión con datos Iris
head(iris)
A <- as.matrix(iris[,2:4])
y <- iris[,1]

coefs <- GradientDescentSumaCuadrados(A,y)
coefs
#' usando función lm (quitamos el término independiente)
lm(data=iris[,-5], formula=Sepal.Length~.-1)

#' Con término independiente: incluimos una columna de unos en la matriz del modelo
A <- cbind(Interc=rep(1,150), A)
coefs <- GradientDescentSumaCuadrados(A,y)
coefs
lm(data=iris[,-5], formula=Sepal.Length~.)


# TO DO. Ejemplo 3: Regresión de Eólica nacional España con datos de un grupo grande de parques (+500, +1000?) y control de tiempos
#' Información ejecución
devtools::session_info()
