##' GradDesc-package: paquete wrapper de la biblioteca Purple
##'
##' La documentacion de la clase \code{GradientDescentSumaCuadradosClase} está en: \code{\link{GradientDescentSumaCuadradosClase}}
##'
##' @name GradDesc-package
##' @docType package
##' @useDynLib GradDesc
##' @import Rcpp
##' @title GradDesc
##' @author Gonzalo Mateo, Paulino Tardáguila
NULL

##' Rcpp module: ModuloSumaCuadrados Clase GradientDescentSumaCuadradosClase
##'
##' clase que envuelve la clase de \code{Purple} \code{Purple::GradientDescent}.
##' Usar \code{show} para ver los métodos que tiene implementados.
##'
##'
##' @name GradientDescentSumaCuadradosClase
##' @details
##'  \enumerate{
##'  \item La función \code{setGradientNormGoal} Establece el objetivo de la norma (hasta el cual se minimiza)
##'  \item La función \code{solve} Ejecuta la optimizacion
##'  \item La función \code{starting_point} Cambia el punto de inicio de la optimizacion
##'  \item La función \code{setEvaluationGoal} Establece el objetivo hasta el que se optimiza
##' }
##' @export
##' @examples
##' #' Para mostrar el contenido del objeto
##' show( GradientDescentSumaCuadradosClase )
##' #' Creamos el objeto
##' m=5
##' n=2
##' A=matrix(rnorm(m*n),nrow=m)
##' y =c(2,3,2,1,1)
##'
##' obj = new(GradientDescentSumaCuadradosClase,A,y)
##'
##' punto_inicial_alternativo = c(2,2)
##' obj$starting_point(punto_inicial_alternativo)
##'
##' #' Podemos comprobar que realmente hemos cambiado el punto inicial en la evaluación inicial. (Esto solo sale en la consola de R)
##' evaluacion_inicial = sum((A%*%matrix(punto_inicial_alternativo,nrow=2)-matrix(y,ncol=1))^2)
##' obj$solve()
NULL
loadModule("ModuloSumaCuadrados",TRUE)
