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
##' show( GradientDescentSumaCuadradosClase )
NULL
loadModule("ModuloSumaCuadrados",TRUE)
