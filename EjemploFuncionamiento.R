
library(GradDesc)
m=5
n=2
A=matrix(rnorm(m*n),nrow=m)
y =c(2,3,2,1,1)

GradientDescentSumaCuadrados(A,y)

y_matrix =matrix(y,ncol=1)
solve(t(A)%*%A,t(A)%*%y_matrix)


# Ejemplo 2: Regresión con datos Iris
head(iris)
A <- as.matrix(iris[,2:4])
y <- iris[,1]

coefs <- GradientDescentSumaCuadrados(A,y)

# usando función lm
lm(data=iris[,-5], formula=Sepal.Length~.-1)

# Con término independiente: incluimos una columna de unos en la matriz del modelo
A <- cbind(Interc=rep(1,150), A)
coefs <- GradientDescentSumaCuadrados(A,y)
lm(data=iris[,-5], formula=Sepal.Length~.)


# TO DO. Ejemplo 3: Regresión de Eólica nacional España con datos de un grupo grande de parques (+500, +1000?) y control de tiempos


