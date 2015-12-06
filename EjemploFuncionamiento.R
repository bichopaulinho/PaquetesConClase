
library(GradDesc)
m=5
n=2
A=matrix(rnorm(m*n),nrow=m)
y =c(2,3,2,1,1)

GradientDescentSumaCuadrados(A,y)

y_matrix =matrix(y,ncol=1)
solve(t(A)%*%A,t(A)%*%y_matrix)
