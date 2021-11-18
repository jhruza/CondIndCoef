#testing

n=400

Z = matrix(runif(2 * n,-1,1), ncol = 2)

X = matrix(runif(2 * n,-1,1), ncol = 2)
Y = X^2+Z


cond.ind.coef(X,Y)
cond.ind.coef(X,Y,Z)
cond.ind.coef(X,Y, method = "mmd")
cond.ind.coef(X,Y,Z, method = "mmd")

Y = matrix(rnorm(2*n), ncol = 2)

cond.ind.coef(X,Y,Z)
cond.ind.coef(X,Y,Z, method = "mmd")

