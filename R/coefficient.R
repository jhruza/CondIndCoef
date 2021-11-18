#' @import drf

multi.which <- function(X, x) {
  if (ncol(X) == 1) {
    idx <- which(X <= x)
  } else{
    idx <- vector()
    for (i in 1:nrow(X)) {
      if (all((X[i, ] <= x))) {
        idx[length(idx) + 1] <- i
      }
    }
  }
  return(idx)
}

#kernel function
ker <- function(x, y) {
  #return(exp(-(rowSums((x-y)^2))*1/2))  #gaussian kernel with sigma = 1/2
  return(rowSums(x * y)) #linear kernel
}



# method is by default 'cdf' which is the sampi coefficient, alternatively one can choose 'pdf' which is the mean of l1 difference of the weights.
#speedupmatrix is a parameter which can be pre-calculated when calling sampi several times with the same X to speed up the calculation such as in variable selection

#' @export
cond.ind.coef <-
  function(X,
           Y,
           Z = NULL,
           method = 'cdf',
           speedupmatrix = NULL) {

    #transform the arguments into matrices in case they are vectors or dataframes
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    #number of data points
    n <- dim(X)[1]


    #separate unconditional case and conditional case
    if (!is.null(Z)) {

#conditional case----------------------------------------------------------------------------

      if (dim(Y)[1] != n || dim(Z)[1] != n) {
        return(warning('X, Y ,Z must be of same length'))
      }

      #transform Z into matrix in case it is vector or dataframe
      Z <- as.matrix(Z)

      #growing forests
      drf.forest.Z <- drf::drf(Z, X)
      drf.forest.Y.Z <- drf::drf(X = cbind(Y, Z), X)

      #calculating weights
      w.Z <-drf::get_sample_weights(drf.forest.Z, newdata = Z)
      w.Z.Y <-drf::get_sample_weights(drf.forest.Y.Z, newdata = cbind(Y, Z))

      if (method == 'cdf') {

        #speedup matrix can be pre computed to speed up computation time
        #if coefficient is called several times with same X
        if (is.null(speedupmatrix)) {
          speedupmatrix <- matrix(0, n, n)
          for (i in 1:n) {
            help.idx <- multi.which(X, X[i, ])
            speedupmatrix[i, help.idx] <- 1
          }
        }

        res <-abs(rowSums(as.matrix(w.Z.Y * speedupmatrix)) - rowSums(as.matrix(w.Z *speedupmatrix)))

      } else if (method == 'pdf') {

        res <- rowSums(as.matrix(abs(w.Z.Y - w.Z))) / n

      } else if (method == 'mmd') {
        #TODO

      } else{
        return(warning('wrong method selected'))
      }

      res <- mean(res)

    }
    else{

#unconditional case----------------------------------------------------------------------

      drf.forest.Y <- drf::drf(Y, X)
      w.Y <- drf::get_sample_weights(drf.forest.Y, newdata = Y)

      #want to estimate mmd(p_x , p_{x |y}) using formula 2.1 of master thesis
      if (method == 'mmd') {

        # nxn matrix with entries K(x_i,x_j)
        kernelmatrix<-matrix(0, nrow = n, ncol = n)
        for (i in 1:n) {
            kernelmatrix[i,]<-ker(matrix(rep(X[i,],n),nrow = n, byrow=TRUE),X)

        }

        res <- vector()
        for (i in 1:n){
          res[length(res)+1]<- 1/n^2*sum(kernelmatrix)+sum( (w.Y[i,]%*%t(w.Y[i,]) *kernelmatrix) )-2/n*sum( as.vector(w.Y[i,])*kernelmatrix )
        }


      } else{
      #cdf difference
        res<-vector()
        for (j in 1:n) {
          idx <- multi.which(X, X[j, ])
          res[length(res) + 1] <- abs(length(idx) / n - sum(w.Y[j, idx]))
        }
      }
    }

    res <- mean(res)

    return(res)
  }

