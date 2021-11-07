

library(drf)

multi.which<-function(X,x){
  if (ncol(X)==1) {
    idx<-which(X<=x)
  }else{
    idx<-vector()
    for (i in 1:nrow(X)) {
      if (all((X[i,] <=x))){
        idx[length(idx)+1]<-i
      }
    }
  }
  return(idx)
}


# method is by default 'cdf' which is the sampi coefficient, alternatively one can choose 'pdf' which is the mean of l1 difference of the weights.
#speedupmatrix is a parameter which can be pre-calculated when calling sampi several times with the same X to speed up the calculation such as in variable selection

#' @export
cond.ind.coef <- function(X,Y,Z=NULL, method='cdf', speedupmatrix=NULL) {
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  n<-dim(X)[1]

  if (!is.null(Z)) {


    Z<-as.matrix(Z)

    if (dim(Y)[1]!= n || dim(Z)[1]!=n) {
      return(warning('X, Y ,Z must be of same length'))
    }


    if (method=='cdf') {

      drf.forest.Z <- drf(Z, X)

      drf.forest.Y.Z <- drf(X=cbind(Y,Z), X)


      w.Z.Y<-get_sample_weights(drf.forest.Y.Z, newdata=cbind(Y,Z))
      w.Z<-get_sample_weights(drf.forest.Z, newdata = Z)


      if (is.null(speedupmatrix)) {
        speedupmatrix<-matrix(0,n,n)
        for (i in 1:n) {
          help.idx<-multi.which(X,X[i,])
          speedupmatrix[i,help.idx]<-1
        }
      }


      res<-abs(rowSums(as.matrix(w.Z.Y*speedupmatrix))-rowSums(as.matrix(w.Z*speedupmatrix)))





    }else if (method=='pdf') {
      drf.forest.Z <- drf(Z, X)

      drf.forest.Y.Z <- drf(X=cbind(Y,Z), X)

      w.Z.Y<-get_sample_weights(drf.forest.Y.Z, newdata=cbind(Y,Z))
      w.Z<-get_sample_weights(drf.forest.Z, newdata = Z)
      res<-rowSums(as.matrix(abs(w.Z.Y-w.Z)))/n
    }else{
      return(warning('wrong method selected'))
    }

    res<-mean(res)

  }
  else{
    drf.forest.Y <- drf(Y, X)
    w.Z.Y<-get_sample_weights(drf.forest.Y, newdata = Y)
    res<-vector()
    for (j in 1:n) {
      idx<-multi.which(X,X[j,])
      res[length(res)+1]<-abs(length(idx)/n - sum(w.Z.Y[j,idx]))
    }

    res<-mean(res)

  }


  return(res)
}
