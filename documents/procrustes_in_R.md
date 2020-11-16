code used from R for procrustes. 
```
# .pGPA
# GPA with partial Procrustes superimposition
# used in gpagen
.pGpa <- function(Y, PrinAxes = FALSE, Proj = FALSE, max.iter = 5){
  iter <- 0
  n <- length(Y); p <- nrow(Y[[1]]); k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+",Ya)/n
  Q <- ss <- n*(1-sum(M^2))
  M <- cs.scale(M)
  iter <- 1
  while(Q > 0.0001){
    iter <- iter+1
    Ya <- apply.pPsup(M, Ya)
    M <- Reduce("+",Ya)/n
    ss2 <- n*(1-sum(M^2))
    Q <- abs(ss-ss2)
    ss <- ss2
    M <- cs.scale(M)
    if(iter > max.iter) break
  }
  if (PrinAxes == TRUE) {
    ref <- M
    rot <- prcomp(ref)$rotation
    for (i in 1:k) if (sign(rot[i, i]) != 1)
      rot[1:k, i] = -rot[1:k, i]
    Ya <- Map(function(y) y%*%rot, Ya)
    M <- cs.scale(Reduce("+", Ya)/n)
  }
  list(coords= Ya, CS=CS, iter=iter, consensus=M, Q=Q, nsliders=NULL)
}


# orp
# projection in GPA
# used in gpagen functions
orp<-function(A){
  if(is.array(A)) {
    dims <- dim(A)
    n <- dims[3]; k <- dims[2]; p <- dims[1]
    Y <- lapply(1:n, function(j) A[,,j])
  } else
    if(is.list(A)){
      Y <- A
      n <- length(A); dims <- dim(A[[1]]); k <- dims[2]; p <- dims[1]
    } else stop("Input must be either a list or array")

  Y1 <- as.vector(center.scale((Reduce("+", Y)/n))$coords)
  oo <- matrix(1,n)%*%Y1
  mat <- t(matrix(unlist(Y),k*p,n))
  Xp <- (mat%*%(diag(1,p*k) - (tcrossprod(Y1)))) +oo
  lapply(1:n, function(j) matrix(Xp[j,],p,k))
}

# .procD.slide
# same as procD.slide, but without progress bar option
# used in pGpa.wSliders
.procD.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  while(Q > 0.0001){
    iter <- iter+1
    if(!is.null(curves)) tans <- Map(function(y) tangents(curves, y, scaled=TRUE), slid0)
    if(is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.procD(y, tn, ref), tans, slid0)
    if(!is.null(surf) & is.null(curves))
      slid <- Map(function(y) semilandmarks.slide.surf.procD(y, surf, ref), slid0)
    if(!is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.surf.procD(y, tn, surf, ref), tans, slid0)
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    ref = cs.scale(Reduce("+", slid0)/n)
    Q <- abs(ss0-ss)
    ss0 <- ss
    if(iter >=max.iter) break
  }
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}

# .BE.slide
# same as BE.slide, but without progress bar option
# used in pGpa.wSliders
# performs sliding iterations using bending energy

.BE.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  while(Q > 0.0001){
    iter <- iter+1
    if(!is.null(curves)) tans <- Map(function(y) tangents(curves, y, scaled=TRUE), slid0)
    L <- Ltemplate(ref)
    if(is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.BE(y, tn, ref, L), tans, slid0)
    if(!is.null(surf) & is.null(curves))
      slid <- Map(function(y) semilandmarks.slide.surf.BE(y, surf, ref, L), slid0)
    if(!is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.surf.BE(y, tn, surf, ref, L), tans, slid0)
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    ref = cs.scale(Reduce("+", slid0)/n)
    Q <- abs(ss0-ss)
    ss0 <- ss
    if(iter >= max.iter) break
  }
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}

.pGpa.wSliders <- function(Y, curves, surf, ProcD = TRUE, PrinAxes = FALSE, Proj = FALSE, max.iter = 5){
  n <- length(Y); p <- nrow(Y[[1]]); k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+", Ya)/n
  if(ProcD == FALSE) gpa.slide <- .BE.slide(curves, surf, Ya, ref=M, max.iter=max.iter) else
    gpa.slide <- .procD.slide(curves, surf, Ya, ref=M, max.iter=max.iter)
  Ya <- gpa.slide$coords
  M <- gpa.slide$consensus
  iter <- gpa.slide$iter
  Q <- gpa.slide$Q
  if (PrinAxes == TRUE) {
    ref <- M
    rot <- prcomp(ref)$rotation
    for (i in 1:k) if (sign(rot[i, i]) != 1)
      rot[1:k, i] = -rot[1:k, i]
    Ya <- Map(function(y) y%*%rot, Ya)
    M <- center.scale(Reduce("+", Ya)/n)$coords
  }
  list(coords= Ya, CS=CS, iter=iter, consensus=M, Q=Q, nsliders=NULL)
}
```

Called using:
```
gpa <- .pGpa(Y, PrinAxes = PrinAxes, max.iter=max.it)
coords <- gpa$coords
  M <- gpa$consensus
  dimnames(M) <- list(p.names, k.names)
    
  if (Proj == TRUE) {
    coords <- orp(coords)
    M <- Reduce("+",coords)/n
    dimnames(M) <- list(p.names, k.names)
  }
```