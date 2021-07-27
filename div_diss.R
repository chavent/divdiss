div_diss <- function (data, D, K = NULL, w=rep(1. / nrow(data), nrow(data))) 
{
  obj <- split_mix(data)
  data_quanti <- obj$data_quanti
  data_quali <- obj$data_quali
  
  # quanti preprocessing
  if (!is.null(data_quanti)) {
    n <- nrow(data_quanti)
    nb_quanti <- ncol(data_quanti)
    rnames <- rownames(data_quanti)
    X_quanti <- as.matrix(data_quanti)
  } else {
    n <- nrow(data_quali)
    X_quanti <-matrix(0,n,0) 
    nb_quanti <- 0
  }     
  # quali preprocessing
  if (!is.null(data_quali)) {
    if (sum(unlist(lapply(data_quali,is.character))) !=0)
      stop("All columns in data_quali must be of class factor")
    res <- format_cat_data(data_quali)
    rnames <- rownames(data_quali) 
    mod_quali <- res$mod_quali
    vec_quali <- res$vec_quali
    cnames_quali <- res$cnames_quali
    X_quali <- as.matrix(res$Y)
    nb_quali <- length(vec_quali)
    vec_order <- unlist(lapply(data_quali,function(x){class(x)[1]=="ordered"})) 
  } else {
    X_quali <-matrix(0,n,0)
    mu_quali <- matrix(0,n,0) 
    cnames_quali <- list()
    nb_quali <- 0
    mod_quali <- list()
    vec_quali <- c()
    vec_order <- c() 
  }    
  
  # preprocessing
  X <- cbind(X_quanti, X_quali) 
  p <- ncol(X)


  # tree construction
  # algorithm's sketch : 
  # start with one leaf consisting in best question chosen on all indices
  # leaves <- {Tree(choose_qestion(Z, all_indices))}
  # each iteration, compute "locally maximal partition"
  # leaves <- (leaves\\{Fmax}) U {Tree(A_l(F_max), Tree(Al_c(F_max))}
  
  dendro <-new.env()
  dendro$v <- choose_question_diss(X, D, c(1:n), vec_quali, w, vec_order)
  leaves <- list(dendro)
  k <- 1 
  kmax <- nrow(unique(X))
  if (is.null(K)) K <- kmax
  if (!(K %in% 2:kmax)) K <- kmax
  height <- rep(NA,K-1)
 
  while (k < K) {
    ind_max <- which.max(lapply(leaves, function(x){x$v$inert}))
    height[k] <- leaves[[ind_max]]$v$inert
    Fmin <- leaves[[ind_max]]
    leaves[[ind_max]] <- NULL
    
    cqg <- choose_question_diss(X, D, Fmin$v$A_l, vec_quali, w, vec_order)
    Fmin$l <- new.env()
    Fmin$l$v <- cqg
    leaves <- append(leaves, Fmin$l)
        
    cqd <- choose_question_diss(X, D, Fmin$v$A_l_c, vec_quali, w, vec_order)
    Fmin$r <- new.env() 
    Fmin$r$v<- cqd
    leaves <- append(leaves, Fmin$r)
    
    k <- k+1
  }
  
  cluster <-list(tree=dendro)
  cluster$X <- X
  cluster$mod_quali <- mod_quali
  cluster$vec_quali <- vec_quali
  cluster$nb_quanti <- nb_quanti
  cluster$rnames <- rnames
  cluster$cnames_quanti <- colnames(data_quanti)
  cluster$cnames_quali <- cnames_quali
  cluster$kmax <- kmax
  cluster$height <- height
  cluster$Total <- inertdiss(D,wt=w)
  cluster$B <- sum(height)/cluster$Total
  ret_leaves <- list_leaves(cluster)

  cluster$description <- make_description(ret_leaves$leaves,cluster)
  names(cluster$description) <- paste("C", 1:K, sep = "")
  cluster$clusters <- lapply(ret_leaves$leaves, function(x) {rnames[x$class]})
  names(cluster$clusters) <- paste("C", 1:K, sep = "")
  cluster$which_cluster <- ret_leaves$classes
  cluster$data_quanti <- data_quanti
  cluster$data_quali <- data_quali
  class(cluster) <- "divclust"
  return(cluster)
}

