generate_RPmatrix <- function(m,p,psi=0.1) {
  vals <- sample(c(-1,0,1), size=m*p, prob = c(psi/2,1-psi,psi/2), replace=TRUE)
  RM <- Matrix(vals,nrow=m,ncol=p,sparse=TRUE)
  zero_row <- !apply(RM,1,function(ro) sum(abs(ro))>0)
  if (any(zero_row)) { # choose 1 random 1 entry in each zero_rows
    for (j in which(zero_row)) {
      RM[j,sample(1:p,1)] <- 1
    }
  }
  return(RM)
}

generate_RPmatrixCW <- function(m,p) {
  RM <- Matrix(c(0),nrow=m,ncol=p,sparse=TRUE)
  for (j in 1:p) {
    RM[sample(1:m,1),j] <- sample(c(-1,1),1)
  }
  zero_row <- !apply(RM,1,function(ro) sum(abs(ro))>0)
  if (any(zero_row)) { # choose 1 random 1 entry in each zero_rows
    for (j in which(zero_row)) {
      RM[j,sample(1:p,1)] <- sample(c(-1,1),1)
    }
  }
  return(RM)
}

generate_RPM <- function(m,
                         p,
                         coef=sample(c(-1,1),p,replace = TRUE)) {
  goal_dims <- sample(1:m,p,replace = TRUE)
  counter <- 0
  # remove zero rows
  for (goal_dim in 1:m) {
    if (sum(goal_dims==(goal_dim-counter))==0) {
      goal_dims[goal_dims>goal_dim-counter] <- goal_dims[goal_dims>goal_dim-counter]-1
      counter <- counter + 1
    }
  }
  
  RM <- Matrix::Matrix(c(0),nrow=m-counter,ncol=p,sparse=TRUE)
  RM@i <- as.integer(goal_dims-1)
  RM@p <- 0:p
  RM@x <- coef
  
  return(RM)
}
