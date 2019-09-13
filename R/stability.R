#' Stability Function
#'
#' Performs a stability analysis
#' @param data A numeric data.frame with fenotypic means values of MET, and genotypes and environments by rows and columns, respectively.
#' @param type type of biplot to produce. fox, kang,nahu, thsu, AMMIindex
#' @export
#'
#'@importFrom stats cor median reorder sd
#'@importFrom calibrate textxy origin


stability<-function(data, type="fox"){

  if(type=="fox"){
    a <- nrow(data)
    b <- ncol(data)
    data.m <- as.data.frame(data)
    l.data <- length (data.m)
    ranks <- matrix(NA,a,b)
    y <- numeric()
    k <- numeric()
    for(i in 1:nrow(data.m)){
      for (j in 1:ncol(data.m)){

        ranks[,j] <- rank(-data.m[,j])
      }

      y <- which(ranks[i,] <= 3)
      k[i] <- length(y)


    }
    means <- round(as.numeric(rowMeans(data)),digits=4)
    result <- as.data.frame(cbind(rownames(data),means,k))
    colnames(result) <- c("Gen","Mean", "TOP")

    rank.y <- apply(-data,2,rank)
    ranks.sum.y <- apply(rank.y,1,sum)
    sd.rank = round(apply(rank.y,1,sd),digits=4)
    ranks.y = data.frame(rank.y,ranks.sum.y,sd.rank)
    colnames(ranks.y) = c(colnames(rank.y),"Sum", "Sd")
    cor.rank.y <- round(cor(ranks.y, method="pearson"), digits = 4)
    cor.rank.y.sp <- round(cor(ranks.y, method="spearman"), digits = 4)
    geral.list <- list("Fox"=result,"Ranks"=ranks.y,"Correlations Pearson"=cor.rank.y,"Correlations Spearman"=cor.rank.y.sp)


    return(geral.list)
}


if(type=="kang"){
  a <- nrow(data)
  b <- ncol(data)
  data.m <- as.matrix(data)
  l.data <- length (data.m)
  yield <- numeric()
  vars <- numeric()

  for (i in 1:nrow(data.m)){

    yield[i] <- sum(data.m[i,])
    vars[i] <- var(data.m[i,])
  }

  rank.y <- rank(-yield)
  rank.v <- rank(vars)
  rank.sum <- rank.y + rank.v
  rank.sum <- as.data.frame(rank.sum)

  means <- round(as.vector(rowMeans(data)),digits=4)
  result <- as.data.frame(cbind(rownames(data),means,rank.y,rank.v,rank.sum))
  colnames(result) <- c("Gen","Mean","Rank.Y","Rank.VAR","rank.sum")

  rank.y <- apply(-data,2,rank)
  ranks.sum.y <- apply(rank.y,1,sum)
  sd.rank = round(apply(rank.y,1,sd),digits=4)
  ranks.y = data.frame(rank.y,ranks.sum.y,sd.rank)
  colnames(ranks.y) = c(colnames(rank.y),"Sum", "Sd")
  cor.rank.y <- round(cor(ranks.y, method="pearson"), digits = 4)
  cor.rank.y.sp <- round(cor(ranks.y, method="spearman"), digits = 4)
  geral.list <- list("Kang"=result,"Ranks"=ranks.y,"Correlations Pearson"=cor.rank.y,"Correlations Spearman"=cor.rank.y.sp)


  return(geral.list)

}



if(type=="thsu"){
  a <- nrow(data)
  b <- ncol(data)
  data.m <- as.matrix(data)
  l.data <- length (data.m)
  data.r <- matrix(NA,a,b)
  ranks <- matrix(NA,a,b)
  ranks.y <- matrix(NA,a,b)
  np.1 <- matrix(NA,a,b)
  np.2 <- matrix(NA,a,b)
  np.3 <- matrix(NA,a,b)
  np.4 <- matrix(NA,a,b)
  np1 <- numeric()
  np2 <- numeric()
  np3 <- numeric()
  np4 <- numeric()
  k <- 2

  for (i in 1:nrow(data.m)){
    for (j in 1:ncol(data.m)){
      data.r[i,j] <- (data.m[i,j]) - (mean(data.m[i,]))
    }
  }

  for (j in 1:ncol(data.r)){

    ranks[,j] <- rank(-data.r[,j])
  }

  for (j in 1:ncol(data.m)){

    ranks.y[,j] <- rank(-data.m[,j])
  }

  for (i in 1:nrow(data)){
    for (j in 1:ncol(data)){

      np.1[i,j] <- abs((ranks[i,j])-median(ranks[i,]))
      np.2[i,j] <- abs((ranks[i,j])-median(ranks[i,]))/(median(ranks.y[i,]))
      np.3[i,j] <- ((ranks[i,j]-mean(ranks[i,]))^2) / b
    }
    np1[i] <- round((1/b) * (sum(np.1[i,])),digits=4)
    np2[i] <- round((1/b) * (sum(np.2[i,])),digits=4)
    np3[i] <- round((sqrt(sum(np.3[i,]))) / (mean(ranks.y[i,])),digits=4)
  }

  for (i in 1:nrow(data)){
    for (j in 1:(b-1)){

      np.4[i,j] <- abs(ranks[i,j] - ranks[i,k])
      while(k < b)
        k <- k + 1
    }
    np4[i] <- round((2/(b*(b-1))) * (sum((np.4[i,j]) / (mean(ranks.y[i,])))),digits=4)
  }

  means <- round(as.numeric(rowMeans(data)),digits=4)
  result <- as.data.frame(cbind(rownames(data),means,np1,np2,np3,np4))
  colnames(result) <- c("Gen","Mean","N1","N2","N3","N4")

  rank.y <- apply(-data,2,rank)
  ranks.sum.y <- apply(rank.y,1,sum)
  sd.rank = round(apply(rank.y,1,sd),digits=4)
  ranks.y = data.frame(rank.y,ranks.sum.y,sd.rank)
  colnames(ranks.y) = c(colnames(rank.y),"Sum", "Sd")
  cor.rank.y <- round(cor(ranks.y, method="pearson"), digits = 4)
  cor.rank.y.sp <- round(cor(ranks.y, method="spearman"), digits = 4)
  geral.list <- list("ThSu"=result,"Ranks"=ranks.y,"Correlations Pearson"=cor.rank.y,"Correlations Spearman"=cor.rank.y.sp)


  return(geral.list)
  }

  if(type=="nahu"){
    a <- nrow(data)
    b <- ncol(data)
    data.m <- as.matrix(data)
    l.data <- length (data.m)
    data.r <- matrix(NA,a,b)
    S2.1 <- matrix(NA,a,b)
    S3.1 <- matrix(NA,a,b)
    S6.1 <- matrix(NA,a,b)
    S2 <- numeric()
    S3 <- numeric()
    S6 <-  numeric()
    S.1 <- matrix(NA,a,b)
    S1 <- numeric()
    k <-  2

    for (i in 1:nrow(data.m)){
      for(j in 1:ncol(data.m)){
        data.r[i,j] <- (data.m[i,j]) - (mean(data.m[i,]))+(mean(data.m))
      }
    }
    ranks <- matrix(NA,a,b)

    for (j in 1:ncol(data.r)){

      ranks[,j] <- rank(data.r[,j])
    }

    ranks.y <- matrix(NA,a,b)

    for (j in 1:ncol(data.m)){

      ranks.y[,j] <- rank(data.m[,j])
    }

    r.means <- rowMeans(ranks)
    r.means.y <- rowMeans(ranks.y)

    for (i in 1:nrow(data)){
      for (j in 1:ncol(data)){
        S2.1[i,j] <- (ranks[i,j]-r.means[i])^2
        S3.1[i,j] <- (ranks.y[i,j]-r.means.y[i])^2
        S6.1[i,j] <- (abs(ranks.y[i,j]-r.means.y[i]))
      }
      S2[i]<-round(((sum(S2.1[i,])) / (b-1)),digits=4)
      S3[i]<-round((sum(S3.1[i,]) / (r.means.y[i])),digits=4)
      S6[i]<-round((sum(S6.1[i,]) / (r.means.y[i])),digits=4)
    }
    for (i in 1:nrow(data)){
      for (j in 1:(b-1)){

        S.1[i,j] <- abs(ranks[i,j] - ranks[i,k])
        while(k < b)
          k <- k + 1
      }
      S1[i] <- round((2*(sum(S.1[i,j]))) / (b*(b-1)), digits = 4)

    }

    means <- round(as.numeric(rowMeans(data)),digits=4)
    result <- as.data.frame(cbind(rownames(data),means,S1,S2,S3,S6))
    colnames(result) <- c("Gen","Mean","S1","S2","S3","S6")



    return(result)
  }
}
