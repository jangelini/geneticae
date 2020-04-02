#'AMMI biplots with \pkg{ggplot2}
#'
#'Produces the AMMI biplot as an object of class 'ggplot'. It is possible to
#'customize it so that the stylistic attributes are to the user's liking.
#'
#'@param Data a data frame
#'@param genotype name of the column that contains the genotypes
#'@param environment name of the column that contains the environments
#'@param response name of the column that contains the response
#'@param rep name of the column that contains the replications.If this argument
#'  is NULL, there is no replications in the data.
#'@param Ncomp number of principal components that will be used in the analysis
#'@param type method. Either "AMMI", "rAMMI", "hAMMI", "gAMMI", "lAMMI" or
#'  "ppAMMI".
#'@param colGen colour for genotype attributes on biplot. Defaults to "gray"
#'@param colEnv colour for environment attributes on biplot. Defaults to
#'  "darkred"
#'@param colSegment colour for segment or circle lines. Defaults to "gray"
#'@param colHull colour for hull when type=6. Defaults to "gray"
#'@param sizeGen text size for genotype labels. Defaults to 4
#'@param sizeEnv text size for environment labels. Defaults to 4
#'@param largeSize text size to use for larger labels where type=5, used for the
#'@param titles logical. If TRUE then include automatically generated titles
#'@param footnote logical. If TRUE then include automatically generated footbote
#'  two selected genotypes, and where type=6, used for the outermost genotypes.
#'  Defaults to 4.5
#'
#'@return A biplot of class  \code{ggplot}
#'@references Rodrigues PC, Monteiro A and Lourenco VM (2015). \emph{A robust
#'  AMMI model for the analysis of genotype-by-environment data}.
#'  Bioinformatics, 32, 58–66.
#'@export
#'
#'@examples
#'
#'library(geneticae)
#'# Data without replication
#'data(yan.winterwheat)
#'GGE1 <- rAMMI(yan.winterwheat, genotype="gen",environment="env", response="yield",
#'              type = "AMMI")
#'GGE1
#'
#'# Data with replication
#'data(plrv)
#'GGE2 <- rAMMI(plrv, genotype="Genotype",environment="Locality", response="Yield",
#'              rep="Rep", type = "AMMI")
#'GGE2
#'
#'
#'@import dplyr
#'@importFrom MASS rlm
#'@importFrom pcaMethods robustSvd
#'@importFrom rrcov PcaHubert PcaGrid PcaLocantore PcaProj
#'@importFrom stats lm residuals
#'@importFrom dplyr group_by summarise rename pull
#'
rAMMI<-function(Data, genotype="gen", environment="env", response="Y", rep=NULL,Ncomp = 2, type = "AMMI",
                colGen="gray47",colEnv="darkred",colSegment="gray30",colHull="gray30",
                sizeGen=4,sizeEnv=4,largeSize=4.5, titles=TRUE, footnote=TRUE){

  if (missing(Data)) stop("Need to provide Data data frame or matrix")
  if(any(is.na(Data))){stop("Missing data in input data frame, run the imputation function first to complete the data set")}
  stopifnot(
    class(Data) == "data.frame",
    class(Ncomp) == "numeric",
    type %in% c("AMMI", "rAMMI", "hAMMI", "gAMMI", "lAMMI", "ppAMMI"),
    class(rep)%in% c("character", "NULL"),
    class(colGen) == "character",
    class(colEnv) == "character",
    class(colSegment) == "character",
    class(colHull) == "character",
    class(sizeGen) == "numeric",
    class(sizeEnv) == "numeric",
    class(largeSize) == "numeric",
    class(titles)  == "logical",
    class(footnote)  == "logical"
  )

  if(!is.null(rep)){
    Data <-
      Data %>%
      group_by(!!sym(genotype), !!sym(environment)) %>%
      summarise(y=mean(!!sym(response)))
    response = "y"
  }

  gen <- as.factor(pull(Data, genotype))
  env <- as.factor(pull(Data, environment))
  y <- pull(Data, response)

  # Asi funciona
  # gen <- as.factor(Data[, 1])
  # env <- as.factor(Data[, 2])
  # y <- as.vector(Data[, 3])
  Ngen <- nlevels(gen)
  Nenv <- nlevels(env)
  #Nreps <- nlevels(reps)
  labelgen <- unique(gen)
  labelenv <- unique(env)
  Ngen <- nlevels(gen)
  Nenv <- nlevels(env)
  Max <- min(Ngen - 1, Nenv - 1)
  lm.x <- lm(y ~ gen + env, data=Data)             # AMMI model - stage 1
  residuals.x <- matrix(residuals(lm.x), nrow = Ngen , ncol = Nenv)
  rlm.x <- rlm(y ~ gen + env, data=Data)           # Robust AMMi model - stage 1
  rresiduals.x <- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)

  if (type == "AMMI"){              # H-AMMI model
    svd.x<- svd((residuals.x))               # AMMI model - stage 2
    eigenvalues<-svd.x$d
    biplot.x.u<- -svd.x$u[,1:2] * matrix(rep(svd.x$d[1:2]^0.5, Ngen), ncol=2, byrow=T)
    biplot.x.v<- -svd.x$v[,1:2] * matrix(rep(svd.x$d[1:2]^0.5, Nenv), ncol=2, byrow=T)
    }

   if (type == "rAMMI"){
     rsvd.x<-robustSvd(rresiduals.x)      # r-AMMI model
     eigenvalues<- rsvd.x$d
     biplot.x.u<- cbind(rsvd.x$u[,1], -rsvd.x$u[,2]) * matrix(rep(rsvd.x$d[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.x.v<- cbind(rsvd.x$v[,1], -rsvd.x$v[,2]) * matrix(rep(rsvd.x$d[1:2]^0.5, Nenv), ncol=2, byrow=T)
   }

   if (type == "hAMMI"){              # H-AMMI model
     modeloHubert<- PcaHubert(rresiduals.x, mcd=FALSE)
     eigenvalues<- modeloHubert@eigenvalues
     loadings.Hx<- modeloHubert@loadings
     scores.Hx<- modeloHubert@scores
     biplot.x.u<- matrix(c(scores.Hx[,1], -scores.Hx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.x.v<- matrix(c(loadings.Hx[,1], -loadings.Hx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues[1:2]^0.5, Nenv), ncol=2, byrow=T)
   }

   if (type == "gAMMI"){             # G-AMMI model
     modeloGrid<- PcaGrid(rresiduals.x)
     eigenvalues<- modeloGrid@eigenvalues
     loadings.gx<- modeloGrid@loadings
     scores.gx<- modeloGrid@scores
     biplot.x.u<- scores.gx[,1:2] * matrix(rep(eigenvalues[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.x.v<- loadings.gx[,1:2] * matrix(rep(eigenvalues[1:2]^0.5, Nenv), ncol=2, byrow=T)
   }

   if (type == "lAMMI"){             # L-AMMI model
     modeloLocantore<- PcaLocantore(rresiduals.x)
     eigenvalues<- modeloLocantore@eigenvalues
     loadings.lx<- modeloLocantore@loadings
     scores.lx<- modeloLocantore@scores
     biplot.x.u<- matrix(c(-scores.lx[,1], scores.lx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.x.v<- matrix(c(-loadings.lx[,1], loadings.lx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues[1:2]^0.5, Nenv), ncol=2, byrow=T)
   }

   if (type == "ppAMMI"){             # PP-AMMI model
     modeloProj<- PcaProj(rresiduals.x)
     eigenvalues<- modeloProj@eigenvalues
     loadings.px<- modeloProj@loadings
     scores.px<- modeloProj@scores
     biplot.x.u<- matrix(c(-scores.px[,1], -scores.px[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.x.v<- matrix(c(-loadings.px[,1], -loadings.px[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues[1:2]^0.5, Nenv), ncol=2, byrow=T)
     }


  coordgenotype=biplot.x.u
  coordenviroment=biplot.x.v


  labelaxes <- paste("Component ",1:ncol(diag(eigenvalues)), sep = "")
  labelaxes=labelaxes[1:2]

  vartotal = round(as.numeric(sum(eigenvalues^2)),5)
  varexpl = c(round(as.numeric((eigenvalues[1]^2/vartotal) *100), 2),
              round(as.numeric((eigenvalues[2]^2/vartotal) *100), 2))

  plotdata<-data.frame(rbind(data.frame(coordgenotype,type="genotype",label=labelgen),
                             data.frame(coordenviroment,type="environment",label=labelenv)))
  colnames(plotdata)[1:2]<-c("Component1","Component2")
  plotdata$type<-factor(plotdata$type)


  AMMI1<-ggplot(data=plotdata,aes(x=Component1,y=Component2,group="type"))+
    theme_classic() +
    theme(plot.caption = element_text(size=12,hjust=0)) +
    scale_color_manual(values=c(colGen,colEnv)) +
    scale_size_manual(values=c(sizeGen,sizeEnv))+
    xlab(paste(labelaxes[1],format(varexpl[1],nsmall=2), "%", sep = " "))+
    ylab(paste(labelaxes[2], format(varexpl[2],nsmall=2),"%", sep = " "))+
    geom_hline(yintercept=0)+geom_vline(xintercept=0)

  AMMI2 <- AMMI1 + geom_segment(xend = 0, yend = 0, col = alpha(colEnv, 0.5),
                                data = subset(plotdata, type == "environment")) +
    geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)

  if(titles==TRUE){AMMI2<-AMMI2+ggtitle(type)}
  if(footnote==TRUE){
    footnotetxt=paste("\n", type,"biplot showing components 1 and 2 explaining ",sum(varexpl),"% of the total variation")

    AMMI2<-AMMI2+ labs(caption = footnotetxt)+theme(plot.caption = element_text(size=10,hjust=0,face="italic"))
  }


  return(AMMI2)
}
