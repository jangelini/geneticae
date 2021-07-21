#'AMMI biplots with \pkg{ggplot2}
#'
#'Produces clasic or robust AMMI biplot as an object of class 'ggplot'. It is
#'possible to customize it so that the stylistic attributes are to the user's
#'preferences.
#'
#'@param Data a dataframe with genotypes, environments, repetitions (if any) and
#'  the phenotypic trait of interest. There is no restriction on the order in
#'  which these variables should be presented in the dataframe, and also other
#'  variables that will not be used in the analysis can be included.
#'@param genotype column name containing genotypes.
#'@param environment column name containing environments.
#'@param response column name containing phenotypic trait.
#'@param rep column name containing replications. If this argument is NULL,
#'  there is no replication available on the data.
#'@param Ncomp number of principal components that will be used in the analysis.
#'@param type method. Either "AMMI", "rAMMI", "hAMMI", "gAMMI", "lAMMI" or
#'  "ppAMMI". Default to "AMMI".
#'@param colGen genotype attributes colour. Defaults to "gray".
#'@param colEnv environment attributes colour. Defaults to "darkred".
#'@param sizeGen genotype labels text size. Defaults to 4.
#'@param sizeEnv environment labels text size. Defaults to 4.
#'@param axis_expand multiplication factor to expand the axis limits by to
#'  enable fitting of labels. Defaults to 1.2.
#'@param titles logical, if this argument is TRUE a plot title is generated.
#'  Defaults to TRUE.
#'@param footnote logical, if this argument is TRUE a footnote is generated.
#'  Defaults to TRUE.
#'@param limits logical. If TRUE then automatically rescale axes. Defaults to
#'  TRUE.
#'@param axes logical, if this argument is TRUE x and y axes going through the
#'  origin. Defaults to TRUE.
#'@param axislabels logical, if this argument is TRUE labels axes are included.
#'  Defaults to TRUE.
#'
#'@return A biplot of class  \code{ggplot}
#'@references Rodrigues P.C., Monteiro A., Lourenco V.M. (2015). \emph{A robust
#'  AMMI model for the analysis of genotype-by-environment data}. Bioinformatics
#'  32, 58–66.
#'@export
#'
#'@examples
#'
#'library(geneticae)
#'# Data without replication
#'data(yan.winterwheat)
#'BIP_AMMI <- rAMMI(yan.winterwheat, genotype="gen",environment="env", response="yield",
#'              type = "AMMI")
#'BIP_AMMI
#'
#'# Data with replication
#'data(plrv)
#'BIP_AMMI2 <- rAMMI(plrv, genotype="Genotype",environment="Locality", response="Yield",
#'              rep="Rep", type = "AMMI")
#'BIP_AMMI2
#'
#'
#'@importFrom MASS rlm
#'@importFrom pcaMethods robustSvd
#'@importFrom rrcov PcaHubert PcaGrid PcaLocantore PcaProj
#'@importFrom stats lm residuals
#'@importFrom dplyr group_by summarise rename pull
#'
rAMMI<-function(Data, genotype="gen", environment="env", response="Y", rep=NULL,Ncomp = 2, type = "AMMI",
                colGen="gray47",colEnv="darkred",sizeGen=4,sizeEnv=4,titles=TRUE, footnote=TRUE, axis_expand=1.2, limits=TRUE,
                axes=TRUE,axislabels=TRUE){

  if (missing(Data)) stop("Need to provide Data data frame")
  if(any(is.na(Data))){stop("Missing data in input data frame, run the imputation function first to complete the data set")}
  stopifnot(
    class(Data) == "data.frame",
    class(genotype) == "character",
    class(environment) == "character",
    class(response) == "character",
    class(rep)%in% c("character", "NULL"),
    class(Ncomp) == "numeric",
    type %in% c("AMMI", "rAMMI", "hAMMI", "gAMMI", "lAMMI", "ppAMMI"),
    class(rep)%in% c("character", "NULL"),
    class(colGen) == "character",
    class(colEnv) == "character",
    class(sizeGen) == "numeric",
    class(sizeEnv) == "numeric",
    class(titles)  == "logical",
    class(footnote)  == "logical",
    class(axis_expand) == "numeric",
    class(limits) == "logical",
    class(axes)  == "logical",
    class(axislabels)  == "logical"
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

  Ngen <- nlevels(gen)
  Nenv <- nlevels(env)
  labelgen <- unique(gen)
  labelenv <- unique(env)
  Ngen <- nlevels(gen)
  Nenv <- nlevels(env)
  Max <- min(Ngen - 1, Nenv - 1)
  lm.x <- lm(y ~ gen + env, data=Data)             # AMMI model - stage 1
  residuals.x <- matrix(residuals(lm.x), nrow = Ngen , ncol = Nenv)
  rlm.x <- rlm(y ~ gen + env, data=Data)           # Robust AMMi model - stage 1
  rresiduals.x <- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)

  if (type == "AMMI"){                   # AMMI model
    svd.x<- svd((residuals.x))
    singlevalue<-svd.x$d
    loading<-(svd.x$v)
    score<-(svd.x$u%*%diag(svd.x$d))
  }

  if (type == "rAMMI"){
    rsvd.x<-robustSvd(rresiduals.x)      # r-AMMI model
    singlevalue<- rsvd.x$d
    loading<-(rsvd.x$v)
    score<-(rsvd.x$u%*%diag(rsvd.x$d))
  }

  if (type == "hAMMI"){                 # H-AMMI model
    modeloHubert<- PcaHubert(rresiduals.x, mcd=FALSE)
    singlevalue<- sqrt(modeloHubert@eigenvalues)
    loading<- modeloHubert@loadings
    score<- modeloHubert@scores
    n <- NROW(score)
    singlevalue <- singlevalue * sqrt(n)
  }

  if (type == "gAMMI"){             # G-AMMI model
    modeloGrid<- PcaGrid(rresiduals.x)
    singlevalue<- sqrt(modeloGrid@eigenvalues)
    loading<- modeloGrid@loadings
    score<- modeloGrid@scores
    n <- NROW(score)
    singlevalue <- singlevalue * sqrt(n)
  }

  if (type == "lAMMI"){             # L-AMMI model
    modeloLocantore<- PcaLocantore(rresiduals.x)
    singlevalue<- modeloLocantore@eigenvalues
    loading<- modeloLocantore@loadings
    score<- modeloLocantore@scores
    n <- NROW(score)
    singlevalue <- singlevalue * sqrt(n)
  }

  if (type == "ppAMMI"){             # PP-AMMI model
    modeloProj<- PcaProj(rresiduals.x)
    singlevalue<- sqrt(modeloProj@eigenvalues)
    loading<- modeloProj@loadings
    score<- modeloProj@scores
    n <- NROW(score)
    singlevalue <- singlevalue * sqrt(n)
  }

  scaling=0.5
  lambda_scaling <- singlevalue^scaling
  scores_scaling <- t(t(score) / lambda_scaling)
  loading_scaling <- t(t(loading) * lambda_scaling)

  coordgenotype <- scores_scaling
  coordenviroment <- loading_scaling


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
    scale_size_manual(values=c(sizeGen,sizeEnv))


  AMMI2 <- AMMI1 + geom_segment(xend = 0, yend = 0, col = alpha(colEnv, 0.5),
                                data = subset(plotdata, type == "environment")) +
    geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)

  if(titles==TRUE){AMMI2<-AMMI2+ggtitle(type)}
  if(footnote==TRUE){
    footnotetxt=paste("\n", type,"biplot showing components 1 and 2 explaining ",sum(varexpl),"% of the total variation")

    AMMI2<-AMMI2+ labs(caption = footnotetxt)+theme(plot.caption = element_text(size=8,hjust=0,face="italic"))
  }

  if(limits==TRUE){
    xlim<-c(min(plotdata$Component1*axis_expand),max(plotdata$Component1*axis_expand))
    ylim<-c(min(plotdata$Component2*axis_expand),max(plotdata$Component2*axis_expand))
    if(which(c(diff(xlim),diff(ylim))==max(c(diff(xlim),diff(ylim))))==1){
      xlim1<-xlim
      ylim1<-c(ylim[1]-(diff(xlim)-diff(ylim))/2,ylim[2]+(diff(xlim)-diff(ylim))/2)
    }
    if(which(c(diff(xlim),diff(ylim))==max(c(diff(xlim),diff(ylim))))==2){
      ylim1<-ylim
      xlim1<-c(xlim[1]-(diff(ylim)-diff(xlim))/2,xlim[2]+(diff(ylim)-diff(xlim))/2)
    }

    AMMI2<-AMMI2+scale_x_continuous(limits=xlim1,expand=c(0,0))+
      scale_y_continuous(limits=ylim1,expand=c(0,0))+
      coord_fixed(1)
  }

  if(axes==TRUE){
    AMMI2<-AMMI2+geom_hline(yintercept=0)+geom_vline(xintercept=0)
  }

  if(axislabels==TRUE){
    AMMI2<-AMMI2+
      xlab(paste(labelaxes[1],format(varexpl[1],nsmall=2), "%", sep = " "))+
      ylab(paste(labelaxes[2], format(varexpl[2],nsmall=2),"%", sep = " "))
  }

  return(AMMI2)
}
