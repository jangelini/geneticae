#'AMMI biplots with \pkg{ggplot2}
#'
#'Produces the AMMI biplot as an object of class 'ggplot'. It is possible to
#'customize it so that the stylistic attributes are to the user's liking.
#'
#'@param x a data frame that contains the genotypes in the first column, the columns in the
#'  next and last the response. In case of replications, the data frame must have the columns
#'  in the following order: genotypes, environments, replications and the response variable.
#'@param Ncomp number of principal components that will be used in the analysis
#'@param type method. Either "AMMI", "rAMMI", "hAMMI", "gAMMI", "lAMMI" or "ppAMMI".
#'@param rep logical. If TRUE the genotype by environment means is calculated.
#'@param colGen colour for genotype attributes on biplot. Defaults to "gray"
#'@param colEnv colour for environment attributes on biplot. Defaults to "darkred"
#'@param colSegment colour for segment or circle lines. Defaults to "gray"
#'@param colHull colour for hull when type=6. Defaults to "gray"
#'@param sizeGen text size for genotype labels. Defaults to 4
#'@param sizeEnv text size for environment labels. Defaults to 4
#'@param largeSize text size to use for larger labels where type=5, used for the
#'@param titles logical. If TRUE then include automatically generated titles
#'@param footnote logical. If TRUE then include automatically generated footbote
#'  two selected genotypes, and where type=6, used for the outermost genotypes.
#'  Defaults to 4.5
#'@return A biplot of class  \code{ggplot}
#'@references Rodrigues PC, Monteiro A and Lourenco VM (2015). \emph{A robust AMMI model for the analysis of
#'genotype-by-environment data}. Bioinformatics, 32, 58–66.
#'@export
#'@import dplyr
#'@importFrom MASS rlm
#'@importFrom pcaMethods robustSvd
#'@importFrom rrcov PcaHubert PcaGrid PcaLocantore PcaProj
#'@importFrom stats lm residuals
#'
rAMMI<-function(x, Ncomp = 2, type = "AMMI", rep=FALSE,
                colGen="gray47",colEnv="darkred",colSegment="gray30",colHull="gray30",
                sizeGen=4,sizeEnv=4,largeSize=4.5, titles=TRUE, footnote=TRUE){

  if (missing(x)) stop("Need to provide x data frame or matrix")

  # stopifnot(
  #   class(x) == "data.frame",
  #   class(Ncomp) == "numerical",
  #   class(type) %in% c("AMMI", "rAMMI", "hAMMI", "gAMMI", "lAMMI", "ppAMMI"),
  #   class(rep)  == "logical",
  #   class(colGen) == "character",
  #   class(colEnv) == "character",
  #   class(colSegment) == "character",
  #   class(colHull) == "character",
  #   class(sizeGen) == "numerical",
  #   class(sizeEnv) == "numerical",
  #   class(largeSize) == "numerical",
  #   class(titles)  == "logical",
  #   class(footnote)  == "logical"
  # )


  if(rep==TRUE){
    x <-
      x %>%
      group_by(x[,1], x[,2]) %>%
      summarise(y = mean(y))

    x <-as.data.frame(x)
    gen <- as.factor(x[,1])
    env <- as.factor(x[, 2])
    y <- as.vector(x[, 3])
    labelgen <- unique(gen)
    labelenv <- unique(env)
    Ngen <- nlevels(gen)
    Nenv <- nlevels(env)
    Max <- min(Ngen - 1, Nenv - 1)
    lm.x<- lm(y~gen+env, data=x)             # AMMI model - stage 1
    residuals.x<- matrix(residuals(lm.x), nrow = Ngen , ncol = Nenv)
    rlm.x<- rlm(y~gen+env, data=x)           # Robust AMMi model - stage 1
    rresiduals.x<- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)
  }

  if(rep==FALSE){
    x <-as.data.frame(x)
    gen <- as.factor(x[,1])
    env <- as.factor(x[, 2])
    y <- as.vector(x[, 3])
    labelgen <- unique(gen)
    labelenv <- unique(env)
    Ngen <- nlevels(gen)
    Nenv <- nlevels(env)
    Max <- min(Ngen - 1, Nenv - 1)
    lm.x<- lm(y~gen+env, data=x)             # AMMI model - stage 1
    residuals.x<- matrix(residuals(lm.x), nrow = Ngen , ncol = Nenv)
    rlm.x<- MASS::rlm(y~gen+env, data=x)           # Robust AMMi model - stage 1
    rresiduals.x<- matrix(residuals(rlm.x), nrow = Ngen, ncol = Nenv)
  }

  if (type == "AMMI"){              # H-AMMI model
    svd.x<- svd((residuals.x))               # AMMI model - stage 2                        # AMMI2 biplot
    biplot.x.u<- -svd.x$u[,1:2] * matrix(rep(svd.x$d[1:2]^0.5, Ngen), ncol=2, byrow=T)
    biplot.x.v<- -svd.x$v[,1:2] * matrix(rep(svd.x$d[1:2]^0.5, Nenv), ncol=2, byrow=T)

    coordgenotype=biplot.x.u
    coordenviroment=biplot.x.v
    eigenvalues = svd.x$d

    labelaxes <- paste("Component ",1:ncol(diag(svd.x$d)), sep = "")
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

     AMMI2<-AMMI1+geom_segment(xend=0,yend=0,col=alpha(colEnv,0.5),data=subset(plotdata,type=="environment"))+
      geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)
     if(titles==TRUE){AMMI2<-AMMI2+ggtitle("AMMI Biplot")}

    }


   if (type == "rAMMI"){
     rsvd.x<-robustSvd(rresiduals.x)      # r-AMMI model
     biplot.rx.u<- cbind(rsvd.x$u[,1], -rsvd.x$u[,2]) * matrix(rep(rsvd.x$d[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.rx.v<- cbind(rsvd.x$v[,1], -rsvd.x$v[,2]) * matrix(rep(rsvd.x$d[1:2]^0.5, Nenv), ncol=2, byrow=T)
     coordgenotype=biplot.rx.u
     coordenviroment=biplot.rx.v
     eigenvalues = rsvd.x$d

     labelaxes <- paste("Component ",1:ncol(diag(rsvd.x$d)), sep = "")
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

     AMMI2<-AMMI1+geom_segment(xend=0,yend=0,col=alpha(colEnv,0.5),data=subset(plotdata,type=="environment"))+
       geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)

     if(titles==TRUE){AMMI2<-AMMI2+ggtitle("rAMMI Biplot")}
   }


   if (type == "hAMMI"){              # H-AMMI model
     modeloHubert<- PcaHubert(rresiduals.x, mcd=FALSE)
     eigenvalues.Hx<- modeloHubert@eigenvalues
     loadings.Hx<- modeloHubert@loadings
     scores.Hx<- modeloHubert@scores
     biplot.Hx.u<- matrix(c(scores.Hx[,1], -scores.Hx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues.Hx[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.Hx.v<- matrix(c(loadings.Hx[,1], -loadings.Hx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues.Hx[1:2]^0.5, Nenv), ncol=2, byrow=T)

     coordgenotype=biplot.Hx.u
     coordenviroment=biplot.Hx.v
     eigenvalues = eigenvalues.Hx

     labelaxes <- paste("Component ",1:ncol(diag(eigenvalues.Hx)), sep = "")
     labelaxes=labelaxes[1:2]

     vartotal = round(as.numeric(sum(eigenvalues)),5)
     varexpl = c(round(as.numeric((eigenvalues[1]/vartotal) *100), 2),
                 round(as.numeric((eigenvalues[2]/vartotal) *100), 2))

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

     AMMI2<-AMMI1+geom_segment(xend=0,yend=0,col=alpha(colEnv,0.5),data=subset(plotdata,type=="environment"))+
       geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)
     if(titles==TRUE){AMMI2<-AMMI2+ggtitle("hAMMI Biplot")}

    }
   if (type == "gAMMI"){             # G-AMMI model
     modeloGrid<- PcaGrid(rresiduals.x)
     eigenvalues.gx<- modeloGrid@eigenvalues
     loadings.gx<- modeloGrid@loadings
     scores.gx<- modeloGrid@scores
     biplot.gx.u<- scores.gx[,1:2] * matrix(rep(eigenvalues.gx[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.gx.v<- loadings.gx[,1:2] * matrix(rep(eigenvalues.gx[1:2]^0.5, Nenv), ncol=2, byrow=T)

     coordgenotype=biplot.gx.u
     coordenviroment=biplot.gx.v
     eigenvalues = eigenvalues.gx[1:2]

     labelaxes <- paste("Component ",1:ncol(diag(eigenvalues.gx)), sep = "")
     labelaxes=labelaxes

     vartotal = round(as.numeric(sum(eigenvalues)),5)
     varexpl = c(round(as.numeric((eigenvalues[1]/vartotal) *100), 2),
                 round(as.numeric((eigenvalues[2]/vartotal) *100), 2))

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

     AMMI2<-AMMI1+geom_segment(xend=0,yend=0,col=alpha(colEnv,0.5),data=subset(plotdata,type=="environment"))+
       geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)
     if(titles==TRUE){AMMI2<-AMMI2+ggtitle("gAMMI Biplot")}
     }
   if (type == "lAMMI"){             # L-AMMI model
     modeloLocantore<- PcaLocantore(rresiduals.x)
     eigenvalues.lx<- modeloLocantore@eigenvalues
     loadings.lx<- modeloLocantore@loadings
     scores.lx<- modeloLocantore@scores
     biplot.lx.u<- matrix(c(-scores.lx[,1], scores.lx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues.lx[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.lx.v<- matrix(c(-loadings.lx[,1], loadings.lx[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues.lx[1:2]^0.5, Nenv), ncol=2, byrow=T)
     coordgenotype=biplot.lx.u
     coordenviroment=biplot.lx.v
     eigenvalues = eigenvalues.lx

     labelaxes <- paste("Component ",1:ncol(diag(eigenvalues.lx)), sep = "")
     labelaxes=labelaxes

     vartotal = round(as.numeric(sum(eigenvalues)),5)
     varexpl = c(round(as.numeric((eigenvalues[1]/vartotal) *100), 2),
                 round(as.numeric((eigenvalues[2]/vartotal) *100), 2))

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

     AMMI2<-AMMI1+geom_segment(xend=0,yend=0,col=alpha(colEnv,0.5),data=subset(plotdata,type=="environment"))+
       geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)

     if(titles==TRUE){AMMI2<-AMMI2+ggtitle("lAMMI Biplot")}

     }
   if (type == "ppAMMI"){             # PP-AMMI model
     modeloProj<- PcaProj(rresiduals.x)
     eigenvalues.px<- modeloProj@eigenvalues
     loadings.px<- modeloProj@loadings
     scores.px<- modeloProj@scores
     biplot.px.u<- matrix(c(-scores.px[,1], -scores.px[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues.px[1:2]^0.5, Ngen), ncol=2, byrow=T)
     biplot.px.v<- matrix(c(-loadings.px[,1], -loadings.px[,2]), ncol=2, byrow=F) * matrix(rep(eigenvalues.px[1:2]^0.5, Nenv), ncol=2, byrow=T)

     coordgenotype=biplot.px.u
     coordenviroment=biplot.px.v
     eigenvalues = eigenvalues.px

     labelaxes <- paste("Component ",1:ncol(diag(eigenvalues.px)), sep = "")
     labelaxes=labelaxes[1:2]

     vartotal = round(as.numeric(sum(eigenvalues)),5)
     varexpl = c(round(as.numeric((eigenvalues[1]/vartotal) *100), 2),
                 round(as.numeric((eigenvalues[2]/vartotal) *100), 2))

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

     AMMI2<-AMMI1+geom_segment(xend=0,yend=0,col=alpha(colEnv,0.5),data=subset(plotdata,type=="environment"))+
       geom_text(aes(col=type,label=label,size=type),show.legend = FALSE)

     if(titles==TRUE){AMMI2<-AMMI2+ggtitle("ppAMMI Biplot")}
     }


  if(footnote==TRUE){
    footnotetxt=paste("\n Biplot showing components 1 and 2",sum(varexpl),"% of the total variation")

    AMMI2<-AMMI2+ labs(caption = footnotetxt)+theme(plot.caption = element_text(size=10,hjust=0,face="italic"))
  }



  return(AMMI2)
}
