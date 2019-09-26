#' Weighted Gabriel Eigen imputation
#'
#' @param DBmiss a data frame or matrix that contains the genotypes in the rows and the
#'    environments in the columns when there are no replications of the experiment.
#' @param Winf peso inferior
#' @param Wsup peso superior
#'
#' @references FALTA
#' @return Imputed data with Weighted Gabriel Eigein method
#' @keywords internal
#' @importFrom MASS ginv
#'

WGabriel<-function(DBmiss,Winf,Wsup){

  if (missing(DBmiss)) stop("Need to provide DBmiss data frame or matrix")

  # stopifnot(
  #   class(DBmiss) %in% c("matrix", "data.frame"),
  #   class(Winf) == "numerical",
  #   class(Wsup) == "numerical"
  # )


  ######
  ######Funci?n Gabriel Pesado
  ######
  Gabriel.Pesado<-function(X,w){
    nfilasX<-nrow(X)
    ncolX<-ncol(X)
    #nfilasX
    #ncolX
    if(ncolX>nfilasX) {X<-t(X)}
    ######CREANDO UNA MATRIZ A PARTIR DE KRZANOWSKI(1988)
    ######OJO QUE ANTES DE CORRER ESTE PROGRAMA TOCA CAMBIAR EL N?MERO DE COLUMNAS SI SE CAMBIA EL DATASET
    #X<-matrix(scan('C:/SAA/Gabriel(2002)_Pesquisa/Krzanowski(1988).txt'), ncol=5, byrow=T) # Pasa el conjunto de datos a una matriz
    #X # Imprime la matriz de datos
    XMISSING<-X
    is.matrix(X) # Pregunta si los datos que se tiene son tipo matriz
    totalelementos<-nfilasX*ncolX

    ######CREACI?N DE MATRIZ INDICADORA, C?LCULO DEL TOTAL DE DATOS FALTANTES Y C?LCULO INICIAL DE LA MEDIA A PARTIR DE MATRIZ INCOMPLETA
    indicamissing<-is.na(XMISSING)# indica la ocurrencia de faltantes en la matriz
    indicadora<-indicamissing*1
    totalfaltantes<-sum(indicadora)
    medcol<-t(matrix(colMeans(XMISSING, na.rm=TRUE))) # Calcula media por columnas sin tener en cuenta datos faltantes


    ######CREACI?N DE MATRIZ CON LAS POSICIONES DE LOS DATOS FALTANTES
    ######REEMPLAZANDO DATOS FALTANTES CON LA MEDIA COLUMNA INICIAL
    posicfaltantes<-matrix(0,totalfaltantes,2)
    indi3<-0
    for (indi1 in 1:nfilasX){
      for(indi2 in 1:ncolX){
        if (indicadora[indi1,indi2]==1){
          indi3<-indi3+1
          posicfaltantes[indi3,1]=indi1
          posicfaltantes[indi3,2]=indi2
          XMISSING[indi1,indi2]=medcol[1,indi2]
        } #FIN DEL if (indicadora[indi1,indi2]
      }# FIN DEL for(indi2 in 1:ncolX)
    } # FIN DEL for (indi1 in 1:nfilasX)
    #posicfaltantes
    #XMISSING

    ######DECLARACI?N DE ELEMENTOS PREVIOS PARA EL ESQUEMA ITERATIVO
    completed_matrix<-XMISSING
    maxiter<-1000
    iter<-0
    desvest<-sqrt(t(diag(var(XMISSING))))
    #print(desvest)
    Xmissingestand<-scale(XMISSING)
    Prediction_Gabriel<-matrix(0,nfilasX,ncolX)
    epsilon<-1*(10**(-6))
    stabilitycrit<-1
    stabilityiter<-1


    ######ESQUEMA ITERATIVO SEG?N ARCINIEGAS ET AL (2010)

    while (stabilitycrit>epsilon & stabilityiter<=30){


      for (i1 in 1:totalfaltantes){

        posic_linha<-posicfaltantes[i1,1]
        posic_col<-posicfaltantes[i1,2]
        auxiliar1<-rbind(Xmissingestand[posic_linha,],Xmissingestand[-posic_linha,])
        auxiliar2<-cbind(auxiliar1[,posic_col],auxiliar1[,-posic_col])
        #print(i1)
        #print(auxiliar2)
        Y1puntoTRANS<-auxiliar2[1,2:ncolX]
        Ypunto1<-auxiliar2[2:nfilasX,1]
        Ybarra11<-auxiliar2[-1,-1]

        ######Criterio para seleccionar el n?mero de componentes en la imputaci?n
        Ybarra11_Estand<-scale(Ybarra11)
        dvs_Ybarra11_Estand<-svd(Ybarra11_Estand)
        autovalores<-matrix(dvs_Ybarra11_Estand$d)**2
        sum_autovalores<-sum(autovalores)
        CritCalinski1<-(autovalores/sum_autovalores)
        CritCalinski2<-matrix(cumsum(CritCalinski1))
        CritCalinski3<-matrix(1,ncol(Ybarra11_Estand),1)
        for (Calins in 1:ncol(Ybarra11_Estand)){
          if (CritCalinski2[Calins,1]>0.70){CritCalinski3[Calins,1]=0}
        } # FIN de for (Calins in 1:ncol(Ybarra11_Estand)){
        f<-sum(CritCalinski3)+1
        #print(i1)
        #print(f)
        ####### FIN de Criterio para seleccionar el n?mero de componentes en la imputaci?n


        #print(Y1puntoTRANS)
        #print(Ypunto1)
        #print(Ybarra11)
        dvs.Ybarra11<-svd(Ybarra11)
        u_Ybarra11<-dvs.Ybarra11$u
        v_Ybarra11<-dvs.Ybarra11$v
        d_Ybarra11<-diag(dvs.Ybarra11$d)

        pred_estand<-w*Y1puntoTRANS%*%v_Ybarra11[,1:f]%*%MASS::ginv(d_Ybarra11[1:f,1:f])%*%t(u_Ybarra11[,1:f])%*%Ypunto1
        completed_matrix[posic_linha,posic_col]<-medcol[1,posic_col]+(pred_estand*desvest[1,posic_col])



      } ###FIN de for (i1 in 1:totalfaltantes){

      stabilitycrit<-(sum((completed_matrix-XMISSING)**2))/(sum((XMISSING*indicadora)**2))
      stabilityiter<-stabilityiter+1
      XMISSING<-completed_matrix
      Xmissingestand<-scale(completed_matrix)
      medcol<-t(matrix(colMeans(completed_matrix)))
      desvest<-sqrt(t(diag(var(completed_matrix))))

      #print(stabilityiter)
      #print(stabilitycrit)
      ######print (f)
      #print(completed_matrix)



    }### FIN de while (stabilitycrit>epsilon & stabilityiter<=500){

    list(Peso=w,NumeroIterGabriel=stabilityiter,CritConvergGabriel=stabilitycrit,GabrielImput=completed_matrix)
  }### Fin function Gabriel.Pesado



  ###
  ### Estudio de ponderaciones para Gabriel en Denis-Baril(1992) subconjunto completo
  ###


  nfilasDBmiss<-nrow(DBmiss)
  ncolDBmiss<-ncol(DBmiss)
  indicaDBmissing<-is.na(DBmiss)# indica la ocurrencia de faltantes en la matriz
  indicaDBcomplete<-1-indicaDBmissing
  posicobservadosDBmiss<-which(indicaDBmissing != 1, TRUE) #Matriz de posiciones con observados
  posicfaltantesDBmiss<-which(indicaDBmissing != 0, TRUE)  #Matriz de posiciones con faltantes
  TotalFaltantesDBmiss<-nrow(posicfaltantesDBmiss)
  TotalObservadosDBmiss<-nrow(posicobservadosDBmiss)
  MatrizPesos<-matrix(seq(Winf,Wsup,0.1))
  NumTotdePesos<-nrow(MatrizPesos)
  VetorResultDBmiss1<-matrix(0,1,2)

  for (w in 1:NumTotdePesos){
    EstimaObsDBmiss<-matrix(0,nfilasDBmiss,ncolDBmiss)
    for (iobsDB in 1:TotalObservadosDBmiss){
      CVDB<-DBmiss
      CVDB[posicobservadosDBmiss[iobsDB,1],posicobservadosDBmiss[iobsDB,2]]<-NA
      CVDB.Imputa<-Gabriel.Pesado(CVDB,MatrizPesos[w,1])
      EstimaObsDBmiss[posicobservadosDBmiss[iobsDB,1],posicobservadosDBmiss[iobsDB,2]]<-CVDB.Imputa$GabrielImput[posicobservadosDBmiss[iobsDB,1],posicobservadosDBmiss[iobsDB,2]]
      #print(iobsDB)
      #tempo1<-sum(is.na(CVDB))
      #print (tempo1)
      #print(CVDB)
      #print(CVKR.Imp)
      #print(MatrizValidacion)
    } # FIN de for (iobsDB in 1:TotalObservadosDBmiss){
    DBcomplete<-DBmiss
    DBcomplete[is.na(DBcomplete)] <- 0
    RMSPD.Obs<-sqrt((sum((EstimaObsDBmiss-(DBcomplete*indicaDBcomplete))^2))/TotalObservadosDBmiss)

    VetorResultDBmiss<-matrix(0,1,2)
    VetorResultDBmiss[1,1]<-MatrizPesos[w,1]
    VetorResultDBmiss[1,2]<-RMSPD.Obs
    VetorResultDBmiss1<-rbind(VetorResultDBmiss1,VetorResultDBmiss)

    #print(paste("Num.PesoExt",w), sep=" ")
    print(paste("  PesoExterno = ",MatrizPesos[w,1], sep= " "))

  } # FIN de for (w in 1:NumTotdePesos){

  ResultadosPesosDB<-VetorResultDBmiss1[-1,]
  colnames(ResultadosPesosDB)<-c("PesoW","RMSPD.Obs")
  ###ResultadosPesosDB
  MinimoRMSPD<-min(ResultadosPesosDB[,2])
  FilaconMinimoRMSPD<-subset(ResultadosPesosDB,ResultadosPesosDB[,2]==MinimoRMSPD)
  Wsugerido<-FilaconMinimoRMSPD[1,1]
  #WGabriel<-Gabriel.Pesado(DBmiss,Wsugerido)

  #list(Peso=WGabriel$Peso,NumeroIterWGabriel=WGabriel$NumeroIterGabriel,CritConvergWGabriel=WGabriel$CritConvergGabriel,GabrielWImput=WGabriel$GabrielImput)

  ###################
  ###################
  ###################
  Liminf2<-Wsugerido-0.1
  Limsup2<-Wsugerido+0.1
  MatrizPesos2<-matrix(seq(Liminf2,Limsup2,0.01))
  NumTotdePesos2<-nrow(MatrizPesos2)
  VetorResultDBmiss22<-matrix(0,1,2)

  for (w2 in 1:NumTotdePesos2){
    EstimaObsDBmiss2<-matrix(0,nfilasDBmiss,ncolDBmiss)
    for (iobsDB2 in 1:TotalObservadosDBmiss){
      CVDB2<-DBmiss
      CVDB2[posicobservadosDBmiss[iobsDB2,1],posicobservadosDBmiss[iobsDB2,2]]<-NA
      CVDB.Imputa2<-Gabriel.Pesado(CVDB2,MatrizPesos2[w2,1])
      EstimaObsDBmiss2[posicobservadosDBmiss[iobsDB2,1],posicobservadosDBmiss[iobsDB2,2]]<-CVDB.Imputa2$GabrielImput[posicobservadosDBmiss[iobsDB2,1],posicobservadosDBmiss[iobsDB2,2]]
      #print(iobsDB)
      #tempo1<-sum(is.na(CVDB))
      #print (tempo1)
      #print(CVDB)
      #print(CVKR.Imp)
      #print(MatrizValidacion)
    } # FIN de for (iobsDB in 1:TotalObservadosDBmiss){
    DBcomplete2<-DBmiss
    DBcomplete2[is.na(DBcomplete2)] <- 0
    RMSPD.Obs2<-sqrt((sum((EstimaObsDBmiss2-(DBcomplete2*indicaDBcomplete))^2))/TotalObservadosDBmiss)

    VetorResultDBmiss2<-matrix(0,1,2)
    VetorResultDBmiss2[1,1]<-MatrizPesos2[w2,1]
    VetorResultDBmiss2[1,2]<-RMSPD.Obs2
    VetorResultDBmiss22<-rbind(VetorResultDBmiss22,VetorResultDBmiss2)

    #print(w2)
    #print(MatrizPesos2[w2,1])
    print(paste("    PesoInterno = ",MatrizPesos2[w2,1], sep= " "))



  } # FIN de for (w2 in 1:NumTotdePesos){

  ResultadosPesosDB2<-VetorResultDBmiss22[-1,]
  colnames(ResultadosPesosDB2)<-c("PesoW","RMSPD.Obs")
  ###ResultadosPesosDB
  MinimoRMSPD2<-min(ResultadosPesosDB2[,2])
  FilaconMinimoRMSPD2<-subset(ResultadosPesosDB2,ResultadosPesosDB2[,2]==MinimoRMSPD2)
  Wsugerido2<-FilaconMinimoRMSPD2[1,1]
  WGabriel2<-Gabriel.Pesado(DBmiss,Wsugerido2)

  list(Peso=WGabriel2$Peso,NumeroIterWGabriel=WGabriel2$NumeroIterGabriel,CritConvergWGabriel=WGabriel2$CritConvergGabriel,GabrielWImput=WGabriel2$GabrielImput)

}### FIN function WGabriel<-function(DBmiss){
