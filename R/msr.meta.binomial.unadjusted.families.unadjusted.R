msr.meta.binomial.unadjusted.families.unadjusted <-
function (
    data1 , control.options1 , data2 , control.options2 ,
    global.options=list(maxSNP=3 , lim=.05 , meta.method="fisher" ) )
{

  lenx <- length(data2$famid);
  lest <- paste(rep(" ",1000),collapse="");
  pr   <- rep(lest,1);


  maxSNP <- global.options$maxSNP
  nloc <- ncol(data1$genotypes)

  res.list <- list()
  pval <- matrix ( NA , nr=choose(nloc,2),ncol=5 )
  pos <- matrix (rep(0,choose(nloc,2)*2) ,ncol=2);

  ### pair begin ###
  k <- 1
cat(paste("Iteration with 2 SNPs at same time.   System.time = ", Sys.time(),"\n",sep=""))
cat(paste("Number of SNP pairs = ", choose(nloc,2),"\n\n",sep=""))

  for( i in 1:(nloc-1)) {
    for( j in (i+1):nloc) {

if ( (k%%5000)==0 ) { cat(paste("Step =  ",k,"   System.time = ", (Sys.time()), "\n",sep="")) }

      # data set 1 -------------------------------------------------------------
    
      xgeno <- data1$genotypes[,c(i,j)];
      pv1 <-  msr.binomial.haplotype.test.unadjusted ( xgeno , data1$trait , lim =global.options$lim,
                  baseline.hap=control.options1$baseline.hap ,
                  min.count=control.options1$min.count )$global.p.value

      # data set 2 -------------------------------------------------------------
      
      geno<-apply(data2$genotypes[,c(i,j)],1,paste,collapse="");
      zzz<- .C("haptdpn",as.character(data2$famid),as.character(data2$patid),
          as.character(geno),
          as.integer(data2$trait),as.integer(lenx),as.double(global.options$lim),pvres=pr)[[7]];
      x<-strsplit(zzz," ");
      x[[1]]<-x[[1]][x[[1]]!=""];
      x5<-as.numeric(x[[1]][4]);
      x1<-as.numeric(x[[1]][1]);
      if(x5==1){
        pv2 <- 1.0-pchisq(x5,1);
      } else{
        pv2 <- 1.0-pchisq((x5-1)*x1/x5,x5-1);
      }

      # save result ------------------------------------------------------------
      
      pval[k,1] <- pv1
      pval[k,2] <- pv2
      pos[k,]<-c(i,j);
      k <- k+1

    } # end of for j
  } # end of for i
   
   pval[,3] <- meta.func ( pval[,1] , pval[,2] , global.options$meta.method )   
   pval[,4] <- rank (pval[,1])
   pval[,5] <- rank (pval[,2])
   ind <- (order(pval[,3] ))[1:min(global.options$nt,nrow(pval))  ]
   pval <- pval[ind,,drop=F]
   pos <- pos [ind,,drop=F]
   res.list[[2]] <- as.matrix(cbind(pos,pval))
   rownames(res.list[[2]]) <- NULL
   colnames(res.list[[2]]) <- c( paste("SNP",1:2,sep=""),"pval.data1","pval.data2","meta.value",
       "rank.pval.data1","rank.pval.data2")

  i<-3

  while ( (i<=nloc) && (i<=maxSNP) )  {
cat(paste("Iteration with ",i," SNPs at same time.   System.time = ", (Sys.time()),"\n",sep=""))

   #pval <- matrix ( NA , nr=choose(nloc,2),ncol=5 )

        # create data set with all possible SNP combinations
        BestPos <- as.matrix(res.list[[i - 1]] [, 1:(i - 1),drop=F])
        storage.mode(BestPos) <- "integer"
        newdim <- as.integer(c(dim(BestPos)[1]*(nloc-2),dim(BestPos)[2]+1))
        out <- .C("create_pattern_matrix", pattern=as.integer(BestPos) , ndim=dim(BestPos)  ,
                               snps=as.integer(1:nloc) , snplen=nloc ,
                               newpat=as.integer(rep(0,newdim[1]*newdim[2])) ,
                               newpatdim=newdim,len=as.integer(0) )
        BestPos <- (matrix(out$newpat,nr=newdim[1] ,nc=newdim[2],byrow=F))[1:out$len,,drop=F]
        cat(paste("Iteration with ",i," SNPs at same time. ---- ",
            dim(BestPos)[1], " detected SNP combinations -----\n" ,sep=""))

       pval <- NULL

    for ( k in 1:(dim(BestPos)[1]) ) {
    
if ( (k%%5000)==0 ) { cat(paste("Step =  ",k,"   System.time = ", (Sys.time()), "\n",sep="")) }
    
      Pos<-sort(BestPos[k,])

      # data set 1 -------------------------------------------------------------
             
        xgeno <- data1$genotypes[,Pos];
        pv1 <-  msr.binomial.haplotype.test.unadjusted ( xgeno , trait=data1$trait ,
                  lim =global.options$lim,
                  baseline.hap=control.options1$baseline.hap ,
                  min.count=control.options1$min.count )$global.p.value
      # data set 2 -------------------------------------------------------------

      geno<-apply(data2$genotypes[,Pos],1,paste,collapse="");
      zzz<- .C("haptdpn",as.character(data2$famid),as.character(data2$patid),
          as.character(geno),as.integer(data2$trait),as.integer(lenx),
          as.double(global.options$lim),pvres=pr)[[7]];
      x<-strsplit(zzz," ");
      x[[1]]<-x[[1]][x[[1]]!=""];
      x5<-as.numeric(x[[1]][4]);
      x1<-as.numeric(x[[1]][1]);
      if(x5==1){
        pv2 <- 1.0-pchisq(x5,1);
      } else{
        pv2 <- 1.0-pchisq((x5-1)*x1/x5,x5-1);
      }

      # save result ------------------------------------------------------------

      #pval <- rbind( pval, c(pv1 , pv2 ,  meta.fisher ( c(pv1,pv2) ) ))
      pval <- rbind( pval, c(pv1 , pv2 ,  NA , NA , NA ))  # via ranks
      

    } # end of for k
                                                                                
   pval[,3] <- meta.func ( pval[,1] , pval[,2] , global.options$meta.method )      
   pval[,4] <- rank (pval[,1])
   pval[,5] <- rank (pval[,2])
   ind <- (order(pval[,3]))[1:min(global.options$nt,dim(pval)[1]) ]
   pval <- pval[ind,,drop=F]
   BestPos <- BestPos [ind,,drop=F]

   res.list[[i]] <- as.matrix(cbind(BestPos,pval))

   rownames(res.list[[i]]) <- NULL
   colnames(res.list[[i]]) <- c( paste("SNP",1:i,sep=""),"pval.data1","pval.data2","meta.value",
       "rank.pval.data1","rank.pval.data2")
   i<-i+1;
  }  # while

  return (res.list)
} # end of msr.meta.binomial.unadjusted.families.unadjusted

