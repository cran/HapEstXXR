msr.binomial.forward.adjusted <-
function (
    snps, trait, adj.var , lim = 0.05 , maxSNP = 3,
    nt = 10, pair.begin = FALSE, pattern.begin.mat=NA ,select.criteria = "p.value" ,
    baseline.hap="max" , min.count=10 )
{

    N <- dim(snps)[1]
    ns <- dim(snps)[2]
    nt <- as.integer(nt)

    ## Begin stepwise regression

    res <- list(NA)
    
    if ( any(!is.na( pattern.begin.mat )) ) {
      if ( pair.begin==TRUE ) {
        stop ( paste ( "choose either",
            " begin with defined pattern or pair.wise." , sep="" ) )
      }
      if ( (dim(pattern.begin.mat)[2] >= maxSNP) ||
           (dim(pattern.begin.mat)[2] >= dim(snps)[2] ) )
        stop ( "dimension of begin.pattern.mat not adequate." )
    }
    
    
    if ( pair.begin==FALSE ) {

    ############################################################################
    # start with all two pair haplotypes!
    # begin single SNP test

    cat(paste("Iteration with 1 SNP  at same time.   System.time = ", Sys.time(),"\n",sep=""))

    nind <- df <- pval <- rep (NA,dim(snps)[2])

    single.test <- single.snp.test ( snps , trait , prt=F , type="binomial"  )

            nind  <- as.integer(single.test$N)
            df    <- as.integer(rep(1,length(single.test$N)))
            pval  <- as.numeric(single.test$p.value)

    i <- 1

      # save results from single SNP or pairwise tests

        if ( nt<=length(pval) ) {
          ii <- (order(as.numeric(pval)))[1:nt]
        } else {
          ii <- (order(as.numeric(pval)))
        }
      res[[1]] <- data.frame (as.integer(single.test$SNP[ii]),
          rep("binomial",length(ii)) ,
          as.integer(nind[ii]) , as.integer(df[ii]) ,
          as.numeric(pval[ii]) , stringsAsFactors=F  )
      colnames(res[[1]]) <- c(paste("snp", 1:i, sep = "") , "type" , "nSubj" , "df" , "p.value" )
      rownames(res[[1]]) <- 1:length(ii)
    } else {
    
    ############################################################################
    # start with all two pair haplotypes!    
    # begin pair wise
    
         cat(paste("Iteration with 2 SNPs at same time.   System.time = ", Sys.time(),"\n",sep=""))

        # construct all pairs

        Z <- 1:ns
        X <- rep(Z, rep.int(length(Z), length(Z)))
        Y <- rep(Z, times = ceiling(length(X)/length(Z)))
        cont <- ifelse(Y > X, T, F)
        snp.pos <- cbind(X[cont], Y[cont])
        rm(X, Y, Z)

        nind <- df <- pval <- rep (NA,dim(snp.pos)[1])

        # Analysis of all pairs
cat(paste("Number of SNP pairs = ", dim(snp.pos)[1],"\n\n",sep=""))

        for (j in 1:dim(snp.pos)[1]) {

if ( (j%%5000)==0 ) { cat(paste("Step =  ",j,"   System.time = ", (Sys.time()), "\n",sep="")) }
          geno.pair <- snps[, snp.pos [j, ], drop = FALSE]

          hap.test <- msr.binomial.haplotype.test.adjusted ( geno.pair, trait, adj.var , lim =lim,
              baseline.hap=baseline.hap , min.count=min.count )

         # found no haplotypes with probability over lim

         if ( all( is.na( hap.test$haplotypes) ) ) {

            cat ( "Step " ,
                ": SNPs ", snp.pos [j, ] ," all inferred haplotypes with probability " ,
                "below threshold (lim = ",lim,")\n" , sep=""  )

            nind[j] <- NA
            df[j]   <- NA
            pval[j] <- NA

         } else {

            nind[j]  <- hap.test$nSubj
            df[j]    <- hap.test$df
            pval[j]  <- hap.test$global.p.value

         }
         
       }

       i <- 2

      # save results from single SNP or pairwise tests

        if ( nt<=length(pval) ) {
          ii <- (order(as.numeric(pval)))[1:nt]
        } else {
          ii <- (order(as.numeric(pval)))
        }

      res[[i]] <- data.frame (snp.pos[ii, ,drop=FALSE], "binomial" ,
      nind[ii] , df[ii] , pval[ii] , stringsAsFactors=F , row.names=1:length(ii) )
      colnames(res[[i]]) <- c(paste("snp", 1:i, sep = "") , "type" , "nSubj" , "df" , "p.value" )

      }

    #################################################
    # Next iteration: from 3,4 ... (or colnumber of pattern.begin.mat) to maxSNP

    i <- i+1

    # i No. of SNPs for every haplotype pattern
    while ( i <= maxSNP ) {

        cat(paste("Iteration with ",i," SNPs at same time.   System.time = ", (Sys.time()),"\n",sep=""))
        # create matrix with all possible combinations
        snp.pos <- as.matrix(res[[i - 1]] [, 1:(i - 1),drop=F])

        storage.mode(snp.pos) <- "integer"
        newdim <- as.integer(c(dim(snp.pos)[1]*(ns-2),dim(snp.pos)[2]+1))
        out <- .C("create_pattern_matrix", pattern=as.integer(snp.pos) , ndim=dim(snp.pos)  ,
                               snps=as.integer(1:ns) , snplen=ns ,
                               newpat=as.integer(rep(0,newdim[1]*newdim[2])) ,
                               newpatdim=newdim,len=as.integer(0) )
        snp.pos <- (matrix(out$newpat,nr=newdim[1] ,nc=newdim[2],byrow=F))[1:out$len,,drop=F]
        cat(paste("Iteration with ",i," SNPs at same time. ---- ", 
            dim(snp.pos)[1], " detected SNP combinations -----\n" ,sep=""))

        nind <- df <- pval <- rep(NA,dim(snp.pos)[1])
        
        # evaluate all combinations in the step before
        k <- 0
        for (j in 1:(dim(snp.pos)[1])) {
        
            Pos <- as.integer(snp.pos[j,])
            if ( (j%%5000)==0 ) { cat(paste("Step =  ",j,"   System.time = ", (Sys.time()), "\n",sep="")) }
            geno <- matrix(snps[, Pos], N,i)

            hap.test <- msr.binomial.haplotype.test.adjusted ( geno,trait, adj.var , lim =lim,
                baseline.hap=baseline.hap , min.count=min.count )

             # found no haplotypes with probability over lim

             if ( all( is.na( hap.test$haplotypes) ) ) {
                txt.pos <- paste ( Pos, collapse=" ")
                cat ( "Step " , i ,
                ": SNPs ", txt.pos ," all inferred haplotypes with probability " ,
                "below threshold (lim = ",lim,")\n" , sep=""  )
             } else {               
               nind[j] <- hap.test$nSubj 
               df[j]   <- hap.test$df
               pval[j] <- hap.test$global.p.value
             }
        } # for j

        # save results from single SNP tests of pairwise test

          if ( nt<=length(pval) ) {
            ii <- (order(as.numeric(pval)))[1:nt]
          } else {
            ii <- (order(as.numeric(pval)))
          }

        if ( length(ii)<1 ) {

          stop ( paste ( "in stepwise: abort of the programm, because find no ",
             "haplotype distributions with frequent haplotypes over lim=",lim ,
             sep="" ) )
          return(res)

        }
        
        res[[i]] <- data.frame (snp.pos[ii, ,drop=FALSE], "binomial" ,
            nind[ii] , df[ii] , pval[ii] , stringsAsFactors=F , row.names =1:length(ii) )
        colnames(res[[i]]) <- c(paste("snp", 1:i, sep = "") , "type" , "nSubj" , "df" , "p.value" )

        i <- i + 1
    } # while i

    return(res)
    
} ## end of stepwise.forward ###################################################

