single.snp.test.gaussian <-
function ( snps, trait, adj.var=NULL , prt=T  ) {

  snps <- as.matrix(snps)
  N <- dim(snps)[1]
  ns <- dim(snps)[2]

  adjusted <- FALSE
  if (!all(is.null(adj.var))) { adjusted <- TRUE }


    pval <- rep(-1,ns)
    nind <- rep(-1,ns)
    aic  <- rep(-1,ns)
    beta <- rep(-1,ns)
    beta.stderr <- rep(-1,ns)
    y <- trait
    
  for ( i in 1:ns ) {
    x <- as.matrix(as.numeric(alleleRto1 ( snps[,i] )),ncol=1)
    miss.value <- which(is.na(x))
    if ( length(miss.value)>0 ) {
      fit <- multi.snp.test ( y[-miss.value] , x[-miss.value,,drop=F] ,
          x.adj=adj.var[-miss.value] , type="gaussian" )
    } else {
      fit <- multi.snp.test ( y , x , x.adj=adj.var , type="gaussian" )
    }
    nind[i]        <- length((fit$fit.glm1)$residuals)
    pval.model     <- (fit$aov.glm)$`Pr(>F)`
    
    if ( all(is.na(pval.model)) ) {

      pval[i]        <- NA
      aic[i]         <- NA
      beta[i]        <- NA
      beta.stderr[i] <- NA

    } else {
    
      pval[i]        <- pval.model[!is.na(pval.model)]
      aic[i]         <- AIC(fit$fit.glm1)
      beta[i]        <- fit$beta[2,1]
      beta.stderr[i] <- fit$beta[2,2]

    }
    
  }

  res <- data.frame ( snp=1:ns , N=nind , type="gaussian",
      beta.estimate=beta, beta.stderr=beta.stderr,test="F",
      p.value=pval,aic=aic , stringsAsFactors=F)
  colnames(res) <- c("SNP","N","type","beta","se(beta)", "Test" ,"p.value","AIC")
  
  return ( res )
} # end of single.snp.test.gaussian

