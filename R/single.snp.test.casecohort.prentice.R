single.snp.test.casecohort.prentice <-
function ( snps, trait ,
   patid , start.time , stop.time ,  subcohort , stratvar=NA ,
   robust , adj.var , prt  ) {

  eps        <- 0.00005
  patid      <- as.character(patid)
  start.time <- as.numeric(start.time)
  stop.time  <- as.numeric(stop.time)
  stratvar   <- as.numeric(stratvar)

  N <- nrow(snps)
  
  # Prentice:
  start.time <- ifelse( (trait==1)&(subcohort==0) ,  stop.time-eps , start.time )

  
  ysurv <- Surv ( start.time , stop.time , trait )

  nloci <- dim(snps)[2]
  
  res <- as.data.frame(matrix(0,nrow=nloci , ncol=9), stringsAsFactors = F)

  if ( robust==F ) {
    colnames(res) <- c("SNP","N","type","beta","se(beta)", "exp(beta)", "lower.95" , "upper.95" ,"p.value")
  } else {
    colnames(res) <- c("SNP","N","type","beta","robust se (beta)", "exp(beta)", "lower.95" , "upper.95" ,"p.value")
  }

 if (!all(is.na(stratvar))) {

    for ( j in 1:dim(snps)[2] ) {
      x <- as.matrix(as.numeric(alleleRto1 ( snps[,j] )),ncol=1)
      phfit <- summary(coxph(  ysurv  ~  cbind ( x , adj.var ) + cluster(as.factor(patid))+ strata(stratvar) ,
          method="breslow" , robust=robust  ))  # marker and covariates
      if ( robust==F ) {

      res[j,3]  <- "case.cohort.prentice"
      res[j,c(1:2,4:9)] <- c(j,phfit$n,
           as.numeric( formatC( c((phfit$coefficients)[1,c("coef","se(coef)")] ,
           (phfit$conf.int)[1,c("exp(coef)","lower .95","upper .95")] ,
           (phfit$coefficients)[1,c("Pr(>|z|)")] )   ,format="fg"))  )

      } else{
      res[j,3]  <- "case.cohort.prentice"
      res[j,c(1:2,4:9)] <- c(j,phfit$n,
           as.numeric(formatC( c((phfit$coefficients)[1,c("coef","robust se")] ,
           (phfit$conf.int)[1,c("exp(coef)","lower .95","upper .95")] ,
           (phfit$coefficients)[1,c("Pr(>|z|)")] )   ,format="fg"))   )

      }
    }
  } else {
    # no stratvar
  for ( j in 1:dim(snps)[2] ) {
      x <- as.matrix(as.numeric(alleleRto1 ( snps[,j] )),ncol=1)
      phfit <- summary(coxph(  ysurv  ~  cbind ( x , adj.var ) + cluster(as.factor(patid)) ,
          method="breslow" , robust=robust  ))  # marker and covariates
      if ( robust==F ) {

      res[j,3]  <- "case.cohort.prentice"
      res[j,c(1:2,4:9)] <- c(j,phfit$n,
           as.numeric( formatC( c((phfit$coefficients)[1,c("coef","se(coef)")] ,
           (phfit$conf.int)[1,c("exp(coef)","lower .95","upper .95")] ,
           (phfit$coefficients)[1,c("Pr(>|z|)")] )   ,format="fg"))  )

      } else{
      res[j,3]  <- "case.cohort.prentice"
      res[j,c(1:2,4:9)] <- c(j,phfit$n,
           as.numeric(formatC( c((phfit$coefficients)[1,c("coef","robust se")] ,
           (phfit$conf.int)[1,c("exp(coef)","lower .95","upper .95")] ,
           (phfit$coefficients)[1,c("Pr(>|z|)")] )   ,format="fg"))   )

      }
    }
  
  
  }

  return ( res )
  
} # end of single.snp.test.casecohort.prentice

