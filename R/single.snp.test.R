single.snp.test <-
function ( snps, trait, adj.var=NULL ,
   type = c("gaussian", "binomial", "families","survival","casecohort") ,
   famid , patid , fid , mid ,
   start.time , stop.time , subcohort , stratvar=NA , robust=F ,
   marker.label=NA ,
   prt=T  ) {

    type <- match.arg ( type )
    
    switch ( type,
      gaussian =  test <- "F" ,
      binomial =  test <- "Chisq" ,
      families   =  test <- "wTDT"  ,
      survival = { print("not ready!") ; stop("end") } ,
      casecohort = test <- "case.cohort"
    )

    snps <- as.matrix(snps)
    N <- dim(snps)[1]
    ns <- dim(snps)[2]
    
    if (!all(is.na(marker.label)) ) {
      if ( ns != length(marker.label) ) {
        stop("Error in single.snp.test: marker.label don't matched nummber of SNPs.")
      }
    }
    
    if (length(trait) != N) {
        stop("length of trait does not match dimension of snps")
    }

    adjusted <- FALSE
    if (!all(is.null(adj.var))) { adjusted <- TRUE }

    # handle missing values

    miss.value <- which(is.na(trait))
    if (adjusted) {
        adj.var <- as.matrix(adj.var)
        miss.value <- unique(c(miss.value, which(apply(is.na(adj.var), 1, any))))
    }
    if ( length(miss.value) > 0) {
        if (adjusted) {
          adj.var <- adj.var[-miss.value, , drop = FALSE]
         }
        trait <- as.numeric(trait[-miss.value])
        snps  <- snps[-miss.value, ,drop=F]
        N <- dim(snps)[1]
    }

    # Adjustment

    if (adjusted) {
        adj.var <- as.matrix(adj.var)
        if (dim(adj.var)[1] != N)
            stop("Error in stepwise: length of patid does not match number of rows in adj.cov")
    }

  ######################
  # print

  if ( prt==T ) {
    cat ("Start procedure: single.snp.test\n")
    cat ("Individuals:          ",N, " (" , length(miss.value) ," excluded)","\n",sep="")
    cat ("SNPs:                 ",ns,"\n",sep="")
    cat ("Trait type:           ",type,"\n",sep="")
    cat ("Statistic:            ",test,"\n"  ,sep="")
    if (  is.null(dim(adj.var)) ) {
      cat ("Number of covariates: 0\n\n",sep="")
    } else {
      cat ("Number of covariates: ",dim(adj.var)[2],"\n\n" ,sep="")
    }
  }

  ######################
  # single SNP analysis
  
  switch ( type,
           gaussian   = res <- single.snp.test.gaussian ( snps, trait, adj.var=adj.var , prt=prt  ) ,
           binomial   = res <- single.snp.test.binomial ( snps, trait, adj.var=adj.var , prt=prt  ) ,
           families   = { res <- single.snp.test.families ( snps, trait, adj.var=adj.var,
                            famid , patid , fid , mid , prt=prt  ) } ,
           survival   = { print("not supported at the moment!") ; stop("end") } ,
           casecohort = res <- single.snp.test.casecohort ( snps, trait ,
                           patid , start.time , stop.time , subcohort , stratvar , robust=robust ,
                           adj.var=adj.var , prt=prt  ) )

 if (!all(is.na(marker.label))) {
   res[,"SNP"] <- marker.label
  }
 return ( res )

} # end of single.snp.test

