msr <-
function (  famid , patid , fid , mid , snps, trait,
    adj.var=NA , lim =0.05, maxSNP = 3 ,
    nt = 10, pair.begin = FALSE, pattern.begin.mat=NA,
    type = "gaussian" ,
    baseline.hap="max" , min.count=10 , sort=F )
{

    # init input

    snps <- as.matrix(snps)

    N <- dim(snps)[1]
    ns <- dim(snps)[2]

    if (maxSNP > 15)
        stop("maxSNP must smaller than 16")
    if ( maxSNP > ns )
        stop("maxSNP must be smaller or equal to number of SNPs.")
    if (!all(trait[!is.na(trait)] == 0 | trait[!is.na(trait)] == 1) &&
        (type == "binomial"))
        stop("trait should be 0 for controls or 1 for cases")
    if (length(trait) != N)
        stop("number of individuals does not match length of trait")
    if (ns%%dim(snps)[2] != 0)
        stop("odd number of coulmns in snps")
    if ( !(baseline.hap=="max") | baseline.hap=="min" )
        stop ("baseline haplotype not correct specified.")

    # check input test

    trait.type <- match(type, c("binomial","gaussian","families"))

    switch ( type,
      gaussian =  test <- "F" ,
      binomial =  test <- "Chisq" ,
      families = test <- "wTDT" )
      
    # adjustment
      
    adjusted <- FALSE
    if (!all(is.na(adj.var))) {
      adjusted <- TRUE
      adj.var <- as.matrix(adj.var)
    }
      
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
        snps  <- snps[-miss.value, ]

    }

    if ( (trait.type == "families") & (sort==T) ) {

       famorder <- order.families (famid,patid,fid,mid)
       famid    <- famid[famorder]
       patid    <- patid[famorder]
       fid      <- fid[famorder]
       mid      <- mid[famorder]
       snps     <- snps[famorder,]
       trait    <- trait[famorder]

       if ( adjusted ) {
           adj.var <- as.matrix(adj.var)[famorder,]
       }
    }

    # message input
    cat ("Start procedure: stepwise.\n")
    cat ("Individuals:           ",N, " (", N-dim(snps)[1] , " individual(s) exluded)\n",sep="")
    cat ("SNPs:                  ",ns,"\n",sep="")
    cat ("Trait type:            ",type,ifelse (adjusted," (adjusted)" , " (unadjusted)") ,"\n" ,sep="")
    cat ("Statistic:             ",test,"\n" ,sep="" )
    cat ("Max SNPs:              ",maxSNP,"\n",sep=""  )
    cat ("Number best pattern    ",nt,"\n",sep="" )
    cat ("Threshold (lim):       ",lim,"\n",sep="" )

    if (  is.null(dim(adj.var)) ) {
      cat ("Number of covariates:  0\n\n")
    } else {
      cat ("Number of covariates:  ",dim(adj.var)[2],"\n\n",sep="" )
    }

    N <- dim(snps)[1]

    # Adjustment

    if (adjusted) {
        adj.var <- as.matrix(adj.var)
        if (dim(adj.var)[1] != N)
            stop("length of number of individuals does not match number of rows in adj.cov")
    }

################################################################################
#                                                                              #
#  Case control data >> family="binomial"                                      #
#                                                                              #
################################################################################

  if ( type=="binomial" ) {
       # unadjusted
       if ( !adjusted ) {
         results <- msr.binomial.forward.unadjusted ( snps, trait, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat,
            baseline.hap , min.count )
       } else {
         # adjusted
         results <- msr.binomial.forward.adjusted ( snps, trait,  adj.var, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat,
            baseline.hap , min.count )
       }
  }


################################################################################
#                                                                              #
#  Quantitative trait loci (QTL) >> family="gaussian"                          #
#                                                                              #
################################################################################

  if ( type=="gaussian" ) {
       # unadjusted
       if ( !adjusted ) {
         results <- msr.gaussian.forward.unadjusted ( snps, trait, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat,
            baseline.hap , min.count )
       } else {

         # adjusted
         results <- msr.gaussian.forward.adjusted ( snps, trait,  adj.var, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat,
            baseline.hap , min.count )
       }
  } # gaussian

################################################################################
#                                                                              #
#  Family data weigthed TDT statistic (TDT) >> type="families"                 #
#                                                                              #
################################################################################

  if ( type=="families" ) {
    if ( !adjusted ) {
      # unadjusted
      results <- msr.families.unadjusted ( famid , patid , fid , mid ,
          trait , snps , pair.begin=pair.begin , lim = lim, maxSNP = maxSNP, nt = nt )
    } else {
      # adjusted
      stop ( "Error in msr: weighted TDT for adjustment sets is not avaiable." )
    }
  }  # families

   return ( results )
   
}
