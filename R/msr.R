msr <-
function (  famid , patid , fid , mid , snps, trait, adj.var=NA , lim =0.05, maxSNP = 3,
    nt = 10, pair.begin = FALSE, pattern.begin.mat=NA, select.criteria = "p.value" ,
    type = "gaussian" , method="forward" ,
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

    trait.type <- match(type, c("binomial","gaussian","families","survival"))

    # check input select.criteria

    select.criteria.type <- match(select.criteria, c("p.value","aic"))
    select.method <- match(method, c("forward","backward","brute.search"))

    if (is.na(select.criteria.type)) { stop("Invalid select.criteria") }
    if (is.na(select.method)) { stop("Invalid method") }

    if ( method=="backward" && ns>15 ) {
      stop ( "Backward selection appropriate to maximal 15 SNPs.\n")
    }
    if ( method=="brute.search" && ns>15 ) {
      stop ( "Extensive selection appropriate to maximal 15 SNPs.\n")
    }

    if (is.na(trait.type)) { stop("Invalid trait type") }

    switch ( type,
      gaussian =  test <- "F" ,
      binomial =  test <- "Chisq" ,
      families   =  test <- "wTDT"  ,
      survival = { print("not ready!") ; stop("end") }  )
      
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
    cat ("Method:                ",method,  "\n",sep="" )
    cat ("Max SNPs:              ",maxSNP,"\n",sep=""  )
    cat ("Number best pattern    ",nt,"\n",sep="" )
    cat ("Threshold (lim):       ",lim,"\n",sep="" )

    if (  is.null(dim(adj.var)) ) {
      cat ("Number of covariates:  0\n\n")
    } else {
      cat ("Number of covariates:  ",dim(adj.var)[2],"\n\n",sep="" )
    }
    if ( method=="brute.search" ) {
      cat ( (2^ns-1) ,  " SNP combination are obtained.\n\n" , sep=""  )
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

     if ( method == "forward" ) {

       # unadjusted
       if ( !adjusted ) {
         results <- msr.binomial.forward.unadjusted ( snps, trait, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat, select.criteria  ,
            baseline.hap , min.count )
       } else {
       
         # adjusted
         results <- msr.binomial.forward.adjusted ( snps, trait,  adj.var, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat, select.criteria  ,
            baseline.hap , min.count )
            
       }
     } # forward
  }


################################################################################
#                                                                              #
#  Quantitative trait loci (QTL) >> family="gaussian"                          #
#                                                                              #
################################################################################



  if ( type=="gaussian" ) {

     if ( method == "forward" ) {

       # unadjusted
       if ( !adjusted ) {
         results <- msr.gaussian.forward.unadjusted ( snps, trait, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat, select.criteria  ,
            baseline.hap , min.count )
       } else {

         # adjusted
         results <- msr.gaussian.forward.adjusted ( snps, trait,  adj.var, lim, maxSNP,
            nt , pair.begin, pattern.begin.mat, select.criteria  ,
            baseline.hap , min.count )

       }
     } # forward
  }


#    ############################################################################
#
#    ## forward selection
#
#    if ( method=="forward" ) {
#      results <- stepwiseforward (
#        patid=patid, snps=snps, trait=trait, adj.var=adj.var ,
#        adjusted=adjusted, lim=lim , maxSNP=maxSNP,
#        nt=nt, pair.begin=pair.begin, pattern.begin.mat=pattern.begin.mat,
#        select.criteria=select.criteria,
#        type=type, test=test , method=method ,
#        infer.Haps.separately=infer.Haps.separately ,
#        baseline.hap=baseline.hap , rest=rest ,
#        min.count=min.count )
#    }
#
#    ## backward selection
#
#    if ( method=="backward" ) {
#      results <- stepwisebackward (
#        patid=patid, snps=snps, trait=trait, adj.var=adj.var ,
#        adjusted=adjusted, lim=lim , maxSNP=maxSNP,
#        nt=nt, pair.begin=pair.begin, select.criteria=select.criteria,
#        type=type, test=test , method=method ,
#        infer.Haps.separately=infer.Haps.separately ,
#        baseline.hap=baseline.hap , rest=rest ,
#        min.count=min.count )
#    }
#
#    ## extensive selection
#    if ( method=="extensive" ) {
#      results <- stepwiseextensive (
#        patid=patid, snps=snps, trait=trait, adj.var=adj.var ,
#        adjusted=adjusted, lim=lim , maxSNP=maxSNP,
#        nt=nt, pair.begin=pair.begin, select.criteria=select.criteria,
#        type=type, test=test , method=method ,
#        infer.Haps.separately=infer.Haps.separately ,
#        baseline.hap=baseline.hap , rest=rest ,
#        min.count=min.count )
#    }


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
  }

   return ( results )
} ## end of msr ################################################################

