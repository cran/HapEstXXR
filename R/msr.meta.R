msr.meta <-
function ( data1 , control.options1 , data2 , control.options2 ,
    global.options = list( maxSNP=3 , lim=.05 ,
    meta.method=c("fisher","rank") ) )
{

  cat ("R/stepwise.meta starts.\n",sep="")

  if ( !any(global.options$meta.method %in% c("fisher","rank")) ) {
    stop ("Error in msr.meta: meta.method should be fisher or rank.")
  }

  nloc  <- dim ( data1$genotypes )[2]
  nind1 <- dim ( data1$genotypes )[1]
  nind2 <- dim ( data2$genotypes )[1]

  data1$genotypes <- as.matrix  ( data1$genotypes )
  data1$trait     <- as.numeric ( data1$trait )
  data2$genotypes <- as.matrix  ( data2$genotypes )
  data2$trait     <- as.numeric ( data2$trait )
  
  maxSNP <- global.options$maxSNP
  lim    <- global.options$lim
  nt     <- global.options$nt
  
  msr.meta.error.check ( data1 , control.options1 , data2 ,
        control.options2 , global.options )

  ### print options
  cat("Opotions:\n")
  cat("Dataset 1: Number of individuals = " , nind1  , " " ,sep="")
  if ( control.options1$type == "families" ) {
    cat ( "(", length(unique(data1$famid)) , " families)\n" , sep="")
  } else {
    cat ( "\n" )
  }
  cat("Dataset 2: Number of individuals = " , nind2 ,sep="")
  if ( control.options2$type == "families" ) {
    cat ( "(", length(unique(data2$famid)) , " families)\n" , sep="")
  } else {
    cat ( "\n" )
  }
  cat("Number of SNPs = " , nloc  , "\n" ,sep="")
  cat("maxSNP = " , maxSNP  , "\n" ,sep="")
  cat("Threshold skipping haplotypes lim = " , lim  , "\n" ,sep="")
  cat("depth nt = " , nt  , "\n" ,sep="")
  cat("\n" ,sep="")
  
  test <- function (x,a) { all(x %in% a); }


  ##############################################################################
  #
  # Different data sets
  #
  #
  
  ##############################################################################
  #
  # (1) binomial (case control data) <---->  families
  #
  ##############################################################################


   if (  ((control.options1$type=="binomial") & (control.options2$type=="families")) |
         ((control.options1$type=="families") & (control.options2$type=="binomial")) ) {
         
      # first data set should be BINOMIAL and second data set should be FAMILIES !!!
      if ( control.options1$type!="binomial" ) {
        cat ("Starting from now: data set 1 and data set 2 schould be interchanged.\n")
        data.tmp <- data1
        control.options.tmp <- control.options1
        data1 <- data2
        control.options1 <- control.options2
        data1 <- data.tmp
        control.options1 <- control.options.tmp
        rm(data.tmp,control.options.tmp)
      }

     adjusted.binomial <- FALSE
     adjusted.families <- FALSE
     if ( exists("data1$adj.var") ) {
       if ( !all(is.na(data1$adj.var)) ) {
         adjusted.binomial <- TRUE
       }
     }

     if ( (adjusted.binomial==F) & (adjusted.families==F) ) {
      result <- msr.meta.binomial.unadjusted.families.unadjusted ( data1 ,
          control.options1 , data2 , control.options2 , global.options )
    }
    if ( (adjusted.binomial==T) & (adjusted.families==F) ) {
      #  result <- msr.meta.binomial.adjusted.families.unadjusted ( data1 ,
      #     control.options1 , data2 , control.options2 , global.options )
    }

  } # end of (1) ###############################################################
  


  ##############################################################################
  #
  # (2) ...
  #
  ##############################################################################

  return(result=result)

} # end of msr.meta

