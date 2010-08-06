msr.meta.error.check <-
function ( data1 , control.options1 , data2 ,
        control.options2 , global.options )
{

  if ( !any(control.options1$type %in%
      c("binomial","gaussian","families","survival")  ) ) {
    stop ( "type of data 1 must be binomial, gaussian, families, or survival." )
  }
  if ( !any(control.options2$type %in%
    c("binomial","gaussian","families","survival")  ) ) {
    stop ( "type of data 2 must be binomial, gaussian, families, or survival." )
  }

  nloc1 <- dim(data1$genotypes)[2]
  nind1 <- dim(data1$genotypes)[1]
  nloc2 <- dim(data2$genotypes)[2]
  nind2 <- dim(data2$genotypes)[1]

  if ( nloc1 != nloc2 ) { stop ("both data sets must contain same number of marker.") }

  maxSNP <- global.options$maxSNP
  lim    <- global.options$lim
  nt     <- global.options$nt

  ### error
  if ( maxSNP>nloc1 ) stop ( "maxSNP>number of avaiable SNPs." )
  if ( length(data1$trait) != nind1 ) stop ( "unexpected length of trait in data 1." )
  if ( length(data2$trait) != nind2 ) stop ( "unexpected length of trait in data 2." )
  if ( (lim<0) || (lim>1) ) stop("Error: Limit of lim: 0 <= lim <= 1.")
  if ( is.null(nt)  | is.na(nt)  )  stop ( "undefined nt" )
  if ( is.null(lim) | is.na(lim) )  stop ( "undefined threshold lim" )

  if ( (global.options$meta.meth!="fisher") & (global.options$meta.meth!="rank") )  {
  
    stop ( "undefined meta method." )
  }

  if ( control.options1$type=="families" ) {
    if ( length(data1$famid) != nind1 ) stop ( "unexpected length of famid in data1." )
    if ( length(data1$patid) != nind1 ) stop ( "unexpected length of patid in data1." )
    if ( length(data1$fid)   != nind1 ) stop ( "unexpected length of fid in data1." )
    if ( length(data1$mid)   != nind1 ) stop ( "unexpected length of mid in data1." )
  }
  
  if ( control.options2$type=="families" ) {
    if ( length(data2$famid) != nind2 ) stop ( "unexpected length of famid in data2." )
    if ( length(data2$patid) != nind2 ) stop ( "unexpected length of patid in data2." )
    if ( length(data2$fid)   != nind2 ) stop ( "unexpected length of fid in data2." )
    if ( length(data2$mid)   != nind2 ) stop ( "unexpected length of mid in data2." )
  }

} # end of msr.meta.error.check

