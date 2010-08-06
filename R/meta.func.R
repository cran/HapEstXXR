meta.func <-
function ( pval1 , pval2 , meta.meth ) {

   if ( meta.meth=="fisher" ) {
    return (  meta.fisher ( pval1 , pval2 )  ) 
   }
   if ( meta.meth=="rank" ) {
    return ( meta.rank ( pval1 , pval2 ) )
   }
   
} ### end of meta.func   

