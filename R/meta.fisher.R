meta.fisher <-
function ( p1 , p2 ) {

 if ( length(p1) != length(p2) ) stop ("Error in meta.fisher: object lengths differ.")
 return (  as.numeric(apply ( data.frame(p1,p2) , 1 ,
     function(p) {1-pchisq(-2 * sum (log(p)),2*length(p) )} )) )
 
} # end of meta.fisher

