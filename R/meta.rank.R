meta.rank <-
function ( p1 , p2 ) {

  if ( length(p1) != length(p2) ) {
    stop ("Error in meta.rank: object lengths differ.")
  }
  return ( rank(p1) + rank(p2)  )
 
} # end of meta.rank

