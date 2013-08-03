msr.gaussian.haplotype.test.unadjusted <-
function ( snps, trait, adj.var=NA , lim = 0.05,
    baseline.hap="max" , min.count=10 )
{
          # Infer haplotypes

           hapest <-  itegeppXXR ( snps , des = 0 , lim = lim )

           desres <- as.matrix(hapest$desres)

          if ( all( is.na( hapest$hap) ) ) {
               return ( list ( haplotypes=NA , nSubj=NA , df=NA ,
                   global.p.value=NA ) )
            }
          
          # check, if "rest" < min.count
          
          if ( any (colnames(desres)=="R") ) {
          
             # skipping "rest" from design matrix when "colsum < min.count"

            if ( sum ( (desres[,colnames(desres)=="R"] ) ) < min.count ) {

              desres <- desres [,colnames(desres)!="R",drop=FALSE]
              
            }
          
          }

          # select baseline haplotype and drop from design matrix
          # if there is only one haplotype then don't drop.

          if ( length(hapest$freq) > 1 ) {
             baseline <- which.max ( hapest$freq )
             desres <- desres[,-c(which(hapest$hap[baseline]==colnames(desres)))  ,drop=FALSE]
           }

          # generalized linear models
          
              # no covariates 
              
              fit.glm0 <- glm ( trait ~ 1 , family="gaussian" )
              fit.glm1 <- glm ( trait ~ . , data=as.data.frame(desres) , family="gaussian")
              aov.glm  <- anova ( fit.glm0 , fit.glm1, test = "F" )

            ## find overall p-value (goodness of fit)

            df.model   <- (aov.glm)$Df
            nind       <- length(fit.glm1$residuals)
            df         <- as.integer(df.model[!is.na(df.model)] )
            f          <- as.numeric((aov.glm)$F[!is.na((aov.glm)$F)])
            pval.model <- (aov.glm)$`Pr(>F)`           
            pval       <- as.numeric( pval.model[!is.na(pval.model)])
            aic        <- AIC(fit.glm1)
                                               
   return ( list ( haplotypes=1 , nSubj=nind ,
        df=df , F=f , global.p.value=pval , AIC=aic  )  )            
       
}
