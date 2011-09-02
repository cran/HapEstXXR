msr.gaussian.haplotype.test.adjusted <-
function ( snps, trait, adj.var , lim = 0.05,
    baseline.hap="max" , min.count=10 )
{

          adj.var <- as.matrix(adj.var)

          # Infer haplotypes

           hapest <-  itegeppXXR ( snps , des = 0 , lim = lim )

           desres <- as.matrix(hapest$desres)

          if ( all( is.na( hapest$hap) ) ) {
               return ( list ( haplotypes=NA , nSubj=NA , df=NA , F=NA ,
                   global.p.value=NA , AIC=NA  ) )
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
              
              fit.glm0 <- glm ( trait ~ adj.var , family="gaussian" )
              fit.glm1 <- glm ( trait ~ . , data=as.data.frame(cbind(desres,adj.var)) , family="gaussian")
              aov.glm  <- anova ( fit.glm0 , fit.glm1, test = "F" )

            ## find overall p-value (goodness of fit)

            df.model   <-  (aov.glm)$Df
            nind       <-  length(fit.glm1$residuals)
            df         <- as.integer(df.model[!is.na(df.model)] )
            pval.model <- (aov.glm)$`Pr(>F)`
            pval       <- ifelse ( any(!is.na(pval.model)) , as.numeric( pval.model[!is.na(pval.model)]) , NA)
            f          <- ifelse ( any(!is.na(pval.model)) , as.numeric((aov.glm)$F[!is.na((aov.glm)$F)])  , NA)              
            aic        <- ifelse ( any(!is.na(pval.model)) , as.numeric(AIC(fit.glm1))  , NA) 
                                               
   return ( list ( haplotypes=1 , nSubj=nind ,
        df=df , F=f , global.p.value=pval , AIC=aic  )  )
        
}
