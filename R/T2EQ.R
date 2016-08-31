T2EQ <-
function( X , Y , eq_margin , alpha=0.05 , print.results=TRUE){
  
  if (print.results) {
    cat("\n")
    cat("*******************************************************************************" , "\n\n")
    cat("********************* \t The T\U00B2-test for equivalence \t***********************" , "\n")
    cat("*************** \t \t( see Wellek (2010) \t\t***************"                    , "\n\n")
    cat("*******************************************************************************" , "\n\n")
  }
  
  m_X         <- nrow(X)                                                # sample size first sample
  n_Y         <- nrow(Y)                                                # sample size second sample
  N           <- m_X + n_Y                                              # total sample size pooled sample
  d           <- ncol(X)                                                # number of variables / dimension
  Mean_X      <- colMeans(X)                                            # Mean of first sample
  Mean_Y      <- colMeans(Y)                                            # Mean of second sample 
  
  S_X         <- var(X)                                                 # empirical covariance matrix of first sample
  S_Y         <- var(Y)                                                 # empirical covariance matrix of second sample
  S           <- 1/(N-2) * (  ( m_X - 1 ) * S_X + ( n_Y - 1 ) * S_Y  )  # pooled covariance matrix
  
  #cat(S,Mean_Y,Mean_X)
  
  M_DIST      <- (Mean_X - Mean_Y) %*% (solve(S,Mean_X - Mean_Y )  )    # estimated Mahalanobis distance
  T2          <- m_X * n_Y / N * M_DIST                                 # Hotelling's T^2
  
  numerator   <- N - 1 - d                                              # numerator of scaling factor
  denominator <- (N-2) * d                                              # denominator of scaling factor
  ncp         <- m_X * n_Y / N * eq_margin                              # noncentrality parameter
  
  if (print.results) {
    cat(" Summary statistics:"                                                           , "\n\n",    
        "Sample size of first sample X:"                   ,"\t" ,m_X                    ,   "\n",    
        "Sample size of second sample Y:"                  ,"\t" ,n_Y                    ,   "\n",
        "Sample size of pooled sample (X,Y):"              ,"\t" ,N                      ,   "\n",
        "Dimension (Number of variables):"                 ,"\t" ,d                      ,   "\n",
        "Mean of first sample:"                        ,"\t\t\t" ,Mean_X                 ,   "\n",
        "Mean of second sample:"                         ,"\t\t" ,Mean_Y                 ,   "\n",
        "Empirical covariance matrix S_X:"                 ,"\t" ,round(S_X,digits=2)    ,   "\n",
        "Empirical covariance matrix S_Y:"                 ,"\t" ,round(S_Y,digits=2)    ,   "\n",
        "Pooled empirical covariance matrix S:"            ,"\t" ,round(S,digits=2)      ,   "\n",
        "Estimated Mahalanobis distance:"                  ,"\t" ,M_DIST                 ,   "\n",      
        "Equivalence margin:"                          ,"\t\t\t" ,eq_margin              ,   "\n",      
        "Hotelling's T\U00B2:"                         ,"\t\t\t" ,T2                     ,   "\n",   
        "Noncentrality parameter:"                       ,"\t\t" ,ncp                    ,   "\n",   
        "Significance level:"                          ,"\t\t\t" ,alpha                  ,   "\n"
    )
  }
  
  quantile		<- qf( p=alpha , df1=d , df2=N-1-d , ncp=ncp )            # quantile of noncentral F-distribution
  teststat		<- T2 * numerator / denominator 	                        # teststatistic
  
  if 	( teststat < quantile ) 	{
    erg_text	<-	'Equivalence comparison successful'	                  # Test result text 
    erg       <- 1                                                      # Test result  
  }
  else 	   {
    erg_text	<-	'Equivalence comparison not successful'		            # Test result text 
    erg       <- 0                                                      # Test result  
  }
  
  p_value			<- pf( q=teststat , df1=d , df2=N-1-d , ncp=ncp )          # p-value of T^2-test for equivalence  
  
  if (print.results) {
    cat(" Teststatistic:"                                        ,"\t\t\t" ,teststat ,   "\n",    
        "Quantile of noncent. F-distribution:"                        ,"\t",quantile ,   "\n",
        "Decision in favor (1) or against (0) equivalence: "               ,erg      ,   "\n\n",
        "Test result:"                                                               ,   "\n\n",
        "\t p-value of the T\U00B2-test for equivalence: p ="              ,p_value  ,   "\n\n",
        "\t\t"                                                             ,erg_text ,   "\n\n",
        "******************************************************************************"  ,   "\n\n"
    )
  }
  result.summary <- data.frame(p.value = p_value , testresult.num = erg , testresult.text = erg_text )
  result.summary	          # Return of Test results  
}
