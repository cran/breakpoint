\name{breakpoint-package}
\alias{breakpoint-package}
\alias{breakpoint}
\docType{package}
\title{Multiple Break-Point Detection via the Cross-Entropy Method}
\description{
The breakpoint package implements variants of the Cross-Entropy (CE) method proposed in Priyadarshana and Sofronov (2014, 2012a and 2012b) to estimate both the number and the corresponding locations of break-points in biological sequences of continuous and discrete measurements. The proposed method is primarily built to detect multiple break-points in genomic sequences. However, it can be easily extended and applied to other problems.
}
\details{
\tabular{ll}{
Package: \tab breakpoint\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-12-21\cr
License: \tab GPL 2.0\cr
}
"breakpoint"" package provides estimates on both the number as well as the corresponding locations of break-points. The algorithms utilize the Cross-Entropy (CE) method, which is a model based stochastic optimization procedure to obtain the estimates on location. Model selection procedures based on penalized likelihood methods are used to obtain the number of break-points. In analyzing continuous data, it uses the modified BIC introduced by Zhang & Siegmund (2007). In discrete data analysis it uses the general BIC. Current implementation of the methodology works as an exact search method in estimating the number of break-points. 

In breakpoint package v1.0, parallel computation was performed using the snow and doSNOW R packages in Windows operating systems and doMC and parallel R packages in Unix/Mac OS X operating systems as described in Priyadarshana and Sofronov (2014). In this version of the breakpoint package (v1.1), we utilize only the parallel and doParallel R packages to carry out parallel computation in both Windows and Unix/Mac OS X operating systems. However, in Windows operating systems it still uses a snow like parallel computation technique.
}
\author{
 Priyadarshana, W.J.R.M. and Sofronov, G.
  
  Maintainer: Priyadarshana, W.J.R.M. <madawa.weerasinghe@mq.edu.au>
}
\references{
 Priyadarshana, W. J. R. M., Sofronov G. (2014). Multiple Break-Points Detection in array CGH Data via the Cross-Entropy Method, IEEE/ACM Transactions on Computational Biology and Bioinformatics, no. 1, pp. 1, PrePrints, doi:10.1109/TCBB.2014.2361639, ISSN: 1545-5963.

  Priyadarshana, W. J. R. M. and Sofronov, G. (2012a). A Modified Cross- Entropy Method for Detecting Multiple Change-Points in DNA Count Data. In Proc. of the IEEE Conference on Evolutionary Computation (CEC), 1020-1027, DOI: 10.1109/CEC.2012.6256470.
  
  Priyadarshana, W. J. R. M. and Sofronov, G. (2012b). The Cross-Entropy Method and Multiple Change-Points Detection in Zero-Inflated DNA read count data. In: Y. T. Gu, S. C. Saha (Eds.) The 4th International Conference on Computational Methods (ICCM2012), 1-8, ISBN 978-1-921897-54-2.
  
  Rubinstein, R., and Kroese, D. (2004) The Cross-Entropy Method: A Unified Approach to Combinatorial Optimization, Monte-Carlo Simulation and Machine Learning. Springer-Verlag, New York.
  
  Zhang, N.R., and Siegmund, D.O. (2007) A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data. Biometrics, 63, 22-32.
}
