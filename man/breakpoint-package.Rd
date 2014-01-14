\name{breakpoint-package}
\alias{breakpoint-package}
\alias{breakpoint}
\docType{package}
\title{Multiple Break-Point Detection in Continuous Measurements}
\description{
The breakpoint package implements a variant of the cross-entropy (CE) method proposed in Priyadarshana and Sofronov (2014) to estimate both the number as well as the corresponding locations of break-points in biological sequences of continuous measurements. The methodology primarily built to detect multiple break-points in genomic sequences of continuous measurements. However, it can be easily extended and applied to other problems.   
}
\details{
  \tabular{ll}{
    Package: \tab breakpoint\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2013-09-16\cr
    License: \tab GPL 2.0\cr
  }
  
"breakpoint"" package provides estimates on both the number as well as the corresponding locations of break-points in continuous scale measurements. The algorithm utilizes the cross-entropy (CE) method, which is a model based stochastic optimization procedure to obtain the estimates. Current implementation of the methodology works as an exact search method in estimating the number of break-points. It selects the best solution from the solution space based on the modified BIC introduced in Zhang & Seigmund (2007). A parallel implementation of the algorithm can be carried-out in Unix/Linux/MAC OS X and Windows operating systems with the use of "doMC", "parallel", "snow" and "doSNOW" packages. }
\author{
  Priyadarshana, W.J.R.M. and Sofronov, G.
  
  Maintainer: Priyadarshana, W.J.R.M. <madawa.weerasinghe@mq.edu.au>
}
\references{
  Priyadarshana, W.J.R.M. and Sofronov, G. (2014) Multiple Break-Points Detection in array CGH Data via the Cross-Entropy Method. (Submitted)
  
  Priyadarshana,W.J.R.M., and Sofronov, G. (2012) A Modified Cross Entropy Method for Detecting Multiple Change Points in DNA Count Data. In Proc.IEEE World Congress on Computational Intelligence (CEC2012), 1020-1027.
  
  Rubinstein, R., and Kroese, D. (2004) The Cross-Entropy Method: A Unified Approach to Combinatorial Optimization, Monte-Carlo Simulation and Machine Learning. Springer-Verlag, New York.
  
  Zhang,N.R., and Seigmund, D.O. (2007) A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data. Biometrics, 63, 22-32.
}
