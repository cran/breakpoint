\name{profilePlot}
\alias{profilePlot}
\title{
Mean profile plot}
\description{
Plotting function to obtain mean profile plot of the data based on the estimates of the break-points through CE method. An R object created from the CE.Normal, CE.NB or CE.ZINB is required. User can alter the axes names. 
}
\usage{
profilePlot(obj, data, x.label = "Data Sequence", y.label = "Value")
}
\arguments{
  \item{obj}{
R object created from CE.Normal, CE.NB or CE.ZINB.
}
  \item{data}{
data to be analysed. A single column array or a data frame.
}
  \item{x.label}{
x axis label. Default is "Data Sequence".
}
  \item{y.label}{
y axis label. Default is "Value".
}
}

\author{
Priyadarshana, W.J.R.M. <madawa.weerasinghe@mq.edu.au>
}

\seealso{
 \code{\link{CE.Normal}}, \code{\link{CE.NB}}, \code{\link{CE.ZINB}}.
}
\examples{
data(ch1.GM03563)
\dontrun{
## CE with four parameter beta distribution ##
obj1 <- CE.Normal(ch1.GM03563, distyp = 1, parallel =TRUE)
profilePlot(obj1)

## CE with truncated normal distribution ##
obj2 <- CE.Normal(ch1.GM03563, distyp = 2, parallel =TRUE)
profilePlot(obj2)
}}
