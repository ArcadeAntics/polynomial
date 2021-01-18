\name{Hermite}
\alias{Hermite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hermite Polynomials
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Hermite(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
par(mar = c(4.5, 4.5, 1.7, 0.4))
plot(
    xlim = c(-3, 6), ylim = c(-10, 20),
    panel.first = grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = ~H ~ ""[n](x), main = "Hermite Polynomials"
)
n <- 0:5
col <- c("#D55E00", "#009E73", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00")
for (i in seq_along(n)) {
    x <- as.polynomial(Hermite(n[i]))
    lines(x, col = col[i], n = 1001, lwd = 2, xlim = par("usr")[1:2])
}
legend("topright",
    legend = paste0("n = ", n),
    fill = col,
    bty = "n")\dontshow{remove(n, col, i, x)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
