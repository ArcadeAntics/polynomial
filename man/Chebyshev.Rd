\name{Chebyshev}
\alias{Chebyshev}
\alias{Chebyshev1}
\alias{Chebyshev2}
\alias{ChebyshevT}
\alias{ChebyshevU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chebyshev Polynomials
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chebyshev1(n)
Chebyshev2(n)
ChebyshevT(n)
ChebyshevU(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
non-negative integer; order of the Chebyshev polynomial
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
par(mar = c(4.5, 4.5, 1.7, 4.5), xpd = FALSE)
plot(
    xlim = c(-1, 1), ylim = c(-1, 1),
    panel.first = grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = ~T ~ ""[n](x), main = "Chebyshev polynomials of the first kind"
)
n <- 0:4
col <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9")
for (i in seq_along(n)) {
    x <- as.polynomial(Chebyshev1(n[i]))
    lines(x, col = col[i], n = 1001, lwd = 2)
}
par(xpd = TRUE)
legend(x = par("usr")[2L], y = par("usr")[4L],
    legend = paste0("n = ", n),
    fill = col,
    bty = "n")





par(mar = c(4.5, 4.5, 1.7, 4.5), xpd = FALSE)
plot(
    xlim = c(-1, 1), ylim = c(-4, 5),
    panel.first = grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = ~U ~ ""[n](x), main = "Chebyshev polynomials of the second kind"
)
n <- 0:4
col <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9")
for (i in seq_along(n)) {
    x <- as.polynomial(Chebyshev2(n[i]))
    lines(x, col = col[i], n = 1001, lwd = 2)
}
par(xpd = TRUE)
legend(x = par("usr")[2L], y = par("usr")[4L],
    legend = paste0("n = ", n),
    fill = col,
    bty = "n")\dontshow{remove(n, col, i, x)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
