

system("R CMD SHLIB main.c")
x <- dyn.load(file.path(getwd(), paste0("main", .Platform$dynlib.ext)))
system2("R CMD SHLIB main.c")


setClass(
    Class = "polynomial",
    contains = c("numeric_or_complex", "oldClass"),
    slots = c(exponents = "integer"),
    prototype = prototype(.S3Class = "polynomial"))










Abel <- function (n, a = 1)
{
    n <- aslength1(as.integer(n))
    if (is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    if (n == 0L)
        1
    else if (n == 1L)
        c(0, 1)
    else {
        a <- aslength1(as.numeric_or_complex(a))
        k <- (n - 1L):0
        c(0, choose(n - 1, k) * (-a * n)^k)
    }
}


as.body <- function (x, ...)
UseMethod("as.body")


as.polynomial <- function (x, exponents = NULL, ...)
UseMethod("as.polynomial")


Bessel <- function (n)
{
    n <- aslength1(as.integer(n))
    if (is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    vapply(seq.int(0L, n), function(k) prod(seq.int(to = n +
        k, length.out = 2L * k)/c(seq_len(k), rep(2L, k))), 0)
}





.Chebyshev <- function (x)
{
    value <- function(n) NULL
    body(value) <- substitute({
        n <- aslength1(as.integer(n))
        if (is.na(n) || n < 0L)
            stop("'n' must be a non-negative integer")
        Ts <- c(list(1, c(0, x)), vector("list", n))
        f <- function(n) if (!is.null(t <- Ts[[n]]))
            t
        else (Ts[[n]] <<- c(0, 2 * f(n - 1L)) - c(f(n - 2L), 0, 0))
        f(n + 1L)
    }, list(x = x))
    value
}


Chebyshev1 <- ChebyshevT <- .Chebyshev(1)


Chebyshev2 <- ChebyshevU <- .Chebyshev(2)


remove(.Chebyshev)





coprime <- function (x, y)
gcd(x, y) == 1L


Cyclotomic <- function (n)
{
    n <- as.integer(n)
    n <- aslength1(n)
    if (is.na(n) || n <= 0L)
        stop("'n' must be a positive integer")
    x <- seq_len(n)
    x <- 2 * x[coprime(x, n)]/n
    x <- lapply(-cospi(x) - 1i * sinpi(x), c, 1)
    value <- x[[1L]]
    f <- function(x, y) {
        l <- length(x)
        value <- rep(0, -1L + length(y) + l)
        for (i in seq_along(x)) value <- value + c(rep(0, i -
            1L), x[i] * y, rep(0, l - i))
        value
    }
    for (xx in x[-1L]) value <- f(value, xx)
    Re(value)
}


degree <- function (x, ...)
UseMethod("degree")


exponents <- function (x, ...)
UseMethod("exponents")


`exponents<-` <- function (x, ..., value)
UseMethod("exponents<-")


gcd <- function (x, y)
{
    x <- as.integer(x)
    y <- as.integer(y)
    lenx <- length(x)
    leny <- length(y)
    if (!lenx || !leny)
        return(integer())
    if (lenx < leny) {
        x <- rep_len(x, leny)
        if (leny%%lenx != 0L)
            warning("longer object length is not a multiple of shorter object length")
    }
    else if (lenx > leny) {
        y <- rep_len(y, lenx)
        if (lenx%%leny != 0L)
            warning("longer object length is not a multiple of shorter object length")
    }
    c(.mapply(function(x, y) {
        if (is.na(x) || is.na(y))
            return(NA_integer_)
        else if (!x)
            return(y)
        else if (!y)
            return(x)
        else if (x < 2L || y < 2L)
            return(1L)
        fx <- seq_len(x)
        fx <- fx[!x%%fx]
        fy <- seq_len(y)
        max(fx[fx %in% fy[!y%%fy]])
    }, list(x, y), NULL), recursive = TRUE, use.names = FALSE)
}


Hermite <- function (n)
{
    a <- 1
    for (i in seq.int(0L, length.out = n))
        a <- c(0, a) - c(a[-1L] * seq_len(i), 0, 0)
    simplifyCoefficients(a = a, k = seq.int(0L, along.with = a))
}


integral <- function (x, ...)
UseMethod("integral")


is0 <- function (x)
!is.na(x) & !x


is.polynomial <- function (x)
inherits(x, "polynomial")


make.exponents <- function (exponents)
{
    exponents <- as.integer(exponents)
    exponents[is.na(exponents) | exponents < 0L] <- 0L
    exponents
}


polynomial <- function (...)
{
    if (!...length())
        return(.polynomial(numeric()))
    x <- lapply(list(...), as.polynomial)
    simplifyExponents(a = c(lapply(x, coef), recursive = TRUE,
        use.names = FALSE), k = c(lapply(x, exponents),
        recursive = TRUE, use.names = FALSE))
}


reduce <- function (x, ...)
UseMethod("reduce")


remove0 <- function (x)
{
    if (!is.polynomial(x))
        x <- as.polynomial(x)
    i <- !is0(x)
    .polynomial(coef(x)[i], exponents(x)[i])
}


simplifyCoefficients <- function (x, a = coef(x), k = exponents(x), na.rm = FALSE, finite = FALSE)
{
    i <- !is0(a)
    if (finite)
        i <- i & is.finite(a)
    else if (na.rm)
        i <- i & !is.na(a)
    .polynomial(a[i], k[i])
}


simplifyExponents <- function (x, a = coef(x), k = exponents(x))
{
    if (anyDuplicated(k)) {
        uk <- unique(k)
        .polynomial(c(lapply(split(a, factor(k, uk)), sum),
            recursive = TRUE, use.names = FALSE), uk)
    }
    else .polynomial(a, k)
}


sortExponents <- function (x, a = coef(x), k = exponents(x))
{
    i <- order(k, decreasing = decreasing)
    .polynomial(a[i], k[i])
}


`.exponentsPoly<-` <- function (x, make.exponents = FALSE, value)
{
    if (!is.polynomial(x))
        x <- as.polynomial(x)
    if (is.null(value)) {
        attr(x, "exponents") <- integer()
        return(x)
    }
    n <- length(x)
    if (is.object(value) || !is.integer(value))
        value <- as.integer(value)
    if (length(value) != n) {
        if (isFALSE(make.exponents))
            stop("invalid 'exponents' length")
        else if (is.na(make.exponents)) {
            attr(x, "exponents") <- integer()
            return(x)
        }
        else if (!isTRUE(make.exponents))
            stop("invalid 'make.exponents'")
        else if ((nv <- length(value)) < n)
            value <- c(value, rep_len(value[nv], n - nv))
        else value <- value[seq_len(n)]
    }
    if (anyNA(value)) {
        if (isFALSE(make.exponents))
            stop("missing values in 'exponents' are not allowed")
        else if (is.na(make.exponents)) {
            attr(x, "exponents") <- integer()
            return(x)
        }
        else if (!isTRUE(make.exponents))
            stop("invalid 'make.exponents'")
        else value <- make.exponents(value)
    }
    else if (any(value < 0L)) {
        if (isFALSE(make.exponents))
            stop("negative values in 'exponents' are not allowed")
        else if (is.na(make.exponents)) {
            attr(x, "exponents") <- integer()
            return(x)
        }
        else if (!isTRUE(make.exponents))
            stop("invalid 'make.exponents'")
        else value <- make.exponents(value)
    }
    attributes(value) <- NULL
    attr(x, "exponents") <- value
    x
}


.polynomial <- function (xx, exponents = integer(), cl = "polynomial")
{
    class(xx) <- cl
    attr(xx, "exponents") <- exponents
    xx
}


.simplifyPoly <- function (x, a = coef(x), k = exponents(x))
{
    if (anyDuplicated(k)) {
        uk <- unique(k)
        a <- c(lapply(split(a, factor(k, uk)), sum), recursive = TRUE,
            use.names = FALSE)
        k <- uk
    }
    i <- !is0(a)
    .polynomial(a[i], k[i])
}










as.body.default <- function (x, ...)
as.body.polynomial(as.polynomial(x), ...)


as.body.polynomial <- function (x, var = "x", ...)
{
    if (!is.call(var))
        var <- as.symbol(var)
    x <- simplifyCoefficients(x, ...)
    if (!length(x))
        return(substitute(rep(0, length(var)), list(var = var)))
    a <- coef(x)
    k <- exponents(x)
    if (length(x) == 1L) {
        num <- if (k == 1L)
            var
        else call("^", var, Re(k))
        body <- if (is.na(a))
            call("*", a, num)
        else if (a == 1)
            num
        else if (a == -1)
            call("-", num)
        else call("*", a, num)
    }
    else {
        nums <- lapply(Re(k), function(y) if (y == 1)
            var
        else if (y == 0)
            NULL
        else call("^", var, y))
        .sign <- function(x) {
            if (is.numeric(x)) {
                value <- sign(x)
                value[is.na(value)] <- 1
            }
            else {
                value <- sign(Re(x))
                i <- is.na(value) | value == 0
                value[i] <- sign(Im(x[i]))
            }
            value
        }
        signs <- .sign(a[-1L])
        parts <- .mapply(function(e1, e2) if (is.null(e2))
            as.numeric_or_complex(e1)
        else if (istrue(e1 == 1))
            e2
        else call("*", as.numeric_or_complex(e1), e2), list(a * c(1, signs),
            nums), NULL)
        i <- signs == -1
        signs[i] <- "-"
        signs[!i] <- "+"
        body <- parts[[1L]]
        for (i in seq_along(signs)) body <- call(signs[i], body,
            parts[[i + 1L]])
    }
    body
}


as.function.polynomial <- function (x, envir = parent.frame(), xname = "x", var = xname, ...)
{
    xname <- as.symbol(xname)
    value <- function(x) NULL
    environment(value) <- envir
    names(formals(value)) <- as.character(xname)
    body(value) <- as.body(x, var = var)
    value
}


as.list.polynomial <- function (x, ...)
.mapply(.polynomial, list(coef(x), exponents(x)), NULL)





as.polynomial.array <- function (x, exponents = NULL, ...)
{
    d <- dim(x)
    if (length(d) == 1L)
        as.polynomial.vector(drop(x), exponents, ...)
    else if (length(d) == 2L)
        as.polynomial.matrix(x, exponents, ...)
    else {
        rn <- dimnames(x)[[1L]]
        dim(x) <- c(d[1L], prod(d[-1L]))
        if (length(rn))
            dimnames(x)[[1L]] <- rn
        as.polynomial.matrix(x, exponents, ...)
    }
}


as.polynomial.character <- function (x, exponents = NULL, ...)
{
    n <- length(x)
    if (!(is.null(exponents) || (is.integer(exponents) && length(exponents) ==
        n))) {
        warning(gettextf("'exponents' is not an integer vector of length %d -- omitting it. Will be an error!",
            n), domain = NA)
        exponents <- NULL
    }
    n <- names(x)
    value <- as.numeric_or_complex(x)
    names(value) <- n
    class(value) <- "polynomial"
    .exponentsPoly(value) <- exponents
    value
}


as.polynomial.complex <- as.polynomial.character


as.polynomial.data.frame <- function (x, exponents = NULL, make.exponents = TRUE, ...)
{
    rn <- if (.row_names_info(x) > 0L)
        row.names(x)
    value <- rep(0, .row_names_info(x, 2L))
    for (e2 in lapply(x, as.numeric_or_complex)) value <- value + e2
    names(value) <- rn
    class(value) <- "polynomial"
    if (is.null(exponents) || length(exponents) != length(value))
        attr(value, "exponents") <- integer()
    else .exponentsPoly(value, make.exponents = make.exponents) <- exponents
    value
}


as.polynomial.default <- function (x, ...)
stop(gettextf("cannot coerce class %s to a polynomial", sQuote(deparse(class(x))[1L])),
    domain = NA)


as.polynomial.integer <- as.polynomial.character


as.polynomial.logical <- as.polynomial.character


as.polynomial.matrix <- function (x, exponents = NULL, make.exponents = TRUE, ...)
{
    d <- dim(x)
    rn <- dimnames(x)[[1L]]
    if (!is.numeric_or_complex(x)) {
        x <- as.numeric_or_complex(x)
        dim(x) <- d
    }
    value <- rowSums(x)
    names(value) <- rn
    class(value) <- "polynomial"
    if (is.null(exponents) || length(exponents) != length(value))
        attr(value, "exponents") <- integer()
    else .exponentsPoly(value, make.exponents = make.exponents) <- exponents
    value
}


as.polynomial.NULL <- function (x, ...)
.polynomial(numeric())


as.polynomial.numeric <- as.polynomial.character


as.polynomial.polynomial <- function (x, exponents = NULL, ...)
{
    cl <- oldClass(x)
    i <- match("polynomial", cl)
    if (i > 1L)
        class(x) <- cl[-(1L:(i - 1L))]
    if (!is.null(exponents)) {
        n <- length(x)
        if (length(exponents) == n)
            .exponentsPoly(x) <- exponents
        else stop(sprintf(ngettext(n, "invalid 'exponents', length %d for a polynomial with %d coefficient",
            "invalid 'exponents', length %d for a polynomial with %d coefficients"),
            length(exponents), n), domain = NA)
    }
    x
}


as.polynomial.vector <- as.polynomial.character





c.polynomial <- function (...)
{
    args <- lapply(list(...), as.polynomial)
    reduce(.polynomial(c(lapply(args, coef), recursive = TRUE),
        c(lapply(args, exponents), recursive = TRUE)))
}


coef.polynomial <- function (object, ...)
c(unclass(object))


degree.polynomial <- function (x, ...)
max(0L, exponents(x)[!is0(x)], na.rm = TRUE)


deriv.polynomial <- function (expr, ...)
{
    k <- exponents(expr)
    reduce(.polynomial(k * coef(expr), k - 1L))
}


exponents.polynomial <- function (x, ...)
{
    value <- attr(x, "exponents")
    if (length(value) || !length(x))
        value
    else seq.int(0L, along.with = x)
}


`exponents<-.polynomial` <- function (x, ..., value)
`.exponentsPoly<-`(x, value = value)


integral.polynomial <- function (x, lower, upper, C, ...)
{
    a <- coef(x)
    k <- exponents(x) + 1L
    if (missing(C))
        diff(predict(.polynomial(a/k, k), c(lower, upper)))
    else .polynomial(c(C, a/k), c(0L, k))
}


lines.polynomial <- function (x, ..., type = "l", add = TRUE)
{
    if (dev.cur() == 1L)
        stop("plot.new has not been called yet")
    plot(x, ..., type = type, add = TRUE)
}


plot.polynomial <- function (x, y = NULL, to = NULL, from = y, n = 101, add = FALSE,
    type = "l", xname = "x", xlab = xname, ylab = NULL, log = NULL,
    xlim = NULL, ...)
{
    if (dev.cur() == 1L && !isFALSE(add)) {
        warning("'add' will be ignored as there is no existing plot")
        add <- FALSE
    }
    addF <- isFALSE(add)
    if (is.null(ylab))
        ylab <- substitute(~subthis, list(subthis = as.body(signif(x, 7L), finite = TRUE)))
    if (is.null(from) || is.null(to)) {
        xl <- if (!is.null(xlim))
            xlim
        else if (!addF) {
            pu <- par("usr")[1:2]
            if (par("xaxs") == "r")
                pu <- extendrange(pu, f = -1/27)
            if (par("xlog"))
                10^pu
            else pu
        }
        else c(0, 1)
        if (is.null(from))
            from <- xl[1L]
        if (is.null(to))
            to <- xl[2L]
    }
    lg <- if (length(log))
        log
    else if (!addF && par("xlog"))
        "x"
    else ""
    y <- x
    if (grepl("x", lg, fixed = TRUE)) {
        if (from <= 0 || to <= 0)
            stop("'from' and 'to' must be > 0 with log=\"x\"")
        x <- exp(seq.int(log(from), log(to), length.out = n))
    }
    else x <- seq.int(from, to, length.out = n)
    y <- predict(y, x)
    if (isTRUE(add))
        lines(x = x, y = y, type = type, ...)
    else plot(x = x, y = y, type = type, xlab = xlab, ylab = ylab,
        xlim = xlim, log = lg, ...)
    invisible(list(x = x, y = y))
}


points.polynomial <- function (x, ..., type = "p", add = TRUE)
{
    if (dev.cur() == 1L)
        stop("plot.new has not been called yet")
    plot(x, ..., type = type, add = TRUE)
}


predict.polynomial <- function (object, newdata = seq.int(0, 1, 0.01), ...)
{
    a <- coef(object)
    k <- exponents(object)
    value <- rep(0, length(newdata))
    i <- is.finite(a)
    if (!all(i)) {
        warning(sprintf(ngettext(sum(!i), "removed %d non-finite coefficient",
            "removed %d non-finite coefficients"), sum(!i)),
            domain = NA)
        a <- a[i]
        k <- k[i]
    }
    i <- a != 0
    a <- a[i]
    k <- k[i]
    for (e2 in .mapply(function(a, k) a * newdata^k, list(a,
        k), NULL)) value <- value + e2
    value
}


print.polynomial <- function (x, ...)
{
    if (length(x)) {
        cat("Coefficients:\n")
        print(coef(x))
        cat("\nExponents:\n")
        print(`names<-`(exponents(x), NULL))
    }
    else cat("No coefficients\n")
    invisible(x)
}


reduce.polynomial <- function (x, na.rm = FALSE, finite = FALSE, sort = FALSE, decreasing = FALSE,
    ...)
{
    a <- coef(x)
    k <- exponents(x)
    i <- !is0(a) & !is.na(k) & k >= 0L
    if (finite)
        i <- i & is.finite(a)
    else if (na.rm)
        i <- i & !is.na(a)
    a <- a[i]
    k <- k[i]
    if (anyDuplicated(k)) {
        uk <- unique(k)
        a <- c(lapply(split(a, factor(k, uk)), sum),
            recursive = TRUE, use.names = FALSE)
        k <- uk
        i <- !is0(a)
        if (finite)
            i <- i & is.finite(a)
        else if (na.rm)
            i <- i & !is.na(a)
        a <- a[i]
        k <- k[i]
    }
    n <- names(x)
    a <- as.numeric_or_complex(a)
    names(a) <- n
    if (sort)
        sortExponents(a = a, k = k, decreasing = decreasing)
    else .polynomial(a, k)
}


`-.polynomial` <- function (e1, e2)
{
    if (nargs() == 1L)
        return(.polynomial(-coef(e1), exponents(e1)))
    c(e1, -e2)
}


`*.polynomial` <- function (e1, e2)
{
    if (nargs() == 1L)
        stop("invalid unary operator")
    e1 <- as.polynomial(e1)
    e2 <- as.polynomial(e2)
    .polynomial(c(lapply(coef(e1), `*`, coef(e2)), recursive = TRUE, use.names = FALSE),
        c(lapply(exponents(e1), `+`, exponents(e2)), recursive = TRUE, use.names = FALSE))
}


`/.polynomial` <- function (e1, e2)
{
    if (nargs() == 1L)
        stop("invalid unary operator")
    stop("currently undefined")
    e1 <- as.polynomial(e1)
    e2 <- as.polynomial(e2)
    reduce(.polynomial(c(lapply(coef(e1), `*`, coef(e2)), recursive = TRUE),
        c(lapply(exponents(e1), `+`, exponents(e2)), recursive = TRUE)))
}


`[.polynomial` <- function (x, i, ...)
{
    if (is.character(i)) {
        xx <- seq_along(x)
        names(xx) <- names(x)
        i <- xx[i]
    }
    .polynomial(coef(x)[i], exponents(x)[i])
}


`[[.polynomial` <- function (x, i, ...)
{
    if (is.character(i)) {
        xx <- seq_along(x)
        names(xx) <- names(x)
        i <- xx[[i]]
    }
    .polynomial(coef(x)[[i]], exponents(x)[[i]])
}


`+.polynomial` <- function (e1, e2)
{
    if (nargs() == 1L)
        return(e1)
    c(e1, e2)
}
