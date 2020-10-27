#' Compute Moran's I statistic
#'
#' Compute Moran's I statistic to measure the strength of spatial associations for each gene. 
#'
#' @examples
#' coordinates <- cbind(rep(1:10, each=10), rep(1:10, 10))
#'
#' # Completely random:
#' x <- matrix(rnorm(1000), ncol=100)
#' stats <- computeMoran(x, coordinates, radius=3)
#' stats
#'
#' # Strongly locally correlated:
#' hotspots <- coordinates[c(10, 25, 60, 85),]
#' for (i in seq_len(nrow(hotspots))) {
#'     in.range <- sqrt(colSums((t(coordinates) - hotspots[i,])^2))
#'     x[1,in.range <= 1] <- x[1,in.range <= 1] + 10
#' }
#' stats2 <- computeMoran(x, coordinates, radius=3)
#' stats2
#' 
#' @export
#' @importFrom beachmat rowBlockApply
#' @importFrom IRanges IntegerList NumericList
#' @importFrom BiocNeighbors findNeighbors
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocGenerics relist
#' @importFrom S4Vectors DataFrame
computeMoran <- function(x, coordinates, prop=0.05, radius=NULL, 
    kernel=c("uniform", "triangular", "biweight", "tricube"),
    iterations=0, BPPARAM=SerialParam())
{
    if (is.null(radius)) {
        out <- apply(coordinates, 2, range) 
        radius <- sqrt(sum(((out[2,] - out[1,]) * prop)^2))
    }

    # Organizing the weights.
    neighbors <- findNeighbors(coordinates, threshold=radius)
    idx.list <- IntegerList(neighbors$index)
    dst.list <- NumericList(neighbors$distance)

    # Compute weights.
    kernel <- match.arg(kernel)
    distances <- unlist(dst.list)
    if (kernel=="uniform") {
        w <- rep(1, length(distances))
    } else if (kernel=="triangular") {
        w <- pmax(0, radius - distances)
    } else if (kernel=="biweight") {
        w <- (1-(distances/radius)^2)^2
    } else if (kernel=="tricube") {
        w <- (1-(distances/radius)^3)^3
    }
    w.list <- relist(w, dst.list)

    keep <- idx.list > seq_along(idx.list)
    indices <- unlist(idx.list[keep])
    weights <- unlist(w.list[keep])
    runlen <- sum(keep)

    output <- rowBlockApply(x, FUN=compute_moran_part, indices=indices, weights=weights, runs=runlen, BPPARAM=BPPARAM)
    raw.stat <- unlist(lapply(output, "[[", i=1))
    squares <- unlist(lapply(output, "[[", i=2))
    quads <- unlist(lapply(output, "[[", i=3))

    # Compute Moran's I. Our set-up is guaranteed to be symmetric, so we
    # technicaly need to double the value to account for the sum in the other
    # triangular half; however, this cancels out with the doubling of the
    # weights, so it works out if we don't do anything at all.
    moran.i <- raw.stat * ncol(x) / sum(weights)

    w.list2 <- w.list[idx.list!=seq_along(idx.list)]
    p <- .moran2pvalue(moran.i, weight.list=w.list2, squares=squares, quads=quads)

    DataFrame(
        statistic=moran.i, 
        p.value=p,
        FDR=p.adjust(p, method="BH")
    )
}

#' @importFrom stats pnorm
.moran2pvalue <- function(stat, weight.list, squares, quads) 
# Stolen from good ol' Wikipedia!
{
    N <- ncol(x)
    imean <- -1/(N - 1)

    raw.weights <- unlist(weight.list)
    W <- sum(raw.weights)

    S1 <- 0.5 * sum((2 * raw.weights)^2)

    S2 <- sum((2 * sum(weight.list))^2)
    
    S3 <- N^-1 * quads / (N^-1 * squares)^2 
    
    S4 <- (N^2 - 3*N + 3)*S1 - N*S2 + 3*W^2

    S5 <- (N^2 - N)*S1 - 2*N*S2 + 6*W^2

    ivar <- (N*S4 - S3*S5)/((N-1)*(N-2)*(N-3)*W^2) - imean^2

    pnorm(stat, mean=imean, sd=sqrt(ivar), lower.tail=FALSE)
}
