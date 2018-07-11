# utils.R
# 3/19/2018
# Mark Hagemann
# Various utility functions

#' Produce a tidy data.frame from a swotlist of matrices
#'
#' @param swotlist A named list of matrices
#' @importFrom purrr map
#' @export
#'
swot_tidy <- function(swotlist) {
  nt <- ncol(swotlist[[1]])
  nx <- nrow(swotlist[[1]])

  outvecs <- map(swotlist, ~as.vector(.))

  times <- rep(1:nt, each = nx)
  locs <- rep(1:nx, nt)

  out <- as.data.frame(outvecs)
  out$time <- times
  out$loc <- locs

  attr(out, "QWBM") <- attr(swotlist, "QWBM")

  out
}


#' Subset dimensions of swot-like data.
#'
#' Returns a swotlist with all of the matrices of the input
#'  subset to keep only the specified times and location.
#'
#' @param swotlist a list of swot-like matrices
#' @param keeptimes Indices of times (columns) to keep. Default is 0, which
#'   keeps all times. Negative indices are allowed.
#' @param keeplocs Indices of locations (rows) to keep. Default is 0, which
#'   keeps all locations. Negative indices are allowed.
#' @export
swot_sset <- function(swotlist, keeptimes = 0L, keeplocs = 0L) {

  nr <- nrow(swotlist$W)
  nc <- ncol(swotlist$W)

  if (length(keeptimes) == 1 && keeptimes == 0)
    keeptimes <- 1:nc
  if (length(keeplocs == 1) && keeplocs == 0L)
    keeplocs <- 1:nr

  out <- lapply(swotlist, `[`, keeplocs, keeptimes, drop = FALSE)

  attr(out, "QWBM") <- attr(swotlist, "QWBM")

  out
}

#' Produce a swotlist from a tidy data.frame
#'
#' The inverse of swot_tidy. Produces a swotlist of matrices from a
#'   tidy data.frame of the sort produced by \code{swot_tidy}
#'
#' @param swotdf a list of swot-like matrices
#'
#' @importFrom dplyr arrange
#' @export
#'
swot_untidy <- function(swotdf) {
  matnames <- setdiff(names(swotdf), c("time", "loc"))
  # browser()
  times <- unique(swotdf$time)
  newtimes <- order(times)

  locs <- unique(swotdf$loc)
  newlocs <- order(locs)

  nc <- max(newtimes)
  nr <- max(newlocs)


  swotdf <- arrange(swotdf, time, loc)
  out <- map(swotdf[matnames],
             ~matrix(., nrow = nr, ncol = nc, byrow = FALSE))

  attr(out, "QWBM") <- attr(swotdf, "QWBM")

  out
}

#' Convert a vector to a suitably-dimensioned matrix
#'
#' Produces a matrix of the same dimensions as \code{pattern}; excactly
#'  one of these dimensions must be equal to the length of \code{vec}.
#'  If \code{pattern} is square, an error is raised.
#'
#' @param vec a vector
#' @param pattern a matrix with the desired dimensions
#'
#' @export
#'
swot_vec2mat <- function(vec, pattern) {
  nr <- nrow(pattern)
  nc <- ncol(pattern)

  if (nr == nc)
    stop("Doesn't work for square matrices")

  if (length(vec) == nr) {
    repvec <- rep(vec, nc)
  } else if (length(vec) == nc) {
    repvec <- rep(vec, each = nr)
  } else {
    stop(paste("vec length must be equal to either number",
         "of rows or columns of pattern."))
  }

  out <- matrix(repvec, nrow = nr, ncol = nc, byrow = FALSE)
  out
}

#' Apply a timelag to the rows of matrices in a swotlist.
#'
#' Produces a swotlist with each matrices' rows shifted by an amount
#'   specified by \code{lag}.
#'
#' @param swotlist a list of swot-like matrices
#' @param lags a vector with the same length as the number of locations (rows)
#'  in the swotlist data, specifying the number of indices to shift each row
#'
#' @importFrom dplyr left_join mutate filter select "%>%"
#' @export
swot_timelag <- function(swotlist, lags) {

  mats <- vapply(swotlist, is.matrix, logical(1))

  swotdf <- swot_tidy(swotlist[mats])

  lagdf <- data.frame(loc = 1:max(swotdf$loc), lag = lags)

  swotdf <- left_join(swotdf, lagdf, by = "loc") %>%
    mutate(time = time + lag) %>%
    filter(time >= min(time) + max(abs(lag)),
           time <= max(time) - max(abs(lag))) %>%
    select(-lag)
  out <- swot_untidy(swotdf)

  attr(out, "QWBM") <- attr(swotlist, "QWBM")

  out
}


# Cross-sectional area functions ------------------------------------------

#' Actual A0 values across locations of a swotlist
#'
#' Returns a vector of unobservable (via SWOT) cross-sectional areas for a given swotlist
#'
#' @param swotlist a list of swot-like matrices
#' @param zero Where to reference the zero value of dA?
#' @export
realA0 <- function(swotlist,
                   zero = c("first", "minimum", "median")) {
  zero = match.arg(zero)

  if (is.null(swotlist$A)) {
    stop("swotlist must contain an A component.\n")
  }

  dA <- swotlist$A - swot_vec2mat(swotlist$A[, 1], swotlist$A)
  if (zero != "none") {
    dA <- rezero_dA(dA, zero = zero)
  }

  out <- (swotlist$A - dA)[, 1]
  out
}

#' Adjust the zero-reference of partial cross-sectional area data.
#'
#' @param dAmat a DAWG-formatted matrix of partial cross-sectional area observations.
#' @param zero Where to reference the zero value of dA?
#' @export
rezero_dA <- function(dAmat, zero = c("first", "minimum", "median")) {
  zero <- match.arg(zero)
  if (zero == "first") {
    shifts <- dAmat[, 1]
  } else if (zero == "minimum") {
    shifts <- apply(dAmat, 1, min, na.rm = TRUE)
  } else if (zero == "median") {
    shifts <- apply(dAmat, 1, median, na.rm = TRUE)
  }
  shiftmat <- matrix(rep(shifts, ncol(dAmat)), ncol = ncol(dAmat))
  out <- dAmat - shiftmat
  out
}

#' Calculate area from A0 (a vector) and dA (a matrix).
#'
#' @param A0vec a vector of A0
#' @param dAmat a matrix of dA observations
#' @param zero If not "same" (the default), passed to rezero_dA before
#'   calculating.
swot_A <- function(A0vec, dAmat, zero = "same") {

  if (!zero == "same") {
    dAmat <- rezero_dA(dAmat = dAmat, zero = zero)
  }

  A0mat <- swot_vec2mat(A0vec, dAmat)
  Amat <- A0mat + dAmat
  Amat
}


#' Partial cross-section area from width and height
#'
#' Calculate partial cross-section area from DAWG-formatted width and height matrices
#'
#' @param w Matrix of widths
#' @param h Matrix of heights
#' @param zero Where to reference the zero value of dA?
#' @export
calcdA_mat <- function (w, h, zero = c("first", "minimum", "median")) {
  zero <- match.arg(zero)
  stopifnot(all(dim(w) == dim(h)))
  dA <- w
  for (i in 1:nrow(dA)) {
    dA[i, ] <- calcdA_vec(w[i, ], h[i, ], zero = zero)
  }
  dA
}


#' Calculate partial cross-section area from width and height vectors
#'
#' @param w vector of widths
#' @param h vector of heights
#' @param zero Where to reference the zero value of dA?
#' @export
calcdA_vec <- function(w, h, zero = c("first", "minimum", "median")) {
  zero <- match.arg(zero)
  words <- order(w)
  warr <- w[words]
  harr <- h[words]
  delh <- c(0, diff(harr))
  delA <- cumsum(warr * delh)
  dA <- 1:length(w)
  dA[words] <- delA

  if (zero == "first") {
    dA <- dA - dA[1]
  } else if (zero == "minimum") {
    dA <- dA - min(dA, na.rm = TRUE)
  } else if (zero == "median") {
    dA <- dA - median(dA, na.rm = TRUE)
  }

  dA
}

#' Purge all times or locations containing NAs
#'
#' Produces a swotlist with the all times (default) or locations that contain
#'  any NAs. This operation is performed for all matrices within the swotlist.
#'
#' @param swotlist a list of swot-like matrices
#' @param purge Which of either "times" or "locs" to purge when removing NAs.
#' @export
swot_purge_nas <- function(swotlist, purge = c("times", "locs")) {
  purge = match.arg(purge)
  getnainds <- function(mat) {
    which(is.na(mat), arr.ind = TRUE)
  }
  inddf <- purrr::map(swotlist, getnainds) %>%
    map(as.data.frame) %>%
    dplyr::bind_rows() %>%
    unique()

  if (nrow(inddf) == 0L) {
    return(swotlist)
  }

  if (purge == "times") {
    out <- swot_sset(swotlist, keeptimes = -unique(inddf[[2]]))
  } else {
    out <- swot_sset(swotlist, keeplocs = -unique(inddf[[1]]))
  }

  attr(out, "QWBM") <- attr(swotlist, "QWBM")

  out
}

#' Plot all variables, or a subset thereof, as timeseries
#'
#' Produces a ggplot of time series faceted by variable. Each facet is a DAWG
#'   plot (a la \code{plot_DAWG()})
#'
#' @param swotlist a list of swot-like matrices
#' @param vars Which variables to plot from within the swotlist.
#'   Defaults to "all".
#'
#' @importFrom ggplot2 ggplot geom_line aes geom_point facet_wrap
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @export
swot_plot <- function(swotlist, vars = "all"){
  if (!(length(vars) == 1 && vars == "all")) {
    swotlist <- swotlist[vars]
  }
  plotdf <- swot_tidy(swotlist) %>%
    gather(key = variable, value = value, -time, -loc) %>%
    mutate(loc = as.factor(loc))

  out <- ggplot(plotdf, aes(x = time, y = value, color = loc)) +
    geom_line() +
    geom_point() +
    facet_wrap(~variable, scales = "free_y")

  out
}

#' Plot all time series of a DAWG-formatted matrix.
#'
#' Produces a ggplot from a space-down, time-across matrix.
#'
#' @param dawgmat a DAWG-formatted (space-down, time-across) matrix.
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line scale_color_gradient
#' @export
plot_DAWG <- function (dawgmat) {
  dawgdf <- as.data.frame(t(dawgmat)) %>%
    setNames(1:nrow(dawgmat)) %>%
    dplyr::mutate(time = 1:ncol(dawgmat)) %>%
    melt(id.vars = "time",
         variable.name = "loc") %>%
    dplyr::mutate(loc = as.numeric(loc))
  ggplot(dawgdf, aes(x = time, y = value, group = loc)) +
    geom_line(aes(color = loc, group = loc)) +
    scale_color_gradient()
}

#' Read a netcdf file to a list.
#'
#' Produces a list of all data contained in the netcdf, with their original names.
#'
#' @param file string providing the location of a netcdf file.
#' @export
nc_list <- function(file) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("The ncdf4 package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc))

  vars <- names(nc$var)

  out <- setNames(lapply(vars, ncdf4::ncvar_get, nc = nc), make.names(vars))
  out
}

#' Read reach-averaged observations and other relevant reach-level data from
#'  a netcdf file
#'
#' Produces a list of swot-like matrices from reach-averaged slots of a netcdf,
#'   renamed for easier reference.
#'
#' @param file string providing the location of a SWOT-like netcdf file.
#' @param good_only if TRUE, subset the data to only include data from reaches
#'  specified as "good" a priori.
#' @export
nc_reach <- function (file, good_only = FALSE) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("The ncdf4 package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  nclist <- nc_list(file)
  t <- as.Date(nclist$Reach_Timeseries.t - 1, origin = "0000-01-01")
  good_reaches <- nclist$River_Info.gdrch
  W <- nclist$Reach_Timeseries.W
  H <- nclist$Reach_Timeseries.H
  S <- nclist$Reach_Timeseries.S
  dA <- calcdA_mat(w = W, h = H)
  A <- nclist$Reach_Timeseries.A
  QWBM <- nclist$River_Info.QWBM[1]
  inbounds <- 1:nrow(W)
  if (good_only)
    inbounds <- good_reaches
  Q <- nclist$Reach_Timeseries.Q
  # Qobs <- apply(Q, 2, median)

  rbnd <- as.vector(nclist$River_Info.rch_bnd)
  nreach <- length(rbnd) - 1
  upbnd = rbnd[1:nreach][inbounds]
  dnbnd = rbnd[2:(nreach + 1)][inbounds]
  x <- (upbnd + dnbnd) / 2

  ptrn <- W[inbounds, ]

  out <- list(W = W[inbounds, ], S = S[inbounds, ], dA = dA[inbounds, ],
              H = H[inbounds, ],
              t = swot_vec2mat(t, ptrn), x = swot_vec2mat(x, ptrn),
              reachid = swot_vec2mat(inbounds, ptrn), Q = Q[inbounds, ],
              A = A[inbounds, ])
  attr(out, "QWBM") <- QWBM
  out
}


#' Create BAM-ready data from a swotlist.
#'
#' Produces a bamdata object to be used with the bamr package.
#'
#' @param swotlist a list of SWOT observables
#' @param Qhat Prior guess of mean discharge. If NULL, will attempt to get
#'   from \code{attr(swotlist, "QWBM")}
#' @export
swot_bamdata <- function(swotlist, Qhat = NULL, ...) {
  if (!requireNamespace("bamr", quietly = TRUE)) {
    stop("The bamr package is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (is.null(Qhat)) {
    Qhat <- attr(swotlist, "QWBM")
    if (is.null(Qhat)) {
      stop("QWBM must be supplied if no QWBM attribute is present in swotlist.\n")
    }
  }

  bd <- bamr::bam_data(w = swotlist$W, s = swotlist$S, dA = swotlist$dA,
                       Qhat = Qhat, ...)
}



