#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom ks Hpi.diag
#' @importFrom stats dist
#' @importFrom stats bw.nrd0
#' @importFrom stats bw.nrd
#' @importFrom utils packageVersion
#' @useDynLib cdcsis, .registration = TRUE
NULL

compute_k <- function(k, z) {
  k <- exp(6 / 5) * exp(dim(z)[2] / 4) * (dim(z)[1])^(1 / 3)
  # k <- exp(1 + dim(z)[2] / 4 + log(dim(z)[1]) / 3)
  k <- floor(k)
  if (k > dim(z)[1]) {
    k <- dim(z)[1]
  }
  k
}

gauss_z_width <- function(z) {
  if (dim(z)[2] == 1) {
    width <- stats::bw.nrd0(as.vector(z))
  } else if (dim(z)[2] <= 3) {
    width <- diag(ks::Hpi.diag(z))
  } else {
    width <- apply(z, 2, stats::bw.nrd)
  }
  width
}

discrete_z_width <- function(z) {
  dist_z <- as.matrix(dist(z))
  if (all(dist_z == 0)) {
    stop("Conditional variable is a constant and no need to check conditional independence!")
  } else {
    if (is.vector(z)) {
      discrete <- length(unique(z)) <= (length(z) / 8)
    } else {
      discrete <- nrow(unique(z)) <= (nrow(z) / 8)
    }
    if (discrete) {
      dist_z <- as.vector(dist_z)
      width <- min(dist_z[dist_z != 0]) / 2
    } else {
      width <- summary(apply(dist_z, 2, function(x) {
        min(x[x != 0])
      }))
      width <- stats::median(width) 
    }
  }
  width
}

check_k_arguments <- function(k, conditional.distance, kernel.type, z) {
  if (conditional.distance || (kernel.type == "rectangle")) {
    if (!missing(k)) {
      stopifnot(0 < k && k <= nrow(z))
    }
  }
}

check_z_hergenerity <- function(z) {
  hergenerity.z <- FALSE
  if (dim(z)[2] >= 2) {
    combn_case <- utils::combn(dim(z)[2], 2)
    hergenerity_p_value <- apply(combn_case, 2, function(index) {
      suppressWarnings(hergenerity_p <- stats::wilcox.test(z[, index[1]], z[, index[2]])[["p.value"]])
      hergenerity_p
    })
    hergenerity.z <- any(stats::p.adjust(hergenerity_p_value) < 0.05)
  }
  hergenerity.z
}

any_discrete_variable <- function(z) {
  z <- as.matrix(z)
  discrete <- apply(z, 2, function(x) {
    length(unique(x)) <= (length(x) / 1.2)
  })
  discrete
}

discrete_z_heterogeneity_width <- function(z, discrete) {
  discrete_index <- which(discrete)
  continuous_index <- which(!discrete)
  discrete_width <- c()
  continuous_width <- c()
  if (length(discrete_index) != 0) {
    discrete_width <- apply(z[, discrete_index, drop = FALSE], 2, discrete_z_width)
  }
  if (length(continuous_index) != 0) {
    hergenerity.z <- check_z_hergenerity(z[, continuous_index, drop = FALSE])
    if (hergenerity.z) {
      continuous_width <- rectangle_z_hergenerity_width2(z)
    } else {
      continuous_width <- rectangle_z_width2(z) 
    }
  }
  
  width <- rep(0, length(discrete))
  width[discrete_index] <- discrete_width
  width[continuous_index] <- continuous_width
  
  width
}

rectangle_z_hergenerity_width2 <- function(z) {
  h <- apply(z, 2, sd)
  num <- dim(z)[1]
  d <- dim(z)[2]
  d_coef <- (4 / (d + 2))^(1 / (d + 4))
  num_coef <- num^(-1 / (3 * d))
  width <- 1.5 * d_coef * num_coef * h
  width
}

rectangle_z_width <- function(z, k) {
  width <- as.matrix(dist(z))
  width <- t(apply(width, 1, sort))
  width <- mean(width[, k]) / sqrt(dim(z)[2])
  width
}

rectangle_z_width2 <- function(z) {
  width <- rectangle_z_hergenerity_width2(z)
  width <- mean(width)
  width
}

.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
  
  # Retrieve Current Version
  this.version <- utils::packageVersion("cdcsis")
  
  ## Print on Screen
  packageStartupMessage("** cdcsis")
  packageStartupMessage("**  - Conditional Feature Screening & Conditional Independence Test.")
  packageStartupMessage("** Version    : ",this.version," (",this.year,")", sep="")
  packageStartupMessage("** Maintainer : Jin Zhu (zhuj37@mail2.sysu.edu.cn)")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
}

.onUnload <- function(libpath) {
  library.dynam.unload("cdcsis", libpath)
}

#' Convert to distance matrix
#'
#' @param x 
#' @param distance 
#' @param index 
#' @return distance matrix
#' @noRd
compute_distance_matrix <- function(x, distance, index) {
  x <- as.matrix(x)
  if (!distance) {
    x <- as.matrix(dist(x))
    if (index != 1) {
      x <- x^index
    }
  }
  x
}

distance_z_width <- function(z, k) {
  width <- t(apply(z, 1, sort))
  width <- width[, k]
  width
}

#' check width arguments
#'
#' @param width 
#' @noRd
check_width_arguments <- function(width) {
  if (any(width <= 0)) {
    stop("width invalid is invalid!")
  }
}

check_xyz_arguments <- function(x) {
  if (anyNA(x)) {
    stop("input dataset include NA!")
  }
  if (any(is.infinite(x))) {
    stop("input dataset include infinite value!")
  }
}

#' check num.threads arguments
#' @param num.threads A integer number
#' @noRd
check_threads_arguments <- function(num.threads) {
  if(is.null(num.threads) | (num.threads < 1)) {
    stop("num.threads arguments is invalid!")
  }
}


#' check sample size of two dataset
#'
#' @param x numeric matrix
#' @param y numeric matrix
#' @noRd
check_sample_size <- function(x, y) {
  num_row_x <- dim(x)[1]
  num_row_y <- dim(y)[1]
  if(num_row_x != num_row_y) {
    stop("dataset with different sample size!")
  }
}


#' check num.threads arguments
#' @param num.threads A integer number
#' @noRd
check_index_arguments <- function(index) {
  if(is.null(index) | (index <= 0) | (index > 2)) {
    stop("index arguments is invalid!")
  }
}

#' check num.threads arguments
#' @param num.threads A integer number
#' @noRd
check_threshold_arguments <- function(threshold, num_col) {
  if(is.null(threshold) | (threshold <= 0) | (threshold > num_col)) {
    stop("threshold arguments is invalid!")
  }
}