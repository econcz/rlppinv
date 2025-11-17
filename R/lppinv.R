#' Solve a linear program via Convex Least Squares Programming (CLSP).
#'
#' @param c numeric vector of length \eqn{p}, optional.
#'   Objective-function coefficients. Included for API parity with
#'   Python's \code{pylppinv}; not used by CLSP.
#'
#' @param A_ub numeric matrix of size \eqn{i \times p}, optional.
#'   Matrix of inequality constraints
#'   \eqn{\mathbf{A}_{ub} \mathbf{x} \le \mathbf{b}_{ub}}.
#'
#' @param b_ub numeric vector of length \eqn{i}, optional.
#'   Right-hand side for the inequality constraints.
#'
#' @param A_eq numeric matrix of size \eqn{j \times p}, optional.
#'   Matrix of equality constraints
#'   \eqn{\mathbf{A}_{eq} \mathbf{x} = \mathbf{b}_{eq}}.
#'
#' @param b_eq numeric vector of length \eqn{j}, optional.
#'   Right-hand side for the equality constraints.
#'
#' @param non_negative logical scalar, default = \code{TRUE}.
#'   If \code{FALSE}, no default nonnegativity bound is applied.
#'
#' @param bounds NULL, \code{numeric(2)}, or list of \code{numeric(2)}.
#'   Bounds on variables. If a single pair \code{c(low, high)} is
#'   given, it is applied to all variables. If \code{NULL},
#'   defaults to \code{c(0, NA)} for each variable (non-negativity).
#'
#' @param replace_value numeric scalar or \code{NA}, default = \code{NA}.
#'   Final replacement value for any variable that violates the bounds by
#'   more than the given tolerance.
#'
#' @param tolerance numeric scalar, default = \code{sqrt(.Machine$double.eps)}.
#'   Convergence tolerance for bounds.
#'
#' @param final logical scalar, default = \code{TRUE}  
#'   If \code{FALSE}, only the first step of the CLSP estimator is performed.
#'
#' @param alpha numeric scalar, numeric vector, or \code{NULL},
#'   Regularization parameter for the second step of the CLSP estimator.
#'
#' @param ... Additional arguments passed to the \pkg{rclsp} solver.
#'
#' @return
#' An object of class \code{"clsp"} containing the fitted CLSP model.
#'
#' @seealso \link[rclsp]{clsp}
#' @seealso \link[CVXR]{CVXR-package}
#'
#' @examples
#' \donttest{
#'   ## Linear Programming via Regularized Least Squares (LPPinv)
#'   ## Underdetermined and potentially infeasible LP system
#'
#'   RNGkind("L'Ecuyer-CMRG")
#'   set.seed(123456789)
#'
#'   # sample (dataset)
#'   A_ub <- matrix(rnorm(50 * 500), nrow = 50L, ncol = 500L)
#'   A_eq <- matrix(rnorm(25 * 500), nrow = 25L, ncol = 500L)
#'   b_ub <- matrix(rnorm(50),       ncol = 1L)
#'   b_eq <- matrix(rnorm(25),       ncol = 1L)
#'
#'   # model (no default non-negativity, unique MNBLUE solution)
#'   model <- lppinv(
#'       A_ub = A_ub,
#'       A_eq = A_eq,
#'       b_ub = b_ub,
#'       b_eq = b_eq,
#'       non_negative = FALSE,
#'       final = TRUE,
#'       alpha = 1.0                                   # unique MNBLUE estimator
#'   )
#'
#'   # coefficients
#'   print("x hat (x_M hat):")
#'   print(round(model$x, 4))
#'
#'   # numerical stability (if available)
#'   if (!is.null(model$kappaC)) {
#'       cat("\nNumerical stability:\n")
#'       cat("  kappaC :", round(model$kappaC, 4), "\n")
#'   }
#'   if (!is.null(model$kappaB)) {
#'       cat("  kappaB :", round(model$kappaB, 4), "\n")
#'   }
#'   if (!is.null(model$kappaA)) {
#'       cat("  kappaA :", round(model$kappaA, 4), "\n")
#'   }
#'
#'   # goodness-of-fit diagnostics (if available)
#'   if (!is.null(model$nrmse)) {
#'       cat("\nGoodness-of-fit:\n")
#'       cat("  NRMSE :", round(model$nrmse, 6), "\n")
#'   }
#'   if (!is.null(model$x_lower)) {
#'       cat("  Diagnostic band (min):", round(min(model$x_lower), 4), "\n")
#'   }
#'   if (!is.null(model$x_upper)) {
#'       cat("  Diagnostic band (max):", round(max(model$x_upper), 4), "\n")
#'   }
#'
#'   # bootstrap NRMSE t-test (if supported by rclsp)
#'   if ("ttest" %in% names(model)) {
#'       cat("\nBootstrap t-test:\n")
#'       tt <- model$ttest(sample_size = 30L,
#'                         seed = 123456789,
#'                         distribution = "normal")
#'       for (nm in names(tt)) {
#'           cat("   ", nm, ": ", round(tt[[nm]], 6), "\n", sep = "")
#'       }
#'   }
#' }
#'
#' @export
lppinv <- function(c=NULL, A_ub=NULL, b_ub= NULL, A_eq=NULL, b_eq=NULL,
                   non_negative=TRUE, bounds=NULL, replace_value=NA_real_,
                   tolerance=sqrt(.Machine$double.eps), final=TRUE,
                   alpha=NULL, ...) {
  dots   <- list(...)
  dots$C <- NULL
  dots$S <- NULL
  dots$b <- NULL
  # assert conformability of constraint sets (A_ub, b_ub) and (A_eq, b_eq)
  if (!((( !is.null(A_ub) && !is.null(b_ub))  ||
         (!is.null(A_eq) && !is.null(b_eq))) &&
        !xor(is.null(A_ub),    is.null(b_ub))  &&
        !xor(is.null(A_eq),    is.null(b_eq))))
    stop(paste0("At least one complete constraint set (A_ub, b_ub) or ",
                "(A_eq, b_eq) must be provided."))
  if (!is.null(A_ub)) {
    A_ub <- as.matrix(A_ub)
    if (length(dim(A_ub)) == 1L) A_ub <-  matrix(A_ub, nrow=1L)
    b_ub <- matrix(as.numeric(b_ub),                   ncol=1L)
    if (nrow(A_ub) != nrow(b_ub))
      stop(sprintf("A_ub and b_ub must have the same number of rows: %d vs %d",
                   nrow(A_ub), nrow(b_ub)))
    n_vars <- ncol(A_ub)                               # number of variables
  }
  if (!is.null(A_eq)) {
    A_eq <- as.matrix(A_eq)
    if (length(dim(A_eq)) == 1L) A_eq <-  matrix(A_eq, nrow=1L)
    b_eq <- matrix(as.numeric(b_eq),                   ncol=1L)
    if (nrow(A_eq) != nrow(b_eq))
      stop(sprintf("A_eq and b_eq must have the same number of rows: %d vs %d",
                   nrow(A_eq), nrow(b_eq)))
    n_vars <- ncol(A_eq)                               # number of variables
  }
  
  # (b) Construct the right-hand side vector
  if        (is.null(bounds))      {
    bounds <- list(c(if (isTRUE(non_negative)) 0 else NA_real_, NA_real_))
  } else if (is.numeric(bounds) && length(bounds) == 2L)
    bounds <- list(bounds)
  if                              (length(bounds) == 1L) {
    bounds <- rep(bounds, n_vars)                      # replicate (low, high)
  } else if (length(bounds) !=   n_vars)
    stop(sprintf("Bounds length %d does not match number of variables %d.",
                 length(bounds), n_vars))
  bounds   <- lapply(bounds, function(v)   {
    if (length(v) != 2L)     stop("Each bounds element must have length 2.")
    vapply(v, function(x)    if (is.null(x) || is.na(x) || length(x) != 1L )
      NA_real_                 else           as.numeric(x), numeric(1L))
  })
  if (isTRUE(non_negative)   &&  any(vapply(bounds, function(v)
    (!is.na(v[1]) && v[1] < 0)  ||  (!is.na(v[2]) && v[2] < 0), logical(1L))))
    stop("Negative lower or upper bounds are not allowed in linear programs.")
  b <- matrix(numeric(0L),   ncol=1L)
  if (!is.null(b_ub)) b <-   rbind(b, b_ub)
  if (!is.null(b_eq)) b <-   rbind(b, b_eq)
  b <- rbind(b, matrix(rbind(vapply(bounds, function(v) if (is.na(v[2]))    Inf
                                    else                          v[2],
                                    numeric(1L)),
                             vapply(bounds, function(v) if (is.na(v[1]))
                               if (isTRUE(non_negative))       0 else -Inf
                               else                          v[1],
                               numeric(1L))),
                       ncol=1L))
  
  # (C), (S) Construct conformable blocks for the design matrix A
  if (!is.null(A_ub) && !is.null(A_eq)) if (ncol(A_ub) != ncol(A_eq))
    stop(sprintf(paste0("A_ub and A_eq must have the same number of columns: ",
                        "%d vs %d"), ncol(A_ub), ncol(A_eq)))
  C <- matrix(numeric(0L), ncol=n_vars)
  S <- matrix(numeric(0L), ncol=0L)
  if      (!is.null(A_ub)) {
    C  <-           A_ub
    S  <- diag(nrow(A_ub))
  }
  if      (!is.null(A_eq)) {
    C  <- rbind(C,  A_eq)
    S  <- rbind(S,  matrix(0,  nrow=nrow(A_eq),        ncol=ncol(S)))
  }
  C <-          rbind(C, diag(n_vars), diag(n_vars))
  S <-          rbind(cbind(S, matrix(0, nrow=nrow(S), ncol=      2 * n_vars)),
                      cbind(   matrix(0, nrow=n_vars,  ncol=          ncol(S)),
                               diag(n_vars),
                               matrix(0, nrow=n_vars,  ncol=          n_vars)),
                      -cbind(   matrix(0, nrow=n_vars,  ncol=ncol(S) + n_vars),
                                diag(n_vars)))
  
  # (result) Perform the estimation
  if (nrow(C) != nrow(S) || nrow(C) != nrow(b))
    stop(sprintf("Row mismatch: C=%s, S=%s, b=%s", paste(dim(C), collapse="x"),
                 paste(dim(S), collapse="x"),
                 paste(dim(b), collapse="x")))
  finite_rows   <- is.finite(b[, 1L])                  # drop rows with +-np.inf
  nonzero_cols  <- colSums(S[finite_rows, ,
                             drop=FALSE]  != 0) > 0    # reduce S width
  result        <- do.call(rclsp::clsp,   c(list(problem="general",
                                                 C=C[finite_rows, , drop=FALSE],
                                                 S=S[finite_rows,
                                                     nonzero_cols,  drop=FALSE],
                                                 b=b[finite_rows, , drop=FALSE],
                                                 tolerance=tolerance,
                                                 final=final, alpha=alpha),
                                            dots))
  
  # (result) Replace out-of-bound values with replace_value
  x        <- matrix(result$x, ncol=1L)
  x_lb     <- vapply(bounds, function(v) if (is.na(v[1]))
    if (isTRUE(non_negative))  0 else -Inf
    else                          v[1], numeric(1L))
  x_ub     <- vapply(bounds, function(v) if (is.na(v[2]))                   Inf
                     else                          v[2], numeric(1L))
  x[(x < x_lb - tolerance) | (x > x_ub + tolerance)] <- replace_value
  result$x <- matrix(x, nrow=nrow(result$x),  ncol=ncol(result$x),  byrow=TRUE)
  
  result
}
