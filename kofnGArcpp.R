library(Rcpp)
Rcpp::sourceCpp("~/Desktop/espingarcia/R/kofnGA.cpp")

# Inputs:
# n             The number of objects to choose from.  The algorithm chooses a subset of integers
#               from 1 to n.
# k             The number of objects to choose.
# OF            The objective function.  The first argument of OF should be an index vector of k
#               integers in the range [1, n].
# popsize       The size of the population, the number of offspring produced each generation.
# keepbest      This argument is used to implement elitism.  The keepbest least fit offspring each
#               generation are replaced by the keepbest most fit members of the previous generation.
# ngen          The number of generations to run.
# tourneysize   The number of individuals to be selected for each tournament.
# mutprob       The probability of mutation for each of the k chosen indices in each individual.
#               An index chosen for mutation simply jumps to any other unused index at random.
# initpop       A popsize-by-k matrix of starting solutions.  Possibly useful if using this
#               function in an adaptive, iterative, or parallel scheme.
# allBetaInd    Truth value of all(beta[Ind]==0)
# IM1, IM2_id, IM3_id, D  Additional arguments to the objective function.
# K.idx         A vector of more than 1 indices or a single index. In either case, the indices must be non-negative
#               integers.
# optimMeasure  1 = D-opt, 2 = A-opt, 3 = Par-spec
# Ind
# npars
# Output is a list with the following elements:
# bestsol       The best solution found.
# bestobj       The objective function value for the best solution found.
# pop           The final population, row-sorted in order of increasing objective function.
# obj           The objective function values corresponding to each row of pop.
# old           A list with elements giving for each generation: the best solution, the best
#               objective value, and the average objective value.

kofnGA.R <- function(n,k,popsize,keepbest,ngen,
                   tourneysize,mutprob,initpop,
                   allBetaInd, IM1, IM2_id, IM3_id, D, K.idx, optimMeasure, Ind, npars) {
  ##=== Basic input checking =====================================================================##
  stopifnot(n %% 1 == 0, n > 0, n >= k,
            k %% 1 == 0, k > 0,
            popsize %% 1 == 0, popsize >= 2, popsize > keepbest, popsize > tourneysize,
            keepbest %% 1 == 0, keepbest >= 0,
            ngen %% 1 == 0, ngen >= 1,
            tourneysize %% 1 == 0, tourneysize >= 2,
            mutprob >= 0, mutprob <= 1,
            all(dim(initpop) == c(popsize,k)),
            optimMeasure %in% c(1,2,3))

  indices <- 1:n;

  # When no population supplied, initialize the population using randomly-chosen indices.
  if (is.null(initpop)) {
    pop <- t(replicate(popsize,sample(indices,k)));
  }
  else {
    pop <- initpop;
  }

  # indexing in RcppArmadillo starts at 0 so we subtract each index by 1.
  pop <- pop - 1;

  # Convert Ind from logical vector to numerical vector of 0's and 1's.
  Ind <- Ind*1;

  # if optimMeasure==3 (i.e. optimMeasure is Par-spec) then K.idx MUST be not null.
  if (optimMeasure==3 & is.null(K.idx))
    stop("For a parameter-specific criterion K.idx must be provided.")

  # if optimMeasure==3 (i.e. optimMeasure is Par-spec) then K.idx MUST be a single index.
  if (optimMeasure==3 & length(K.idx) != 1)
    stop("For a parameter-specific criterion K.idx cannot be a vector of indices.")

  K.type <-"null";

  if (is.null(K.idx)) {
    K_IS_NULL <- 1; # 1 is TRUE
    K.idx = -1; # We assign an integer to K.idx so we can still pass it to the function without issue.
  }
  else if (length(K.idx) != 1) { # K.idx is a vector of indices
    K_IS_NULL <- 0; # 0 is FALSE
    K.type <- "vector";
    K.idx <- K.idx - 1; # indexing in RcppArmadillo starts at 0 so we subtract each index by 1.
  }
  else { # K.idx is a single integer
    K_IS_NULL <- 0; # 0 is FALSE
    K.type <- "int";
    K.idx <- K.idx - 1; # indexing in RcppArmadillo starts at 0 so we subtract each index by 1.
  }

  # The objective function that is called upon by kofnGA.mc (the original R code) changes its computations according
  # to the truth value of all(beta[Ind]==0) and whether K.idx is a vector of more than 1 indices (K.type=="vector"),
  # a single index (K.type=="int"), or null (K.type=="null"). Thus, we define 4 different Rcpp versions of
  # a close equivalent to kofnGA.mc, each of which is designed to handle a specific combination of the truth value
  # of all(beta[Ind]==0) and the "type" of K.idx.

  if (isTRUE(allBetaInd) & K.type == "vector") {
    out <- kofnGAv1(n, k, popsize, keepbest, ngen,
           tourneysize, mutprob, pop, IM1, IM2_id,
           IM3_id, D, K.idx, K_IS_NULL, optimMeasure, Ind, npars);
  }
  else if (isTRUE(allBetaInd) & (K.type == "int" | K.type == "null")) {
    out <- kofnGAv2(n, k, popsize, keepbest, ngen,
             tourneysize, mutprob, pop, IM1, IM2_id,
             IM3_id, D, K.idx, K_IS_NULL, optimMeasure, Ind, npars);
  }
  else if (!isTRUE(allBetaInd) & K.type == "vector") {
    out <- kofnGAv3(n, k, popsize, keepbest, ngen,
             tourneysize, mutprob, pop, IM1, IM2_id,
             IM3_id, D, K.idx, K_IS_NULL, optimMeasure, Ind, npars);
  }
  else if (!isTRUE(allBetaInd) & (K.type == "int" | K.type == "null")) {
    out <- kofnGAv4(n, k, popsize, keepbest, ngen,
             tourneysize, mutprob, pop, IM1, IM2_id,
             IM3_id, D, K.idx, K_IS_NULL, optimMeasure, Ind, npars);
  }

  class(out) <- "GAsearch"
  out
}
