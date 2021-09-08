#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// === Helper Functions ========================================================

/*
 * Class required for random_rank helper function.
 * Sources:
 * https://stackoverflow.com/questions/39153082/rcpp-rank-function-that-does-average-ties
 * https://stackoverflow.com/questions/30822729/create-ranking-for-vector-of-double/30827731#30827731
 */

class Comparator {
private:
  const Rcpp::NumericVector& ref;

  bool is_na(double x) const
  {
    return Rcpp::traits::is_na<REALSXP>(x);
  }

public:
  Comparator(const Rcpp::NumericVector& ref_)
    : ref(ref_)
  {}

  bool operator()(const int ilhs, const int irhs) const
  {
    double lhs = ref[ilhs], rhs = ref[irhs];
    if (is_na(lhs)) return false;
    if (is_na(rhs)) return true;
    return lhs < rhs;
  }
};

/*
 * Equivalent to R's rank function when ties.method="random".
 * Sources:
 * https://stackoverflow.com/a/39170596/16328349
 * https://stackoverflow.com/questions/30822729/create-ranking-for-vector-of-double/30827731#30827731
 */
// [[Rcpp::export]]
Rcpp::NumericVector random_rank(Rcpp::NumericVector x)
{
  R_xlen_t sz = x.size();
  Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
  std::sort(w.begin(), w.end(), Comparator(x));

  Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      r[w[i+k]] = i + k + 1;
    }
  }

  return r;
}

/*
 * Equivalent to R's order function.
 * Sources:
 * https://stackoverflow.com/a/21613591/16328349
 */
arma::urowvec order_fun(arma::rowvec arma_vec) {
  NumericMatrix rcpp_vec = wrap(arma_vec) ;
  NumericVector sorted = clone(rcpp_vec).sort();
  IntegerVector rcpp_ord = match(sorted, rcpp_vec);
  arma::urowvec arma_ord = as<arma::urowvec>(rcpp_ord);
  return arma_ord;
}


/*
 *  Fills a matrix by assigning each row a random sample of cols # of values from the vector indices.
 *  This function is used to imitate the R code t(replicate(popsize,sample(1:popsize,tourneysize)))
 *  which is used in the construction of Tourneys1 and Tourneys2 in kofnGA_mc.R:
 *  Tourneys1[,] <- t(replicate(popsize,sample(1:popsize,tourneysize)))
 *  Tourneys2[,] <- t(replicate(popsize,sample(1:popsize,tourneysize)))
 *  Inputs:
 *  mat       The matrix to be filled.
 *  rows      The number of rows in mat.
 *  cols      The number of columns in mat.
 *  indices   A vector of integers to sample from.
 */
arma::umat rep_transpose(arma::umat mat, int rows, int cols, IntegerVector indices) {
  for (arma::uword i = 0; i < rows; i++) {
    IntegerVector ith_row = Rcpp::sample(indices, cols);
    mat.row(i) = as<arma::urowvec>(ith_row);
  }
  return mat;
}

/*
 *  Equivalent to the R function .optMsure (defined in optimJTC_v2.R) when k (the argument for the R function .optMsure)
 *  k is a single index or null.
 *  Inputs:
 *  vcov          The matrix vcov or invV1 calculated in the corresponding objective function before calling upon our Rcpp equivalent
 *                of optMsuref.
 *  k             Equivalent to K_idx (K_idx is just the Rcpp version of the R variable K.idx). K.idx is a single
 *                index or null in this case thus k is -1 (K.idx is null) or a non-negative index (K.idx is a single index).
 *  K_IS_NULL     Integer value representing whether K.idx is null (1) or not null (0).
 *  optimMeasure  1 = D-opt, 2 = A-opt, 3 = Par-spec.
 */
double optMsuref_int(arma::mat vcov, arma::uword k, int K_IS_NULL, int optimMeasure) {
  switch(optimMeasure) {
  case 1: // "D-opt"
    if (K_IS_NULL) {
      return arma::det(vcov);
    } else {
      return arma::det(vcov(k,k));
    }
    break;

  case 2: // "A-opt"
    if (K_IS_NULL) {
      return arma::trace(vcov);
    } else {
      return vcov(k,k);
    }
    break;

  case 3: // "Par-spec"
    return vcov(k,k);
    break;

  default:
    printf("Error!");
    exit(1);
  }
}

/*
 *  Equivalent to the R function .optMsure (defined in optimJTC_v2.R) when k (the argument for the R function .optMsure)
 *  is a vector of more than 1 indices.
 *  Inputs:
 *  vcov          The matrix vcov or invV1 calculated in the corresponding ObjFun before calling upon optMsuref.
 *  k             Equivalent to K_idx (K_idx is just the Rcpp version of the R variable K.idx). k is a vector
 *                of more than 1 indices in this case.
 *  K_IS_NULL     Integer value representing whether K.idx is null (1) or not null (0).
 *  optimMeasure  1 = D-opt, 2 = A-opt, 3 = Par-spec.
 */
double optMsuref_vec(arma::mat vcov, arma::urowvec k, int K_IS_NULL, int optimMeasure) {
  switch(optimMeasure) {
  case 1: // "D-opt"
    if (K_IS_NULL) {
      return arma::det(vcov);
    } else {
      return arma::det(vcov.submat(k,k));
    }
    break;

  case 2: // "A-opt"
    if (K_IS_NULL) {
      return arma::trace(vcov);
    } else {
      return arma::trace(vcov.submat(k,k));
    }
    break;
  default:
    printf("Error!");
    exit(1);
  }
}

/*
 * Equivalent to the R function Assess_fitness (defined in optimJTC_v2.R) when
 * all(beta[Ind]==0) is true and K.idx is a vector of more than 1 indices.
 * v              The vector which is to be assessed.
 * Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D Rcpp form of additional arguments to the objective function.
 * K_idx          RcppArmadillo version of the R variable K.idx. A vector of more than 1 indices in this case.
 * K_IS_NULL      Integer value representing whether K.idx is null (1) or not null (0).
 * optimMeasure   1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Arma_Ind
 * npars
 *
 */
double  ObjFun_v1_vec(arma::urowvec v, arma::mat Arma_IM1,
                      arma::mat Arma_IM2_id, arma::mat Arma_IM3_id, arma::mat Arma_D,
                      arma::urowvec K_idx, int K_IS_NULL, int optimMeasure, arma::urowvec Arma_Ind, int npars) {

  arma::mat submat1 = Arma_IM2_id.rows(v);
  Rcpp::NumericMatrix Y1 = Rcpp::wrap(submat1);
  Rcpp::NumericVector v1 = Rcpp::colSums(Y1);
  v1.attr("dim") = Dimension(npars);
  NumericMatrix m1 = as<NumericMatrix>(v1);
  arma::mat MAT1 = as<arma::mat>(m1);

  arma::mat submat2 = Arma_IM3_id.rows(v);
  Rcpp::NumericMatrix Y2 = Rcpp::wrap(submat2);
  Rcpp::NumericVector v2 = Rcpp::colSums(Y2);
  v2.attr("dim") = Dimension(npars);
  NumericMatrix m2 = as<NumericMatrix>(v2);
  arma::mat MAT2 = as<arma::mat>(m2);

  arma::mat FIM = Arma_IM1 + MAT1 - MAT2;
  arma::mat F1 = Arma_D.t() * FIM * Arma_D;

  arma::uvec Ind_ids = arma::find(Arma_Ind = 1);
  Arma_Ind.replace(0, 2);
  arma::uvec Not_Ind_ids = arma::find(Arma_Ind = 2);

  arma::mat invF1_noInd = arma::pinv(F1.submat(Not_Ind_ids, Not_Ind_ids));
  arma::mat V1 = F1.submat(Ind_ids,Ind_ids) - F1.submat(Ind_ids,Not_Ind_ids) * invF1_noInd * F1.submat(Not_Ind_ids,Ind_ids);
  arma::mat invV1 = inv(V1);
  return optMsuref_vec(invV1, K_idx, K_IS_NULL, optimMeasure);
}

/*
 * Equivalent to Assess_fitness when all(beta[Ind]==0) is true and
 * K.idx is a single index or null.
 * Inputs:
 * v              The vector which is to be assessed.
 * Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D Rcpp form of additional arguments to the objective function.
 * K_idx          K.idx is a single index or null in this case thus K_idx is -1 (K.idx is null) or a
 *                non-negative index (K.idx is a single index).
 * K_IS_NULL      Integer value representing whether K.idx is null (1) or not null (0).
 * optimMeasure   1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Arma_Ind
 * npars
 */
double  ObjFun_v1_int(arma::urowvec v, arma::mat Arma_IM1,
                      arma::mat Arma_IM2_id, arma::mat Arma_IM3_id, arma::mat Arma_D,
                      int K_idx, int K_IS_NULL, int optimMeasure, arma::urowvec Arma_Ind, int npars) {

  arma::mat submat1 = Arma_IM2_id.rows(v);
  Rcpp::NumericMatrix Y1 = Rcpp::wrap(submat1);
  Rcpp::NumericVector v1 = Rcpp::colSums(Y1);
  v1.attr("dim") = Dimension(npars, npars);
  NumericMatrix m1 = as<NumericMatrix>(v1);
  arma::mat MAT1 = as<arma::mat>(m1);

  arma::mat submat2 = Arma_IM3_id.rows(v);
  Rcpp::NumericMatrix Y2 = Rcpp::wrap(submat2);
  Rcpp::NumericVector v2 = Rcpp::colSums(Y2);
  v2.attr("dim") = Dimension(npars, npars);
  NumericMatrix m2 = as<NumericMatrix>(v2);
  arma::mat MAT2 = as<arma::mat>(m2);

  arma::mat FIM = Arma_IM1 + MAT1 - MAT2;
  arma::mat F1 = Arma_D.t() * FIM * Arma_D;

  Arma_Ind.replace(0, 2);
  arma::uvec Not_Ind_ids = arma::find(Arma_Ind == 2);
  arma::uvec Ind_ids = arma::find(Arma_Ind == 1);

  arma::mat invF1_noInd = arma::pinv(F1.submat(Not_Ind_ids, Not_Ind_ids));
  arma::mat V1 = F1.submat(Ind_ids,Ind_ids) - F1.submat(Ind_ids,Not_Ind_ids) * invF1_noInd * F1.submat(Not_Ind_ids,Ind_ids);
  arma::mat invV1 = inv(V1);

  return optMsuref_int(invV1, K_idx, K_IS_NULL, optimMeasure);
}

/*
 * Equivalent to Assess_fitness when all(beta[Ind]==0) is false and
 * K.idx is a vector of more than 1 indices.
 * Inputs:
 * v              The vector which is to be assessed.
 * Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D Rcpp form of additional arguments to the objective function.
 * K_idx          RcppArmadillo version of the R variable K.idx. A vector of more than 1 indices in this case.
 * K_IS_NULL      Integer value representing whether K.idx is null (1) or not null (0).
 * optimMeasure   1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Arma_Ind
 * npars
 *
 */
double  ObjFun_v2_vec(arma::urowvec v, arma::mat Arma_IM1,
                      arma::mat Arma_IM2_id, arma::mat Arma_IM3_id, arma::mat Arma_D,
                      arma::urowvec K_idx, int K_IS_NULL, int optimMeasure, arma::urowvec Arma_Ind, int npars) {

  arma::mat submat1 = Arma_IM2_id.rows(v);
  Rcpp::NumericMatrix Y1 = Rcpp::wrap(submat1);
  Rcpp::NumericVector v1 = Rcpp::colSums(Y1);
  v1.attr("dim") = Dimension(npars);
  NumericMatrix m1 = as<NumericMatrix>(v1);
  arma::mat MAT1 = as<arma::mat>(m1);

  arma::mat submat2 = Arma_IM3_id.rows(v);
  Rcpp::NumericMatrix Y2 = Rcpp::wrap(submat2);
  Rcpp::NumericVector v2 = Rcpp::colSums(Y2);
  v2.attr("dim") = Dimension(npars);
  NumericMatrix m2 = as<NumericMatrix>(v2);
  arma::mat MAT2 = as<arma::mat>(m2);

  arma::mat FIM = Arma_IM1 + MAT1 - MAT2;
  arma::mat F1 = Arma_D.t() * FIM * Arma_D;

  arma::mat vcov = inv(F1);

  return optMsuref_vec(vcov, K_idx, K_IS_NULL, optimMeasure);
}

/*
 * Equivalent to Assess_fitness when all(beta[Ind]==0) is false and
 * K.idx is a single index or null.
 * Inputs:
 * v              The vector which is to be assessed.
 * IM1, IM2_id, IM3_id, D Additional arguments to the objective function.
 * K_idx          K.idx is a single index or null in this case thus K_idx is -1 (K.idx is null)
 *                or a non-negative index (K.idx is a single index).
 * K_IS_NULL      Integer value representing whether K.idx is null (1) or not null (0).
 * optimMeasure   1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Arma_Ind
 * npars
 */
double  ObjFun_v2_int(arma::urowvec v, arma::mat Arma_IM1,
                      arma::mat Arma_IM2_id, arma::mat Arma_IM3_id, arma::mat Arma_D,
                      int K_idx, int K_IS_NULL, int optimMeasure, arma::urowvec Arma_Ind, int npars) {

  arma::mat submat1 = Arma_IM2_id.rows(v);
  Rcpp::NumericMatrix Y1 = Rcpp::wrap(submat1);
  Rcpp::NumericVector v1 = Rcpp::colSums(Y1);
  v1.attr("dim") = Dimension(npars);
  NumericMatrix m1 = as<NumericMatrix>(v1);
  arma::mat MAT1 = as<arma::mat>(m1);

  arma::mat submat2 = Arma_IM3_id.rows(v);
  Rcpp::NumericMatrix Y2 = Rcpp::wrap(submat2);
  Rcpp::NumericVector v2 = Rcpp::colSums(Y2);
  v2.attr("dim") = Dimension(npars);
  NumericMatrix m2 = as<NumericMatrix>(v2);
  arma::mat MAT2 = as<arma::mat>(m2);

  arma::mat FIM = Arma_IM1 + MAT1 - MAT2;
  arma::mat F1 = Arma_D.t() * FIM * Arma_D;

  arma::mat vcov = inv(F1);
  return optMsuref_int(vcov, K_idx, K_IS_NULL, optimMeasure);
}

// === Main Function ===========================================================

/*
 * Performs the majority of actions done in kofnGA_mc.R which take place after the basic input
 * checking step has been completed. This particular function handles the case where allBetaInd (all(beta[Ind]==0) in R)
 * is true and K_idx (K.idx in R) is a vector of more than 1 indices.
 * Inputs:
 * n             The number of objects to choose from.  The algorithm chooses a subset of integers
 *               from 1 to n.
 * k             The number of objects to choose.
 * OF            The objective function.  The first argument of OF should be an index vector of k
 *               integers in the range [1, n].
 * popsize       The size of the population, the number of offspring produced each generation.
 * keepbest      This argument is used to implement elitism.  The keepbest least fit offspring each
 *               generation are replaced by the keepbest most fit members of the previous generation.
 * ngen          The number of generations to run.
 * tourneysize   The number of individuals to be selected for each tournament.
 * mutprob       The probability of mutation for each of the k chosen indices in each individual.
 *               An index chosen for mutation simply jumps to any other unused index at random.
 * pop           The population.
 * IM1, IM2_id, IM3_id, D  Additional arguments to the objective function.
 * K_idx         (Rcpp version of K.idx) A vector of more than 1 indices.
 * optimMeasure  1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Ind
 * npars
 */
// [[Rcpp::export]]
Rcpp::List kofnGAv1(int n, int k, int popsize,int keepbest,int ngen,
                    int tourneysize,float mutprob, NumericMatrix pop,
                    NumericMatrix IM1, NumericMatrix IM2_id,
                    NumericMatrix IM3_id, NumericMatrix D, NumericVector K_idx,
                    int K_IS_NULL, int optimMeasure, NumericVector Ind, int npars) {

  //=== Create needed objects =================================================

  IntegerVector indicesRcpp = Rcpp::seq_len(n);
  IntegerVector indices = indicesRcpp - 1;
  arma::urowvec indicesArma = as<arma::urowvec>(indices);
  arma::urowvec elitespotsArma;
  arma::urowvec newspotsArma;
  arma::umat Arma_pop = as<arma::umat>(pop);

  if (keepbest > 0) {
    IntegerVector elitespotsRcpp = Rcpp::seq_len(keepbest);
    IntegerVector elitespots = elitespotsRcpp - 1;
    elitespotsArma = as<arma::urowvec>(elitespots);

    IntegerVector newspotsRcpp = Rcpp::seq(keepbest+1, popsize);
    IntegerVector newspots = newspotsRcpp - 1;
    newspotsArma = as<arma::urowvec>(newspots);
  }

  arma::rowvec fitness_old = arma::rowvec(popsize);
  arma::rowvec fitness_new = arma::rowvec(popsize);
  arma::umat offspring = arma::umat(popsize, k, arma::fill::zeros);
  arma::umat old_best = arma::umat(ngen+1, k, arma::fill::zeros);
  arma::vec old_obj = arma::vec(ngen+1);
  arma::vec old_avg = arma::vec(ngen+1);
  Rcpp::List old = Rcpp::List::create(Named("best") = old_best, Named("obj") = old_obj, Named("avg") = old_avg);
  arma::umat Tourneys1 = arma::umat(popsize, tourneysize, arma::fill::zeros);
  arma::umat Tourneys2 = arma::umat(popsize, tourneysize, arma::fill::zeros);

  arma::mat Arma_D = as<arma::mat>(D);
  arma::mat Arma_IM1 = as<arma::mat>(IM1);
  arma::mat Arma_IM2_id = as<arma::mat>(IM2_id);
  arma::mat Arma_IM3_id = as<arma::mat>(IM3_id);
  arma::urowvec Arma_Ind = as<arma::urowvec>(Ind);
  arma::urowvec Arma_K_idx = as<arma::urowvec>(K_idx);


  //=== Initialize the population and evaluate their fitness ================
  for (arma::uword i = 0; i < popsize; i++) {
    fitness_old(i) = ObjFun_v1_vec(Arma_pop.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, Arma_K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);

  }

  arma::umat old_best_list = old["best"];
  NumericVector rcpp_fitness_old = wrap(fitness_old);
  NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
  urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
  arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
  // Note:
  // old_best_list.row(0) = arma::sort(Arma_pop.rows(pop_row_ids));
  // If I only use the above line, a bug occurs and the output is NOT
  // guaranteed to be sorted each time. Hence the two lines below instead:
  old_best_list.row(0) = Arma_pop.rows(pop_row_ids);
  old_best_list.row(0) = arma::sort(old_best_list.row(0));
  old["best"] = old_best_list;

  arma::vec old_obj_list = old["obj"];
  old_obj_list(0) = arma::min(fitness_old);
  old["obj"] = old_obj_list;

  arma::vec old_avg_list = old["avg"];
  old_avg_list(0) = arma::mean(fitness_old);
  old["avg"] = old_avg_list;

  //=== Loop through generations ============================================
  IntegerVector tourney_indices = Rcpp::seq(0,popsize-1);
  arma::umat tourneys1_to_transpose = arma::umat(popsize, tourneysize);
  arma::umat tourneys2_to_transpose = arma::umat(popsize, tourneysize);
  arma::urowvec Parents1 = arma::urowvec(popsize);
  arma::urowvec Parents2 = arma::urowvec(popsize);

  for (uword gen = 0; gen < ngen; gen++) {
    // Choose mating pairs (selection)

    Tourneys1 = rep_transpose(tourneys1_to_transpose, popsize, tourneysize, tourney_indices);
    Tourneys2 = rep_transpose(tourneys2_to_transpose, popsize, tourneysize, tourney_indices);

    NumericVector rcpp_fitness_old = wrap(fitness_old);
    NumericVector fitness_old_copy_rcpp = clone(rcpp_fitness_old);
    arma::rowvec fitness_old_copy = as<arma::rowvec>(fitness_old_copy_rcpp);

    for (arma::uword i = 0; i < Tourneys1.n_rows; i++) {
      // Parents1
      arma::rowvec fitness_old_sliced = fitness_old_copy.cols(Tourneys1.row(i));
      fitness_old_sliced *= -1;

      NumericVector rcpp_fitness_old_copy = wrap(fitness_old_sliced);
      NumericVector rank_rcpp_fitness_old_copy = random_rank(rcpp_fitness_old_copy);
      arma::urowvec prob = as<arma::urowvec>(rank_rcpp_fitness_old_copy);

      vec prob_vec = conv_to<vec>::from(prob);
      arma::urowvec sample_vec1 = Tourneys1.row(i);
      Parents1(i) = Rcpp::RcppArmadillo::sample(sample_vec1, 1, false, prob_vec)(0);

      // Parents2
      arma::rowvec fitness_old_sliced_T2 = fitness_old_copy.cols(Tourneys2.row(i));
      fitness_old_sliced_T2 *= -1;

      NumericVector rcpp_fitness_old_copy_T2 = wrap(fitness_old_sliced_T2);
      NumericVector rank_rcpp_fitness_old_copy_T2 = random_rank(rcpp_fitness_old_copy_T2);
      arma::urowvec prob_T2 = as<arma::urowvec>(rank_rcpp_fitness_old_copy_T2);

      vec prob_vec_T2 = conv_to<vec>::from(prob_T2);
      arma::urowvec sample_vec2 = Tourneys2.row(i);
      Parents2(i) = Rcpp::RcppArmadillo::sample(sample_vec2, 1, false, prob_vec_T2)(0);
    }

    arma::umat chosen = arma::umat(popsize, k);
    arma::vec nchosen = vec(popsize);
    for (arma::uword i = 0; i < popsize; i++) {
      urowvec combo = arma::unique(arma::join_rows(Arma_pop.row(Parents1(i)), Arma_pop.row(Parents2(i))));
      offspring.row(i) = Rcpp::RcppArmadillo::sample(combo, k, false);

      NumericVector ith_row_rcpp = Rcpp::rbinom(k, 1, mutprob);
      chosen.row(i) = as<arma::urowvec>(ith_row_rcpp);
      nchosen(i) = arma::accu(chosen.row(i));

      if (nchosen(i) > 0) {
        arma::urowvec offspring_row_copy = offspring.row(i);
        IntegerVector indices_clone = clone(indices);
        arma::urowvec arma_indices_clone = as<arma::urowvec>(indices_clone);
        arma_indices_clone.shed_cols(offspring_row_copy);
        arma::urowvec toadd = Rcpp::RcppArmadillo::sample(arma_indices_clone, nchosen(i), false);
        arma::urowvec chosen_rows_i = chosen.row(i);
        arma::uvec chosen_rows_i_inds = arma::find(chosen_rows_i == 1);
        urowvec index = {i};
        offspring.submat(index, chosen_rows_i_inds) = toadd;
      }
      fitness_new(i) = ObjFun_v1_vec(offspring.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, Arma_K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);
    }

    if (keepbest == 0) {
      Arma_pop = offspring;
      fitness_old = fitness_new;
    }
    else {
      NumericVector rcpp_fitness_old = wrap(fitness_old);
      NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
      urowvec old_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_old);
      arma::uvec old_keep = arma::find(old_keep_temp <= keepbest);

      NumericVector rcpp_fitness_new = wrap(fitness_new);
      NumericVector rank_rcpp_fitness_new = random_rank(rcpp_fitness_new);
      arma::urowvec new_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_new);
      arma::uvec new_keep = arma::find(new_keep_temp <= (popsize-keepbest));

      Arma_pop.rows(elitespotsArma) = Arma_pop.rows(old_keep);
      Arma_pop.rows(newspotsArma) = offspring.rows(new_keep);
      fitness_old.elem(elitespotsArma) = fitness_old.elem(old_keep);
      fitness_old.elem(newspotsArma) = fitness_new.elem(new_keep);

    }

    // Record progress of the search---------------------------------------------------------------
    arma::umat old_best_list = old["best"];
    rcpp_fitness_old = wrap(fitness_old);
    rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
    urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
    arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
    // Note:
    // old_best_list.row(gen+1) = arma::sort(Arma_pop.rows(pop_row_ids));
    // If I only use the above line, a bug occurs and the output is NOT
    // guaranteed to be sorted each time. Hence the two lines below instead:
    old_best_list.row(gen+1) = Arma_pop.rows(pop_row_ids);
    old_best_list.row(gen+1) = arma::sort(old_best_list.row(gen+1));

    old["best"] = old_best_list;

    arma::vec old_obj_list = old["obj"];
    old_obj_list(gen+1) = arma::min(fitness_old);
    old["obj"] = old_obj_list;

    arma::vec old_avg_list = old["avg"];
    old_avg_list(gen+1) = arma::mean(fitness_old);
    old["avg"] = old_avg_list;
  }

  //=== Package the outputs ======================================================
  //Return the final population, sorted.
  arma::urowvec ord = order_fun(fitness_old);
  IntegerVector rcpp_ord = wrap(ord);
  rcpp_ord = rcpp_ord - 1;
  ord = as<arma::urowvec>(rcpp_ord);

  arma::umat out_pop = arma::umat(popsize, k, arma::fill::zeros);
  for (arma::uword i = 0; i < popsize; i++) {
    out_pop.row(i) = arma::sort(Arma_pop.row(ord(i)));
  }

  out_pop += 1;

  // Return the objfun values for the final population.
  arma::vec out_obj = fitness_old.elem(ord);

  // Return the best solution found, and its objfun value.
  NumericVector rcpp_old_obj = wrap(old_obj);
  arma::rowvec no_na_old_obj = as<arma::rowvec>(na_omit(rcpp_old_obj));
  arma::uvec alltimebest = arma::find(old_obj == arma::min(no_na_old_obj));
  alltimebest = alltimebest.tail(1);

  arma::umat old_best_final = old["best"];
  old_best_final += 1;
  old["best"] = old_best_final;
  arma::umat sub_mat_old_best = old_best_final.rows(alltimebest);

  arma::rowvec old_obj_final = old["obj"];
  arma::rowvec sub_vec_old_obj = old_obj_final.elem(alltimebest);

  Rcpp::List out = Rcpp::List::create(Named("old") = old,
                                      Named("pop") = out_pop,
                                      Named("obj") = out_obj,
                                      Named("bestsol") = sub_mat_old_best,
                                      Named("bestobj") = sub_vec_old_obj);

  return out;
}

/*
 * Performs the majority of actions done in kofnGA_mc.R which take place after the basic input
 * checking step has been completed. This particular function handles the case where allBetaInd (all(beta[Ind]==0) in R)
 * is true and K.idx is a single index or null.
 * Inputs:
 * n             The number of objects to choose from.  The algorithm chooses a subset of integers
 *               from 1 to n.
 * k             The number of objects to choose.
 * OF            The objective function.  The first argument of OF should be an index vector of k
 *               integers in the range [1, n].
 * popsize       The size of the population, the number of offspring produced each generation.
 * keepbest      This argument is used to implement elitism.  The keepbest least fit offspring each
 *               generation are replaced by the keepbest most fit members of the previous generation.
 * ngen          The number of generations to run.
 * tourneysize   The number of individuals to be selected for each tournament.
 * mutprob       The probability of mutation for each of the k chosen indices in each individual.
 *               An index chosen for mutation simply jumps to any other unused index at random.
 * pop           The population.
 * IM1, IM2_id, IM3_id, D  Additional arguments to the objective function.
 * K_idx         K.idx is a single index or null in this case thus K_idx is -1 (K.idx is null)
 *               or a non-negative index (K.idx is a single index).
 * optimMeasure  1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Ind
 * npars
 */
// [[Rcpp::export]]
Rcpp::List kofnGAv2(int n, int k, int popsize,int keepbest,int ngen,
                    int tourneysize,float mutprob, NumericMatrix pop,
                    NumericMatrix IM1, NumericMatrix IM2_id,
                    NumericMatrix IM3_id, NumericMatrix D, int K_idx,
                    int K_IS_NULL, int optimMeasure, NumericVector Ind, int npars) {

  //=== Create needed objects =================================================

  IntegerVector indicesRcpp = Rcpp::seq_len(n);
  IntegerVector indices = indicesRcpp - 1;
  arma::urowvec indicesArma = as<arma::urowvec>(indices);
  arma::urowvec elitespotsArma;
  arma::urowvec newspotsArma;
  arma::umat Arma_pop = as<arma::umat>(pop);

  if (keepbest > 0) {
    IntegerVector elitespotsRcpp = Rcpp::seq_len(keepbest);
    IntegerVector elitespots = elitespotsRcpp - 1;
    elitespotsArma = as<arma::urowvec>(elitespots);

    IntegerVector newspotsRcpp = Rcpp::seq(keepbest+1, popsize);
    IntegerVector newspots = newspotsRcpp - 1;
    newspotsArma = as<arma::urowvec>(newspots);
  }

  arma::rowvec fitness_old = arma::rowvec(popsize);
  arma::rowvec fitness_new = arma::rowvec(popsize);
  arma::umat offspring = arma::umat(popsize, k, arma::fill::zeros);
  arma::umat old_best = arma::umat(ngen+1, k, arma::fill::zeros);
  arma::vec old_obj = arma::vec(ngen+1);
  arma::vec old_avg = arma::vec(ngen+1);
  Rcpp::List old = Rcpp::List::create(Named("best") = old_best, Named("obj") = old_obj, Named("avg") = old_avg);
  arma::umat Tourneys1 = arma::umat(popsize, tourneysize, arma::fill::zeros);
  arma::umat Tourneys2 = arma::umat(popsize, tourneysize, arma::fill::zeros);

  arma::mat Arma_D = as<arma::mat>(D);
  arma::mat Arma_IM1 = as<arma::mat>(IM1);
  arma::mat Arma_IM2_id = as<arma::mat>(IM2_id);
  arma::mat Arma_IM3_id = as<arma::mat>(IM3_id);
  arma::urowvec Arma_Ind = as<arma::urowvec>(Ind);

  //=== Initialize the population and evaluate their fitness ================
  for (arma::uword i = 0; i < popsize; i++) {
    fitness_old(i) = ObjFun_v1_int(Arma_pop.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);

  }

  arma::umat old_best_list = old["best"];
  NumericVector rcpp_fitness_old = wrap(fitness_old);
  NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
  urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
  arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
  // Note:
  // old_best_list.row(0) = arma::sort(Arma_pop.rows(pop_row_ids));
  // If I only use the above line, a bug occurs and the output is NOT
  // guaranteed to be sorted each time. Hence the two lines below instead:
  old_best_list.row(0) = Arma_pop.rows(pop_row_ids);
  old_best_list.row(0) = arma::sort(old_best_list.row(0));
  old["best"] = old_best_list;

  arma::vec old_obj_list = old["obj"];
  old_obj_list(0) = arma::min(fitness_old);
  old["obj"] = old_obj_list;

  arma::vec old_avg_list = old["avg"];
  old_avg_list(0) = arma::mean(fitness_old);
  old["avg"] = old_avg_list;

  //=== Loop through generations ============================================
  IntegerVector tourney_indices = Rcpp::seq(0,popsize-1);
  arma::umat tourneys1_to_transpose = arma::umat(popsize, tourneysize);
  arma::umat tourneys2_to_transpose = arma::umat(popsize, tourneysize);
  arma::urowvec Parents1 = arma::urowvec(popsize);
  arma::urowvec Parents2 = arma::urowvec(popsize);

  for (uword gen = 0; gen < ngen; gen++) {
    // Choose mating pairs (selection)

    Tourneys1 = rep_transpose(tourneys1_to_transpose, popsize, tourneysize, tourney_indices);
    Tourneys2 = rep_transpose(tourneys2_to_transpose, popsize, tourneysize, tourney_indices);

    NumericVector rcpp_fitness_old = wrap(fitness_old);
    NumericVector fitness_old_copy_rcpp = clone(rcpp_fitness_old);
    arma::rowvec fitness_old_copy = as<arma::rowvec>(fitness_old_copy_rcpp);

    for (arma::uword i = 0; i < Tourneys1.n_rows; i++) {
      // Parents1
      arma::rowvec fitness_old_sliced = fitness_old_copy.cols(Tourneys1.row(i));
      fitness_old_sliced *= -1;

      NumericVector rcpp_fitness_old_copy = wrap(fitness_old_sliced);
      NumericVector rank_rcpp_fitness_old_copy = random_rank(rcpp_fitness_old_copy);
      arma::urowvec prob = as<arma::urowvec>(rank_rcpp_fitness_old_copy);

      vec prob_vec = conv_to<vec>::from(prob);
      arma::urowvec sample_vec1 = Tourneys1.row(i);
      Parents1(i) = Rcpp::RcppArmadillo::sample(sample_vec1, 1, false, prob_vec)(0);

      // Parents2
      arma::rowvec fitness_old_sliced_T2 = fitness_old_copy.cols(Tourneys2.row(i));
      fitness_old_sliced_T2 *= -1;

      NumericVector rcpp_fitness_old_copy_T2 = wrap(fitness_old_sliced_T2);
      NumericVector rank_rcpp_fitness_old_copy_T2 = random_rank(rcpp_fitness_old_copy_T2);
      arma::urowvec prob_T2 = as<arma::urowvec>(rank_rcpp_fitness_old_copy_T2);

      vec prob_vec_T2 = conv_to<vec>::from(prob_T2);
      arma::urowvec sample_vec2 = Tourneys2.row(i);
      Parents2(i) = Rcpp::RcppArmadillo::sample(sample_vec2, 1, false, prob_vec_T2)(0);
    }

    arma::umat chosen = arma::umat(popsize, k);
    arma::vec nchosen = vec(popsize);
    for (arma::uword i = 0; i < popsize; i++) {
      urowvec combo = arma::unique(arma::join_rows(Arma_pop.row(Parents1(i)), Arma_pop.row(Parents2(i))));
      offspring.row(i) = Rcpp::RcppArmadillo::sample(combo, k, false);

      NumericVector ith_row_rcpp = Rcpp::rbinom(k, 1, mutprob);
      chosen.row(i) = as<arma::urowvec>(ith_row_rcpp);
      nchosen(i) = arma::accu(chosen.row(i));

      if (nchosen(i) > 0) {
        arma::urowvec offspring_row_copy = offspring.row(i);
        IntegerVector indices_clone = clone(indices);
        arma::urowvec arma_indices_clone = as<arma::urowvec>(indices_clone);
        arma_indices_clone.shed_cols(offspring_row_copy);
        arma::urowvec toadd = Rcpp::RcppArmadillo::sample(arma_indices_clone, nchosen(i), false);
        arma::urowvec chosen_rows_i = chosen.row(i);
        arma::uvec chosen_rows_i_inds = arma::find(chosen_rows_i == 1);
        urowvec index = {i};
        offspring.submat(index, chosen_rows_i_inds) = toadd;
      }
      fitness_new(i) = ObjFun_v1_int(offspring.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);
    }

    if (keepbest == 0) {
      Arma_pop = offspring;
      fitness_old = fitness_new;
    }
    else {
      NumericVector rcpp_fitness_old = wrap(fitness_old);
      NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
      urowvec old_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_old);
      arma::uvec old_keep = arma::find(old_keep_temp <= keepbest);

      NumericVector rcpp_fitness_new = wrap(fitness_new);
      NumericVector rank_rcpp_fitness_new = random_rank(rcpp_fitness_new);
      arma::urowvec new_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_new);
      arma::uvec new_keep = arma::find(new_keep_temp <= (popsize-keepbest));

      Arma_pop.rows(elitespotsArma) = Arma_pop.rows(old_keep);
      Arma_pop.rows(newspotsArma) = offspring.rows(new_keep);
      fitness_old.elem(elitespotsArma) = fitness_old.elem(old_keep);
      fitness_old.elem(newspotsArma) = fitness_new.elem(new_keep);

    }

    // Record progress of the search---------------------------------------------------------------
    arma::umat old_best_list = old["best"];
    rcpp_fitness_old = wrap(fitness_old);
    rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
    urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
    arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
    // Note:
    // old_best_list.row(gen+1) = arma::sort(Arma_pop.rows(pop_row_ids));
    // If I only use the above line, a bug occurs and the output is NOT
    // guaranteed to be sorted each time. Hence the two lines below instead:
    old_best_list.row(gen+1) = Arma_pop.rows(pop_row_ids);
    old_best_list.row(gen+1) = arma::sort(old_best_list.row(gen+1));

    old["best"] = old_best_list;

    arma::vec old_obj_list = old["obj"];
    old_obj_list(gen+1) = arma::min(fitness_old);
    old["obj"] = old_obj_list;

    arma::vec old_avg_list = old["avg"];
    old_avg_list(gen+1) = arma::mean(fitness_old);
    old["avg"] = old_avg_list;
  }

  //=== Package the outputs ======================================================
  //Return the final population, sorted.
  arma::urowvec ord = order_fun(fitness_old);
  IntegerVector rcpp_ord = wrap(ord);
  rcpp_ord = rcpp_ord - 1;
  ord = as<arma::urowvec>(rcpp_ord);

  arma::umat out_pop = arma::umat(popsize, k, arma::fill::zeros);
  for (arma::uword i = 0; i < popsize; i++) {
    out_pop.row(i) = arma::sort(Arma_pop.row(ord(i)));
  }

  out_pop += 1;

  // Return the objfun values for the final population.
  arma::vec out_obj = fitness_old.elem(ord);

  // Return the best solution found, and its objfun value.
  NumericVector rcpp_old_obj = wrap(old_obj);
  arma::rowvec no_na_old_obj = as<arma::rowvec>(na_omit(rcpp_old_obj));
  arma::uvec alltimebest = arma::find(old_obj == arma::min(no_na_old_obj));
  alltimebest = alltimebest.tail(1);

  arma::umat old_best_final = old["best"];
  old_best_final += 1;
  old["best"] = old_best_final;
  arma::umat sub_mat_old_best = old_best_final.rows(alltimebest);

  arma::rowvec old_obj_final = old["obj"];
  arma::rowvec sub_vec_old_obj = old_obj_final.elem(alltimebest);

  Rcpp::List out = Rcpp::List::create(Named("old") = old,
                                      Named("pop") = out_pop,
                                      Named("obj") = out_obj,
                                      Named("bestsol") = sub_mat_old_best,
                                      Named("bestobj") = sub_vec_old_obj);

  return out;
}

/*
 * Performs the majority of actions done in kofnGA_mc.R which take place after the basic input
 * checking step has been completed. This particular function handles the case where allBetaInd (all(beta[Ind]==0) in R)
 * is false and K_idx (K.idx in R) is a vector of more than 1 indices.
 * Inputs:
 * n             The number of objects to choose from.  The algorithm chooses a subset of integers
 *               from 1 to n.
 * k             The number of objects to choose.
 * OF            The objective function.  The first argument of OF should be an index vector of k
 *               integers in the range [1, n].
 * popsize       The size of the population, the number of offspring produced each generation.
 * keepbest      This argument is used to implement elitism.  The keepbest least fit offspring each
 *               generation are replaced by the keepbest most fit members of the previous generation.
 * ngen          The number of generations to run.
 * tourneysize   The number of individuals to be selected for each tournament.
 * mutprob       The probability of mutation for each of the k chosen indices in each individual.
 *               An index chosen for mutation simply jumps to any other unused index at random.
 * pop           The population.
 * IM1, IM2_id, IM3_id, D  Additional arguments to the objective function.
 * K_idx         (Rcpp version of K.idx) A vector of more than 1 indices.
 * optimMeasure  1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Ind
 * npars
 */
// [[Rcpp::export]]
Rcpp::List kofnGAv3(int n, int k, int popsize,int keepbest,int ngen,
                    int tourneysize,float mutprob, NumericMatrix pop,
                    NumericMatrix IM1, NumericMatrix IM2_id,
                    NumericMatrix IM3_id, NumericMatrix D, NumericVector K_idx,
                    int K_IS_NULL, int optimMeasure, NumericVector Ind, int npars) {

  //=== Create needed objects =================================================

  IntegerVector indicesRcpp = Rcpp::seq_len(n);
  IntegerVector indices = indicesRcpp - 1;
  arma::urowvec indicesArma = as<arma::urowvec>(indices);
  arma::urowvec elitespotsArma;
  arma::urowvec newspotsArma;
  arma::umat Arma_pop = as<arma::umat>(pop);

  if (keepbest > 0) {
    IntegerVector elitespotsRcpp = Rcpp::seq_len(keepbest);
    IntegerVector elitespots = elitespotsRcpp - 1;
    elitespotsArma = as<arma::urowvec>(elitespots);

    IntegerVector newspotsRcpp = Rcpp::seq(keepbest+1, popsize);
    IntegerVector newspots = newspotsRcpp - 1;
    newspotsArma = as<arma::urowvec>(newspots);
  }

  arma::rowvec fitness_old = arma::rowvec(popsize);
  arma::rowvec fitness_new = arma::rowvec(popsize);
  arma::umat offspring = arma::umat(popsize, k, arma::fill::zeros);
  arma::umat old_best = arma::umat(ngen+1, k, arma::fill::zeros);
  arma::vec old_obj = arma::vec(ngen+1);
  arma::vec old_avg = arma::vec(ngen+1);
  Rcpp::List old = Rcpp::List::create(Named("best") = old_best, Named("obj") = old_obj, Named("avg") = old_avg);
  arma::umat Tourneys1 = arma::umat(popsize, tourneysize, arma::fill::zeros);
  arma::umat Tourneys2 = arma::umat(popsize, tourneysize, arma::fill::zeros);

  arma::mat Arma_D = as<arma::mat>(D);
  arma::mat Arma_IM1 = as<arma::mat>(IM1);
  arma::mat Arma_IM2_id = as<arma::mat>(IM2_id);
  arma::mat Arma_IM3_id = as<arma::mat>(IM3_id);
  arma::urowvec Arma_Ind = as<arma::urowvec>(Ind);
  arma::urowvec Arma_K_idx = as<arma::urowvec>(K_idx);

  //=== Initialize the population and evaluate their fitness ================
  for (arma::uword i = 0; i < popsize; i++) {
    fitness_old(i) = ObjFun_v2_vec(Arma_pop.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, Arma_K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);

  }

  arma::umat old_best_list = old["best"];
  NumericVector rcpp_fitness_old = wrap(fitness_old);
  NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
  urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
  arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
  // Note:
  // old_best_list.row(0) = arma::sort(Arma_pop.rows(pop_row_ids));
  // If I only use the above line, a bug occurs and the output is NOT
  // guaranteed to be sorted each time. Hence the two lines below instead:
  old_best_list.row(0) = Arma_pop.rows(pop_row_ids);
  old_best_list.row(0) = arma::sort(old_best_list.row(0));
  old["best"] = old_best_list;

  arma::vec old_obj_list = old["obj"];
  old_obj_list(0) = arma::min(fitness_old);
  old["obj"] = old_obj_list;

  arma::vec old_avg_list = old["avg"];
  old_avg_list(0) = arma::mean(fitness_old);
  old["avg"] = old_avg_list;

  //=== Loop through generations ============================================
  IntegerVector tourney_indices = Rcpp::seq(0,popsize-1);
  arma::umat tourneys1_to_transpose = arma::umat(popsize, tourneysize);
  arma::umat tourneys2_to_transpose = arma::umat(popsize, tourneysize);
  arma::urowvec Parents1 = arma::urowvec(popsize);
  arma::urowvec Parents2 = arma::urowvec(popsize);

  for (uword gen = 0; gen < ngen; gen++) {
    // Choose mating pairs (selection)

    Tourneys1 = rep_transpose(tourneys1_to_transpose, popsize, tourneysize, tourney_indices);
    Tourneys2 = rep_transpose(tourneys2_to_transpose, popsize, tourneysize, tourney_indices);

    NumericVector rcpp_fitness_old = wrap(fitness_old);
    NumericVector fitness_old_copy_rcpp = clone(rcpp_fitness_old);
    arma::rowvec fitness_old_copy = as<arma::rowvec>(fitness_old_copy_rcpp);

    for (arma::uword i = 0; i < Tourneys1.n_rows; i++) {
      // Parents1
      arma::rowvec fitness_old_sliced = fitness_old_copy.cols(Tourneys1.row(i));
      fitness_old_sliced *= -1;

      NumericVector rcpp_fitness_old_copy = wrap(fitness_old_sliced);
      NumericVector rank_rcpp_fitness_old_copy = random_rank(rcpp_fitness_old_copy);
      arma::urowvec prob = as<arma::urowvec>(rank_rcpp_fitness_old_copy);

      vec prob_vec = conv_to<vec>::from(prob);
      arma::urowvec sample_vec1 = Tourneys1.row(i);
      Parents1(i) = Rcpp::RcppArmadillo::sample(sample_vec1, 1, false, prob_vec)(0);

      // Parents2
      arma::rowvec fitness_old_sliced_T2 = fitness_old_copy.cols(Tourneys2.row(i));
      fitness_old_sliced_T2 *= -1;

      NumericVector rcpp_fitness_old_copy_T2 = wrap(fitness_old_sliced_T2);
      NumericVector rank_rcpp_fitness_old_copy_T2 = random_rank(rcpp_fitness_old_copy_T2);
      arma::urowvec prob_T2 = as<arma::urowvec>(rank_rcpp_fitness_old_copy_T2);

      vec prob_vec_T2 = conv_to<vec>::from(prob_T2);
      arma::urowvec sample_vec2 = Tourneys2.row(i);
      Parents2(i) = Rcpp::RcppArmadillo::sample(sample_vec2, 1, false, prob_vec_T2)(0);
    }

    arma::umat chosen = arma::umat(popsize, k);
    arma::vec nchosen = vec(popsize);
    for (arma::uword i = 0; i < popsize; i++) {
      urowvec combo = arma::unique(arma::join_rows(Arma_pop.row(Parents1(i)), Arma_pop.row(Parents2(i))));
      offspring.row(i) = Rcpp::RcppArmadillo::sample(combo, k, false);

      NumericVector ith_row_rcpp = Rcpp::rbinom(k, 1, mutprob);
      chosen.row(i) = as<arma::urowvec>(ith_row_rcpp);
      nchosen(i) = arma::accu(chosen.row(i));

      if (nchosen(i) > 0) {
        arma::urowvec offspring_row_copy = offspring.row(i);
        IntegerVector indices_clone = clone(indices);
        arma::urowvec arma_indices_clone = as<arma::urowvec>(indices_clone);
        arma_indices_clone.shed_cols(offspring_row_copy);
        arma::urowvec toadd = Rcpp::RcppArmadillo::sample(arma_indices_clone, nchosen(i), false);
        arma::urowvec chosen_rows_i = chosen.row(i);
        arma::uvec chosen_rows_i_inds = arma::find(chosen_rows_i == 1);
        urowvec index = {i};
        offspring.submat(index, chosen_rows_i_inds) = toadd;
      }
      fitness_new(i) = ObjFun_v2_vec(offspring.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, Arma_K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);
    }

    if (keepbest == 0) {
      Arma_pop = offspring;
      fitness_old = fitness_new;
    }
    else {
      NumericVector rcpp_fitness_old = wrap(fitness_old);
      NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
      urowvec old_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_old);
      arma::uvec old_keep = arma::find(old_keep_temp <= keepbest);

      NumericVector rcpp_fitness_new = wrap(fitness_new);
      NumericVector rank_rcpp_fitness_new = random_rank(rcpp_fitness_new);
      arma::urowvec new_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_new);
      arma::uvec new_keep = arma::find(new_keep_temp <= (popsize-keepbest));

      Arma_pop.rows(elitespotsArma) = Arma_pop.rows(old_keep);
      Arma_pop.rows(newspotsArma) = offspring.rows(new_keep);
      fitness_old.elem(elitespotsArma) = fitness_old.elem(old_keep);
      fitness_old.elem(newspotsArma) = fitness_new.elem(new_keep);

    }

    // Record progress of the search---------------------------------------------------------------
    arma::umat old_best_list = old["best"];
    rcpp_fitness_old = wrap(fitness_old);
    rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
    urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
    arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
    // Note:
    // old_best_list.row(gen+1) = arma::sort(Arma_pop.rows(pop_row_ids));
    // If I only use the above line, a bug occurs and the output is NOT
    // guaranteed to be sorted each time. Hence the two lines below instead:
    old_best_list.row(gen+1) = Arma_pop.rows(pop_row_ids);
    old_best_list.row(gen+1) = arma::sort(old_best_list.row(gen+1));

    old["best"] = old_best_list;

    arma::vec old_obj_list = old["obj"];
    old_obj_list(gen+1) = arma::min(fitness_old);
    old["obj"] = old_obj_list;

    arma::vec old_avg_list = old["avg"];
    old_avg_list(gen+1) = arma::mean(fitness_old);
    old["avg"] = old_avg_list;
  }

  //=== Package the outputs ======================================================
  //Return the final population, sorted.
  arma::urowvec ord = order_fun(fitness_old);
  IntegerVector rcpp_ord = wrap(ord);
  rcpp_ord = rcpp_ord - 1;
  ord = as<arma::urowvec>(rcpp_ord);

  arma::umat out_pop = arma::umat(popsize, k, arma::fill::zeros);
  for (arma::uword i = 0; i < popsize; i++) {
    out_pop.row(i) = arma::sort(Arma_pop.row(ord(i)));
  }

  out_pop += 1;

  // Return the objfun values for the final population.
  arma::vec out_obj = fitness_old.elem(ord);

  // Return the best solution found, and its objfun value.
  NumericVector rcpp_old_obj = wrap(old_obj);
  arma::rowvec no_na_old_obj = as<arma::rowvec>(na_omit(rcpp_old_obj));
  arma::uvec alltimebest = arma::find(old_obj == arma::min(no_na_old_obj));
  alltimebest = alltimebest.tail(1);

  arma::umat old_best_final = old["best"];
  old_best_final += 1;
  old["best"] = old_best_final;
  arma::umat sub_mat_old_best = old_best_final.rows(alltimebest);

  arma::rowvec old_obj_final = old["obj"];
  arma::rowvec sub_vec_old_obj = old_obj_final.elem(alltimebest);

  Rcpp::List out = Rcpp::List::create(Named("old") = old,
                                      Named("pop") = out_pop,
                                      Named("obj") = out_obj,
                                      Named("bestsol") = sub_mat_old_best,
                                      Named("bestobj") = sub_vec_old_obj);

  return out;
}

/*
 * Performs the majority of actions done in kofnGA_mc.R which take place after the basic input
 * checking step has been completed. This particular function handles the case where allBetaInd (all(beta[Ind]==0) in R)
 * is false and K.idx is a single index or null.
 * Inputs:
 * n             The number of objects to choose from.  The algorithm chooses a subset of integers
 *               from 1 to n.
 * k             The number of objects to choose.
 * OF            The objective function.  The first argument of OF should be an index vector of k
 *               integers in the range [1, n].
 * popsize       The size of the population, the number of offspring produced each generation.
 * keepbest      This argument is used to implement elitism.  The keepbest least fit offspring each
 *               generation are replaced by the keepbest most fit members of the previous generation.
 * ngen          The number of generations to run.
 * tourneysize   The number of individuals to be selected for each tournament.
 * mutprob       The probability of mutation for each of the k chosen indices in each individual.
 *               An index chosen for mutation simply jumps to any other unused index at random.
 * pop           The population.
 * IM1, IM2_id, IM3_id, D  Additional arguments to the objective function.
 * K_idx         K.idx is a single index or null in this case thus K_idx is -1 (K.idx is null)
 *               or a non-negative index (K.idx is a single index).
 * optimMeasure  1 = D-opt, 2 = A-opt, 3 = Par-spec.
 * Ind
 * npars
 */
// [[Rcpp::export]]
Rcpp::List kofnGAv4(int n, int k, int popsize,int keepbest,int ngen,
                    int tourneysize,float mutprob, NumericMatrix pop,
                    NumericMatrix IM1, NumericMatrix IM2_id,
                    NumericMatrix IM3_id, NumericMatrix D, int K_idx,
                    int K_IS_NULL, int optimMeasure, NumericVector Ind, int npars) {

  //=== Create needed objects =================================================

  IntegerVector indicesRcpp = Rcpp::seq_len(n);
  IntegerVector indices = indicesRcpp - 1;
  arma::urowvec indicesArma = as<arma::urowvec>(indices);
  arma::urowvec elitespotsArma;
  arma::urowvec newspotsArma;
  arma::umat Arma_pop = as<arma::umat>(pop);

  if (keepbest > 0) {
    IntegerVector elitespotsRcpp = Rcpp::seq_len(keepbest);
    IntegerVector elitespots = elitespotsRcpp - 1;
    elitespotsArma = as<arma::urowvec>(elitespots);

    IntegerVector newspotsRcpp = Rcpp::seq(keepbest+1, popsize);
    IntegerVector newspots = newspotsRcpp - 1;
    newspotsArma = as<arma::urowvec>(newspots);
  }

  arma::rowvec fitness_old = arma::rowvec(popsize);
  arma::rowvec fitness_new = arma::rowvec(popsize);
  arma::umat offspring = arma::umat(popsize, k, arma::fill::zeros);
  arma::umat old_best = arma::umat(ngen+1, k, arma::fill::zeros);
  arma::vec old_obj = arma::vec(ngen+1);
  arma::vec old_avg = arma::vec(ngen+1);
  Rcpp::List old = Rcpp::List::create(Named("best") = old_best, Named("obj") = old_obj, Named("avg") = old_avg);
  arma::umat Tourneys1 = arma::umat(popsize, tourneysize, arma::fill::zeros);
  arma::umat Tourneys2 = arma::umat(popsize, tourneysize, arma::fill::zeros);

  arma::mat Arma_D = as<arma::mat>(D);
  arma::mat Arma_IM1 = as<arma::mat>(IM1);
  arma::mat Arma_IM2_id = as<arma::mat>(IM2_id);
  arma::mat Arma_IM3_id = as<arma::mat>(IM3_id);
  arma::urowvec Arma_Ind = as<arma::urowvec>(Ind);

  //=== Initialize the population and evaluate their fitness ================
  for (arma::uword i = 0; i < popsize; i++) {
    fitness_old(i) = ObjFun_v2_int(Arma_pop.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);

  }

  arma::umat old_best_list = old["best"];
  NumericVector rcpp_fitness_old = wrap(fitness_old);
  NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
  urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
  arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
  // Note:
  // old_best_list.row(0) = arma::sort(Arma_pop.rows(pop_row_ids));
  // If I only use the above line, a bug occurs and the output is NOT
  // guaranteed to be sorted each time. Hence the two lines below instead:
  old_best_list.row(0) = Arma_pop.rows(pop_row_ids);
  old_best_list.row(0) = arma::sort(old_best_list.row(0));
  old["best"] = old_best_list;

  arma::vec old_obj_list = old["obj"];
  old_obj_list(0) = arma::min(fitness_old);
  old["obj"] = old_obj_list;

  arma::vec old_avg_list = old["avg"];
  old_avg_list(0) = arma::mean(fitness_old);
  old["avg"] = old_avg_list;

  //=== Loop through generations ============================================
  IntegerVector tourney_indices = Rcpp::seq(0,popsize-1);
  arma::umat tourneys1_to_transpose = arma::umat(popsize, tourneysize);
  arma::umat tourneys2_to_transpose = arma::umat(popsize, tourneysize);
  arma::urowvec Parents1 = arma::urowvec(popsize);
  arma::urowvec Parents2 = arma::urowvec(popsize);

  for (uword gen = 0; gen < ngen; gen++) {
    // Choose mating pairs (selection)

    Tourneys1 = rep_transpose(tourneys1_to_transpose, popsize, tourneysize, tourney_indices);
    Tourneys2 = rep_transpose(tourneys2_to_transpose, popsize, tourneysize, tourney_indices);

    NumericVector rcpp_fitness_old = wrap(fitness_old);
    NumericVector fitness_old_copy_rcpp = clone(rcpp_fitness_old);
    arma::rowvec fitness_old_copy = as<arma::rowvec>(fitness_old_copy_rcpp);

    for (arma::uword i = 0; i < Tourneys1.n_rows; i++) {
      // Parents1
      arma::rowvec fitness_old_sliced = fitness_old_copy.cols(Tourneys1.row(i));
      fitness_old_sliced *= -1;

      NumericVector rcpp_fitness_old_copy = wrap(fitness_old_sliced);
      NumericVector rank_rcpp_fitness_old_copy = random_rank(rcpp_fitness_old_copy);
      arma::urowvec prob = as<arma::urowvec>(rank_rcpp_fitness_old_copy);

      vec prob_vec = conv_to<vec>::from(prob);
      arma::urowvec sample_vec1 = Tourneys1.row(i);
      Parents1(i) = Rcpp::RcppArmadillo::sample(sample_vec1, 1, false, prob_vec)(0);

      // Parents2
      arma::rowvec fitness_old_sliced_T2 = fitness_old_copy.cols(Tourneys2.row(i));
      fitness_old_sliced_T2 *= -1;

      NumericVector rcpp_fitness_old_copy_T2 = wrap(fitness_old_sliced_T2);
      NumericVector rank_rcpp_fitness_old_copy_T2 = random_rank(rcpp_fitness_old_copy_T2);
      arma::urowvec prob_T2 = as<arma::urowvec>(rank_rcpp_fitness_old_copy_T2);

      vec prob_vec_T2 = conv_to<vec>::from(prob_T2);
      arma::urowvec sample_vec2 = Tourneys2.row(i);
      Parents2(i) = Rcpp::RcppArmadillo::sample(sample_vec2, 1, false, prob_vec_T2)(0);
    }

    arma::umat chosen = arma::umat(popsize, k);
    arma::vec nchosen = vec(popsize);
    for (arma::uword i = 0; i < popsize; i++) {
      urowvec combo = arma::unique(arma::join_rows(Arma_pop.row(Parents1(i)), Arma_pop.row(Parents2(i))));
      offspring.row(i) = Rcpp::RcppArmadillo::sample(combo, k, false);

      NumericVector ith_row_rcpp = Rcpp::rbinom(k, 1, mutprob);
      chosen.row(i) = as<arma::urowvec>(ith_row_rcpp);
      nchosen(i) = arma::accu(chosen.row(i));

      if (nchosen(i) > 0) {
        arma::urowvec offspring_row_copy = offspring.row(i);
        IntegerVector indices_clone = clone(indices);
        arma::urowvec arma_indices_clone = as<arma::urowvec>(indices_clone);
        arma_indices_clone.shed_cols(offspring_row_copy);
        arma::urowvec toadd = Rcpp::RcppArmadillo::sample(arma_indices_clone, nchosen(i), false);
        arma::urowvec chosen_rows_i = chosen.row(i);
        arma::uvec chosen_rows_i_inds = arma::find(chosen_rows_i == 1);
        urowvec index = {i};
        offspring.submat(index, chosen_rows_i_inds) = toadd;
      }
      fitness_new(i) = ObjFun_v2_int(offspring.row(i), Arma_IM1, Arma_IM2_id, Arma_IM3_id, Arma_D, K_idx, K_IS_NULL, optimMeasure, Arma_Ind, npars);
    }

    if (keepbest == 0) {
      Arma_pop = offspring;
      fitness_old = fitness_new;
    }
    else {
      NumericVector rcpp_fitness_old = wrap(fitness_old);
      NumericVector rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
      urowvec old_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_old);
      arma::uvec old_keep = arma::find(old_keep_temp <= keepbest);

      NumericVector rcpp_fitness_new = wrap(fitness_new);
      NumericVector rank_rcpp_fitness_new = random_rank(rcpp_fitness_new);
      arma::urowvec new_keep_temp = as<arma::urowvec>(rank_rcpp_fitness_new);
      arma::uvec new_keep = arma::find(new_keep_temp <= (popsize-keepbest));

      Arma_pop.rows(elitespotsArma) = Arma_pop.rows(old_keep);
      Arma_pop.rows(newspotsArma) = offspring.rows(new_keep);
      fitness_old.elem(elitespotsArma) = fitness_old.elem(old_keep);
      fitness_old.elem(newspotsArma) = fitness_new.elem(new_keep);

    }

    // Record progress of the search---------------------------------------------------------------
    arma::umat old_best_list = old["best"];
    rcpp_fitness_old = wrap(fitness_old);
    rank_rcpp_fitness_old = random_rank(rcpp_fitness_old);
    urowvec rank_fitness_old = as<arma::urowvec>(rank_rcpp_fitness_old);
    arma::uvec pop_row_ids = arma::find(rank_fitness_old == 1);
    // Note:
    // old_best_list.row(gen+1) = arma::sort(Arma_pop.rows(pop_row_ids));
    // If I only use the above line, a bug occurs and the output is NOT
    // guaranteed to be sorted each time. Hence the two lines below instead:
    old_best_list.row(gen+1) = Arma_pop.rows(pop_row_ids);
    old_best_list.row(gen+1) = arma::sort(old_best_list.row(gen+1));

    old["best"] = old_best_list;

    arma::vec old_obj_list = old["obj"];
    old_obj_list(gen+1) = arma::min(fitness_old);
    old["obj"] = old_obj_list;

    arma::vec old_avg_list = old["avg"];
    old_avg_list(gen+1) = arma::mean(fitness_old);
    old["avg"] = old_avg_list;
  }

  //=== Package the outputs ======================================================
  //Return the final population, sorted.
  arma::urowvec ord = order_fun(fitness_old);
  IntegerVector rcpp_ord = wrap(ord);
  rcpp_ord = rcpp_ord - 1;
  ord = as<arma::urowvec>(rcpp_ord);

  arma::umat out_pop = arma::umat(popsize, k, arma::fill::zeros);
  for (arma::uword i = 0; i < popsize; i++) {
    out_pop.row(i) = arma::sort(Arma_pop.row(ord(i)));
  }

  out_pop += 1;

  // Return the objfun values for the final population.
  arma::vec out_obj = fitness_old.elem(ord);

  // Return the best solution found, and its objfun value.
  NumericVector rcpp_old_obj = wrap(old_obj);
  arma::rowvec no_na_old_obj = as<arma::rowvec>(na_omit(rcpp_old_obj));
  arma::uvec alltimebest = arma::find(old_obj == arma::min(no_na_old_obj));
  alltimebest = alltimebest.tail(1);

  arma::umat old_best_final = old["best"];
  old_best_final += 1;
  old["best"] = old_best_final;
  arma::umat sub_mat_old_best = old_best_final.rows(alltimebest);

  arma::rowvec old_obj_final = old["obj"];
  arma::rowvec sub_vec_old_obj = old_obj_final.elem(alltimebest);

  Rcpp::List out = Rcpp::List::create(Named("old") = old,
                                      Named("pop") = out_pop,
                                      Named("obj") = out_obj,
                                      Named("bestsol") = sub_mat_old_best,
                                      Named("bestobj") = sub_vec_old_obj);

  return out;
}
