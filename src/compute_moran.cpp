#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include <algorithm>

double moran_numerator(const double* exprs, const Rcpp::IntegerVector& indices, const Rcpp::NumericVector& weights, const Rcpp::IntegerVector& runs) {
    auto iIt=indices.begin();
    auto wIt=weights.begin();

    auto copy = exprs;
    double val = 0;

    for (auto r : runs) {
        double summed = 0;
        for (int i = 0; i < r; ++i, ++iIt, ++wIt) {
            summed += *wIt * exprs[*iIt - 1];
        }
        val += summed * (*copy);
        ++copy;
    }

    Rprintf("%.5f\n", val);
    return val;
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_moran_part(Rcpp::RObject input, Rcpp::IntegerVector indices, Rcpp::NumericVector weights, Rcpp::IntegerVector runs) {
    auto mat = beachmat::read_lin_block(input);
    const size_t ncells = mat->get_ncol();
    const size_t ngenes = mat->get_nrow();

    if (ncells!=runs.size()) {
        throw std::runtime_error("'runs' must be of length equal to 'ncol(input)'");
    }

    Rcpp::NumericVector stat(ngenes), squares(ngenes), quads(ngenes);
    Rcpp::NumericVector work(ncells);

    for (size_t g=0; g<ngenes; ++g) {
        auto ptr = mat->get_row(g, work.begin());
        if (ptr != work.begin()) {
            std::copy(ptr, ptr+ncells, work.begin());
        }
   
        // Compute mean and subtract it from the expression.
        const double mean = std::accumulate(work.begin(), work.end(), 0.0) / ncells;
        double& sqdiff = squares[g];
        double& quadiff = quads[g];
        for (auto& w : work) {
            w -= mean;
            sqdiff += w * w;
            quadiff += w * w * w * w;
        }

        if (sqdiff==0) {
            stat[g] = R_NaReal;
        } else {
            stat[g] = moran_numerator(work.begin(), indices, weights, runs) / sqdiff;
        }
    }

    return Rcpp::List::create(stat, squares, quads);
}
