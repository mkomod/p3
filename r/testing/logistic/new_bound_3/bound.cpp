#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double tll(arma::vec mu, arma::vec sig, int p) 
{
    const int n = mu.n_rows;
    double res = 0.0;

    for (int i = 0; i < n; ++i) 
    {
	double a = sig(i) / sqrt(2.0 * M_PI) * 
	    exp(- 0.5 * mu(i)*mu(i) / (sig(i)*sig(i))) +
	    mu(i) * R::pnorm(mu(i) / sig(i), 0, 1, 1, 0);

	double b = 0.0;
	for(int j = 1; j <= (2 * p - 1); ++j) {
	    b += pow((-1.0), (j-1)) / j * exp(
		    mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(-mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
	    ) +
	    pow((-1.0), (j-1)) / j * exp(
		    -mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
	    );
	}

	res += a + b;
    }
    
    return res;
}

