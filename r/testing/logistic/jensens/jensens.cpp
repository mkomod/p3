#include "RcppEnsmallen.h"

// [[Rcpp::depends(RcppEnsmallen)]]

using namespace arma;

// [[Rcpp::export]]
vec mvnMGF(const mat &X, const vec &mu, const vec &sig) 
{
    return exp(X * mu + 0.5 * (X % X) * (sig % sig));
}


class jen_update_mu_fn
{
    public:
	jen_update_mu_fn(const vec &y, const mat &X, const vec &mu,
		const vec &s, const double lambda, const uvec &G, 
		const vec &P) :
	    y(y), X(X), mu(mu), s(s), lambda(lambda), G(G), P(P)
	{};

	double EvaluateWithGradient(const mat &mG, mat &grad)
	{
	    const vec PP = P % mvnMGF(X.cols(G), mG, s(G));

	    double res = accu(log1p(PP) - y % (X.cols(G) * mG)) +
		lambda * sqrt(accu(s(G) % s(G) + mG % mG)); 
	    
	    vec dPPmG = vec(mG.size(), arma::fill::zeros);

	    for (uword j = 0; j < mG.size(); ++j) {
		dPPmG(j) = accu( X.col(G(j)) % PP / (1 + PP) );
	    }

	    grad = 

	    return res;
	};

    private:
	const vec &y;
	const mat &X;
	const vec &mu;
	const vec &s;
	const double lambda;
	const uvec &G;
	const vec &P;
};


// [[Rcpp::export]]
mat test_minus(mat a, vec b) {
    return a.each_col() - b;
}



vec jen_update_mu()
{
    
}


vec jen_update_s()
{

} 


vec jen_update_g()
{

}



