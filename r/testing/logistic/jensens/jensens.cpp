#include "RcppEnsmallen.h"

// [[Rcpp::depends(RcppEnsmallen)]]

using namespace arma;

// [[Rcpp::export]]
vec mvnMGF(const mat &X, const vec &mu, const vec &sig) 
{
    return exp(X * mu + 0.5 * (X % X) * (sig % sig));
}

// [[Rcpp::export]]
inline vec compute_P_G(const mat &X, const vec &mu, const vec &s, const vec &g, 
	const uvec &G)
{
    return ( (1 - g(G(0))) + g(G(0)) * mvnMGF(X.cols(G), mu(G), s(G)) );
}

// [[Rcpp::export]]
vec compute_P(const mat &X, const vec &mu, const vec &s, const vec &g, 
	const uvec &groups)
{
    vec P = vec(X.n_rows, arma::fill::ones);
    const uvec ugroups = unique(groups);

    for (uword group : ugroups) {
	uvec G = find(groups == group);
	P %= compute_P_G(X, mu, s, g, G);
    }
    return P;
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

	    grad = dPPmG -
		X.cols(G).t() * y +
		lambda * mG * pow(dot(s(G), s(G)) + dot(mG, mG), -0.5);
	    
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
vec jen_update_mu(const vec &y, const mat &X, const vec &mu, const vec &s,
	const double lambda, const uvec &G, const vec &P)
{
    ens::L_BFGS opt;
    opt.MaxIterations() = 50;
    jen_update_mu_fn fn(y, X, mu, s, lambda, G, P);

    arma::vec mG = mu(G);
    opt.Optimize(fn, mG);

    return mG;
}


class jen_update_s_fn
{
    public:
	jen_update_s_fn(const vec &y, const mat &X, const vec &mu,
		const double lambda, const uvec &G, const vec &P) :
	    y(y), X(X), mu(mu), lambda(lambda), G(G), P(P)
	{};

	double EvaluateWithGradient(const mat &u, mat &grad)
	{
	    const vec sG = exp(u);

	    const vec PP = P % mvnMGF(X.cols(G), mu(G), sG);

	    double res = accu(log1p(PP)) -
		accu(log(sG)) +
		lambda * sqrt(accu(sG % sG + mu(G) % mu(G)));
	    
	    vec dPPsG = vec(sG.size(), arma::fill::zeros);

	    for (uword j = 0; j < sG.size(); ++j) {
		dPPsG(j) = accu( sG(j) * (X.col(G(j)) % X.col(G(j))) % 
			PP / (1 + PP) );
	    }

	    // df/duG = df/dsG * dsG/du
	    grad = (dPPsG -
		1.0 / sG +
		lambda * sG * pow(dot(sG, sG) + dot(mu(G), mu(G)), -0.5)) % sG;
	    
	    // Rcpp::Rcout << res;
	    return res;
	};

    private:
	const vec &y;
	const mat &X;
	const vec &mu;
	const double lambda;
	const uvec &G;
	const vec &P;
};


// [[Rcpp::export]]
vec jen_update_s(const vec &y, const mat &X, const vec &mu, const vec &s,
	const double lambda, const uvec &G, const vec &P)
{
    ens::L_BFGS opt;
    opt.MaxIterations() = 50;
    jen_update_s_fn fn(y, X, mu, lambda, G, P);

    arma::vec u = log(s(G));
    opt.Optimize(fn, u);

    return exp(u);
}


// [[Rcpp::export]]
double jen_update_g(const vec &y, const mat &X, const vec &mu, const vec &s,
	const double lambda, const double w, const uvec &G, const vec &P)
{
    const double mk = G.size();
    const double Ck = mk * log(2.0) + 0.5*(mk-1.0)*log(M_PI) + 
	lgamma(0.5*(mk + 1.0));

    const vec PP = P % mvnMGF(X.cols(G), mu(G), s(G));

    const double res =
	log(w / (1 - w)) + 
	0.5 * mk - 
	Ck +
	mk * log(lambda) +
	0.5 * accu(log(2.0 * M_PI * s(G) % s(G))) -
	lambda * sqrt(dot(s(G), s(G)) + dot(mu(G), mu(G))) +
	dot(y, (X.cols(G) * mu(G))) -
	accu(log1p(PP)) + 
	accu(log1p(P));

    return 1.0/(1.0 + exp(-res));
}

