#include "RcppEnsmallen.h"

// [[Rcpp::depends(RcppEnsmallen)]]

using namespace arma;

class jaak_update_mu_fn
{
    public:
	jaak_update_mu_fn(const vec &y, const mat &X, const mat &XAX,
		const vec &mu, const vec &s, const vec &g, const double lambda,
		const uvec &G, const uvec &Gc) :
	    y(y), X(X), XAX(XAX), mu(mu), s(s), g(g), lambda(lambda), 
	    G(G), Gc(Gc)
	{};

	double EvaluateWithGradient(const mat &mG, mat &grad)
	{
	    const double res = 0.5 * dot(mG, XAX(G, G) * mG) +
		dot(mG, XAX(G, Gc) * (g(Gc) % mu(Gc))) +
		dot(0.5 - y, X.cols(G) * mG) +
		lambda * sqrt(accu(s(G) % s(G) + mG % mG)); 

	    // 0.5 * t(m_G) %*% XAX[G, G] %*% m_G +
	    // t(m_G) %*% XAX[G, Gc] %*% (g[Gc] * m[Gc]) +
	    // sum((0.5 - y) * X[ , G] %*% m_G) +
	    // lambda * (sum(s[G]^2 + m_G^2))^(1/2)

	    grad = XAX(G, G) * mG +
		XAX(G, Gc) * (g(Gc) % mu(Gc)) +
		X.cols(G).t() * (0.5 - y) +
		lambda * mG * pow(dot(s(G), s(G)) + dot(mG, mG), -0.5);
	    
	    return res;
	};

    private:
	const vec &y;
	const mat &X;
	const mat &XAX;
	const vec &mu;
	const vec &s;
	const vec &g;
	const double lambda;
	const uvec &G;
	const uvec &Gc;
};


// [[Rcpp::export]]
vec jaak_update_mu(const vec &y, const mat &X, const mat &XAX,
	const vec &mu, const vec &s, const vec &g, const double lambda,
	const uvec &G, const uvec &Gc)
{
    ens::L_BFGS opt;
    opt.MaxIterations() = 50;
    jaak_update_mu_fn fn(y, X, XAX, mu, s, g, lambda, G, Gc);

    vec mG = mu(G);
    opt.Optimize(fn, mG);

    return mG;
}


class jaak_update_s_fn
{
    public:
	jaak_update_s_fn(const vec &y, const mat &XAX, const vec &mu, 
		const double lambda, const uvec &G) :
	    y(y), XAX(XAX), mu(mu), lambda(lambda), G(G)
	{};

	double EvaluateWithGradient(const mat &u, mat &grad)
	{
	    const vec sG = exp(u);

	    const double res = 0.5 * accu(diagvec(XAX(G, G)) % sG % sG) -
		accu(log(sG)) +
		lambda * sqrt(accu(sG % sG + mu(G) % mu(G))); 

	    // 0.5 * sum(diag(XAX[G, G]) * s_G^2) -
	    // sum(log(s_G)) +
	    // lambda * (sum(s_G^2 + m[G]^2))^(1/2)

	    grad = (
		diagvec(XAX(G, G)) % sG -
		1 / sG +
		lambda * sG * pow(dot(sG, sG) + dot(mu(G), mu(G)), -0.5)
	    ) % sG;

	    return res;
	};

    private:
	const vec &y;
	const mat &XAX;
	const vec &mu;
	const double lambda;
	const uvec &G;
};


// [[Rcpp::export]]
vec jaak_update_s(const vec &y, const mat &XAX, const vec &mu, 
	const vec &s, const double lambda, const uvec &G)
{
    ens::L_BFGS opt;
    opt.MaxIterations() = 50;
    jaak_update_s_fn fn(y, XAX, mu,lambda, G);

    vec u = log(s(G));
    opt.Optimize(fn, u);

    return exp(u);
}


// [[Rcpp::export]]
double jaak_update_g(const vec &y, const mat &X, const mat &XAX,
	const vec &mu, const vec &s, const vec &g, const double lambda,
	const double w, const uvec &G, const uvec &Gc)
{
    const double mk = G.size();
    const double Ck = mk * log(2.0) + 0.5*(mk-1.0)*log(M_PI) + 
	lgamma(0.5*(mk + 1.0));

    const double res =
	log(w / (1 - w)) + 
	0.5 * mk - 
	Ck +
	mk * log(lambda) +
	0.5 * accu(log(2.0 * M_PI * s(G) % s(G))) -
	lambda * sqrt(dot(s(G), s(G)) + dot(mu(G), mu(G))) +
	dot((y - 0.5), X.cols(G) * mu(G)) -
	0.5 * dot(mu(G), XAX(G, G) * mu(G)) -
	0.5 * accu(diagvec(XAX(G, G)) % s(G) % s(G)) -
	dot(mu(G), XAX(G, Gc) * (g(Gc) % mu(Gc)));

    return 1.0 / (1.0 + exp(-res));
}


// [[Rcpp::export]]
vec jaak_update_l(const mat &X, const vec &mu, const vec &s, const vec &g) 
{
    return sqrt(pow(X * (g % mu), 2) + (X % X) * (g % s % s));
}

vec sigmoid(const vec &x) {
    return 1/(1 + exp(-x));
}

// [[Rcpp::export]]
vec a(const vec &x)
{
    return (sigmoid(x) - 0.5) / x;
}


