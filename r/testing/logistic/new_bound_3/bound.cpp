#include "RcppEnsmallen.h"
#include <bitset>

// [[Rcpp::depends(RcppEnsmallen)]]

double ell(arma::mat Xm, arma::mat Xs, arma::vec g, double thresh, int l);
arma::vec dell_dm(arma::mat X, arma::mat Xm, arma::mat Xs, arma::vec g,
	arma::uvec G, double thresh, int l);
arma::vec dell_ds(arma::mat X, arma::mat Xm, arma::mat Xs, arma::vec s,
	arma::vec g, arma::uvec G, double thresh, int l);
double tll(arma::vec mu, arma::vec sig, int l);
arma::vec dt_dm(arma::mat X, arma::vec mu, arma::vec sig, arma::uvec G, int l);
arma::vec dt_ds(arma::mat X, arma::vec s, arma::vec mu, arma::vec sig,
	arma::uvec G, int l);


// Rcpp::List fit(arma::vec y, arma::mat X, arma::uvec groups, const double lambda, 
//     const double a0, const double b0, arma::vec mu, arma::vec s, arma::vec g,
//     double thresh, unsigned int terms, unsigned int niter, double tol, bool verbose)
// {
//     const arma::uword n = X.n_rows;
//     const arma::uword p = X.n_cols;
//     const double w = a0 / (a0 + b0);
    
//     // compute commonly used expressions
//     const mat xtx = X.t() * X;
//     const double yty = dot(y, y);
//     const vec yx = (y.t() * X).t();
    
//     // init
//     const uvec ugroups = arma::unique(groups);

//     // if not constrained we are using a full covariance for S
//     std::vector<mat> Ss;
//     if (!diag_cov) {
// 	for (uword group : ugroups) {
// 	    uvec G = find(groups == group);	
// 	    Ss.push_back(arma::diagmat(s(G)));
// 	}
//     }

//     vec mu_old, s_old, g_old;
//     double tau_a = tau_a0, tau_b = tau_b0, e_tau = tau_a0 / tau_b0;

//     uword num_iter = niter;
//     bool converged = false;
//     std::vector<double> elbo_values;

//     for (unsigned int iter = 1; iter <= niter; ++iter)
//     {
// 	mu_old = mu; s_old = s; g_old = g;

// 	// update expected value of tau^2
// 	e_tau = tau_a / tau_b;

// 	// update mu, sigma, gamma
// 	for (arma::uword group : ugroups)
// 	{
// 	    uvec G  = arma::find(groups == group);
// 	    uvec Gc = arma::find(groups != group);
	    
// 	    uword gi = arma::find(ugroups == group).eval().at(0);
// 	    mat &S = Ss.at(gi);

// 	    mu(G) = update_mu(G, Gc, xtx, yx, mu, sqrt(diagvec(S)), g, e_tau, lambda);
// 	    s(G)  = update_S(G, xtx, mu, S, s(G), e_tau, lambda);
// 	    double tg = update_g(G, Gc, xtx, yx, mu, S, g, e_tau, lambda, w);
// 	    for (uword j : G) g(j) = tg;
// 	}
	
// 	// check for break, print iter
// 	Rcpp::checkUserInterrupt();
// 	if (verbose) Rcpp::Rcout << iter;
	
// 	// compute the ELBO if option enabled
// 	// if (track_elbo && (iter % track_elbo_every == 0)) {
// 	    // double e = diag_cov ?
// 		// elbo_linear_c(yty, yx, xtx, groups, n, p, mu, s, g, tau_a, tau_b,
// 		    // lambda, a0, b0, tau_a0, tau_b0, track_elbo_mcn, false) :
// 		// elbo_linear_u(yty, yx, xtx, groups, n, p, mu, Ss, g, tau_a, tau_b,
// 		    // lambda, a0, b0, tau_a0, tau_b0, track_elbo_mcn, false);

// 	    // elbo_values.push_back(e);
// 	}

// 	// check convergence
// 	if (sum(abs(mu_old - mu)) < tol &&
// 	    sum(abs(s_old - s))   < tol &&
// 	    sum(abs(g_old - g))   < tol) 
// 	{
// 	    if (verbose)
// 		Rcpp::Rcout << "\nConverged in " << iter << " iterations\n";

// 	    num_iter = iter;
// 	    converged = true;
// 	    break;
// 	}
//     }
    
//     // compute elbo for final eval
//     // if (track_elbo) {
// 	// double e = diag_cov ?
// 	    // elbo_linear_c(yty, yx, xtx, groups, n, p, mu, s, g, tau_a, tau_b,
// 		// lambda, a0, b0, tau_a0, tau_b0, track_elbo_mcn, false) :
// 	    // elbo_linear_u(yty, yx, xtx, groups, n, p, mu, Ss, g, tau_a, tau_b,
// 		// lambda, a0, b0, tau_a0, tau_b0, track_elbo_mcn, false);
// 	// elbo_values.push_back(e);
//     // }

//     return Rcpp::List::create(
// 	Rcpp::Named("mu") = mu,
// 	Rcpp::Named("sigma") = s,
// 	Rcpp::Named("S") = Ss,
// 	Rcpp::Named("gamma") = g,
// 	Rcpp::Named("tau_a") = tau_a,
// 	Rcpp::Named("tau_b") = tau_b,
// 	Rcpp::Named("converged") = converged,
// 	Rcpp::Named("iterations") = num_iter,
// 	Rcpp::Named("elbo") = elbo_values
//     );
// }


class update_m_fn
{
    public:
	update_m_fn(arma::vec y, arma::mat X, arma::vec m,
		arma::vec s, arma::vec g, double lambda, arma::uword group,
		arma::uvec G, arma::mat Xm, arma::mat Xs,
		double thresh, int l) :
	    y(y), X(X), m(m), s(s), g(g), lambda(lambda),
	    group(group), G(G), Xm(Xm), Xs(Xs), thresh(thresh), l(l)
	    { }

	double EvaluateWithGradient(const arma::mat &mG, arma::mat &grad) 
	{ 
	    g(group) = 1;
	    const arma::vec xm = X.cols(G) * mG;
	    Xm.col(group) = xm;

	    const double res = ell(Xm, Xs, g, thresh, l) -
		dot(y, xm) +
		lambda * sqrt(sum(s(G) % s(G) + mG % mG));

	    grad = dell_dm(X, Xm, Xs, g, G, thresh, l) -
		X.cols(G).t() * y +
		lambda * mG * pow(dot(s(G), s(G)) + dot(mG, mG), -0.5);

	    return res;
	}

    private:
	const arma::vec y;
	const arma::mat X;
	const arma::vec m;
	const arma::vec s;
	arma::vec g;
	const double lambda;
	const arma::uword group;
	const arma::uvec G;
	arma::mat Xm;
	const arma::mat Xs;
	const double thresh;
	const int l;
};


// [[Rcpp::export]]
arma::vec opt_m_cpp(arma::vec y, arma::mat X, arma::vec m,
	arma::vec s, arma::vec g, double lambda, arma::uword group,
	arma::uvec G, arma::mat Xm, arma::mat Xs,
	double thresh, int l)  
{
    ens::L_BFGS opt;
    opt.MaxIterations() = 50;
    update_m_fn fn(y, X, m, s, g, lambda, group, G, Xm, Xs, thresh, l);

    arma::vec mG = m(G);
    opt.Optimize(fn, mG);

    return mG;
}


// [[Rcpp::export]]
double fns(arma::vec sG, arma::mat Xm, arma::mat Xs, arma::vec m,
	arma::vec g, double thresh, int l, double lambda, arma::uvec G) 
{
    const double res = ell(Xm, Xs, g, thresh, l) -
	accu(log(sG)) +
	lambda * sqrt(dot(sG, sG) + dot(m(G), m(G)));

    return res;
}



class update_s_fn
{
    public:
	update_s_fn(arma::vec y, arma::mat X, arma::vec m,
		arma::vec s, arma::vec g, double lambda, arma::uword group,
		arma::uvec G, arma::mat Xm, arma::mat Xs,
		double thresh, int l) :
	    y(y), X(X), m(m), s(s), g(g), lambda(lambda),
	    group(group), G(G), Xm(Xm), Xs(Xs), thresh(thresh), l(l)
	    { }

	double EvaluateWithGradient(const arma::mat &u, arma::mat &grad) 
	{ 
	    arma::vec sG = exp(u);

	    g(group) = 1;
	    const arma::vec xs = (X.cols(G) % X.cols(G)) * (sG % sG);
	    Xs.col(group) = xs;

	    const double res = ell(Xm, Xs, g, thresh, l) -
		accu(log(sG)) +
		lambda * sqrt(dot(sG, sG) + dot(m(G), m(G)));

	    grad = (dell_ds(X, Xm, Xs, s, g, G, thresh, l) -
		1.0 / sG +
		lambda * sG * pow(dot(sG, sG) + dot(m(G), m(G)), -0.5)) % sG;

	    return res;
	}

    private:
	const arma::vec y;
	const arma::mat X;
	const arma::vec m;
	const arma::vec s;
	arma::vec g;
	const double lambda;
	const arma::uword group;
	const arma::uvec G;
	const arma::mat Xm;
	arma::mat Xs;
	const double thresh;
	const int l;
};


// [[Rcpp::export]]
arma::vec opt_s_cpp(arma::vec y, arma::mat X, arma::vec m,
	arma::vec s, arma::vec g, double lambda, arma::uword group,
	arma::uvec G, arma::mat Xm, arma::mat Xs,
	double thresh, int l)  
{
    ens::L_BFGS opt;
    opt.MaxIterations() = 50;
    update_s_fn fn(y, X, m, s, g, lambda, group, G, Xm, Xs, thresh, l);

    arma::vec u = log(s(G));
    opt.Optimize(fn, u);

    return exp(u);
}


// [[Rcpp::export]]
double opt_g_cpp(arma::vec y, arma::mat X, arma::vec m,
	arma::vec s, arma::vec g, double lambda, arma::uword group,
	arma::uvec G, arma::mat Xm, arma::mat Xs,
	double thresh, int l, const double w)  
{
    const double mk = G.size();
    const double Ck = mk * log(2.0) + 0.5*(mk-1.0)*log(M_PI) + 
	lgamma(0.5*(mk + 1.0));

    g(group) = 1;
    const double S1 = ell(Xm, Xs, g, thresh, l);

    g(group) = 0;
    const double S0 = ell(Xm, Xs, g, thresh, l);

    const double res =
	log(w / (1- w)) + 
	0.5 * mk - 
	Ck +
	mk * log(lambda) +
	0.5 * accu(log(2.0 * M_PI * s(G) % s(G))) -
	lambda * sqrt(dot(s(G), s(G)) + dot(m(G), m(G))) -
	S1 + S0 + dot(y, (X.cols(G) * m(G)));
    
    return 1.0 / (1.0 + exp(-res));
}


// [[Rcpp::export]]
double ell(arma::mat Xm, arma::mat Xs, arma::vec g, double thresh, int l)
{
    const arma::uvec mid = find(g >= thresh && g <= (1.0 - thresh));
    const arma::uvec big = find(g > (1.0 - thresh));
    const int msize = mid.size();

    arma::vec mu = sum(Xm.cols(big), 1);
    arma::vec sig = sum(Xs.cols(big), 1);

    double res = 0.0;
    if (msize == 0) 
    {
	res = tll(mu, sqrt(sig), l);
    } 
    else 
    {
	double tot = 0.0;
	for (int i = 0; i < pow(2, msize); ++i) 
	{
	    auto b = std::bitset<12>(i);
	    double prod_g = 1.0;

	    arma::vec mu_new = mu;
	    arma::vec sig_new = sig;

	    for (int j = 0; j < msize; ++j) 
	    {
		if ((i >> j) & 1) 
		{
		    mu_new += Xm.col(mid(j));
		    sig_new += Xs.col(mid(j));
		    prod_g *= g(mid(j));
		} 
		else 
		{
		    prod_g *= (1 - g(mid(j)));
		}
	    }
	    
	    tot += prod_g * tll(mu_new, sqrt(sig_new), l);
	}

	res = tot;
    }

    return  res;
}


// [[Rcpp::export]]
arma::vec dell_dm(arma::mat X, arma::mat Xm, arma::mat Xs, arma::vec g,
	arma::uvec G, double thresh, int l)
{
    const int mk = G.n_rows;
    arma::vec res = arma::vec(mk, arma::fill::zeros);

    const arma::uvec mid = find(g >= thresh && g <= (1.0 - thresh));
    const arma::uvec big = find(g > (1.0 - thresh));
    const int msize = mid.size();

    arma::vec mu = sum(Xm.cols(big), 1);
    arma::vec sig = sum(Xs.cols(big), 1);

    if (msize == 0) 
    {
	res = dt_dm(X, mu, sqrt(sig), G, l);
    } 
    else 
    {
	for (int i = 0; i < pow(2, msize); ++i) 
	{
	    auto b = std::bitset<12>(i);
	    double prod_g = 1.0;

	    arma::vec mu_new = mu;
	    arma::vec sig_new = sig;

	    for (int j = 0; j < msize; ++j) 
	    {
		if ((i >> j) & 1) 
		{
		    mu_new += Xm.col(mid(j));
		    sig_new += Xs.col(mid(j));
		    prod_g *= g(mid(j));
		} 
		else 
		{
		    prod_g *= (1 - g(mid(j)));
		}
	    }
	    res += prod_g * dt_dm(X, mu_new, sqrt(sig_new), G, l);
	}
    }

    return  res;
}


// [[Rcpp::export]]
arma::vec dell_ds(arma::mat X, arma::mat Xm, arma::mat Xs, arma::vec s,
	arma::vec g, arma::uvec G, double thresh, int l) 
{
    const int mk = G.n_rows;
    arma::vec res = arma::vec(mk, arma::fill::zeros);

    const arma::uvec mid = find(g >= thresh && g <= (1.0 - thresh));
    const arma::uvec big = find(g > (1.0 - thresh));
    const int msize = mid.size();

    arma::vec mu = sum(Xm.cols(big), 1);
    arma::vec sig = sum(Xs.cols(big), 1);

    if (msize == 0) 
    {
	res = dt_ds(X, s, mu, sqrt(sig), G, l);
    } 
    else 
    {
	double tot = 0.0;
	for (int i = 0; i < pow(2, msize); ++i) 
	{
	    auto b = std::bitset<12>(i);
	    double prod_g = 1.0;

	    arma::vec mu_new = mu;
	    arma::vec sig_new = sig;

	    for (int j = 0; j < msize; ++j) 
	    {
		if ((i >> j) & 1) 
		{
		    mu_new += Xm.col(mid(j));
		    sig_new += Xs.col(mid(j));
		    prod_g *= g(mid(j));
		} 
		else 
		{
		    prod_g *= (1 - g(mid(j)));
		}
	    }
	    res += prod_g * dt_ds(X, s, mu_new, sqrt(sig_new), G, l);
	}
    }

    return  res;
}


// [[Rcpp::export]]
double tll(arma::vec mu, arma::vec sig, int l)
{
    const int n = mu.n_rows;
    double res = 0.0;

    for (int i = 0; i < n; ++i) 
    {
	double a = sig(i) / sqrt(2.0 * M_PI) * 
	    exp(- 0.5 * mu(i)*mu(i) / (sig(i)*sig(i))) +
	    mu(i) * R::pnorm(mu(i) / sig(i), 0, 1, 1, 0);

	double b = 0.0;
	for(int j = 1; j <= (2 * l - 1); ++j) {
	    b += pow((-1.0), (j-1)) / j * (
		exp(
		    mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(-mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
		) +
		exp(
		    -mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
		)
	    );
	}

	res += a + b;
    }
    
    return res;
}


// [[Rcpp::export]]
arma::vec dt_dm(arma::mat X, arma::vec mu, arma::vec sig, arma::uvec G, int l) 
{
    // gradient wrt. mu
    const int n = mu.n_rows;
    const int mk = G.n_rows;
    arma::vec res = arma::vec(mk, arma::fill::zeros);

    for (int i = 0; i < n; ++i) 
    {
	double dt_dmu = 0.0;

	dt_dmu += sig(i) / sqrt(2.0 * M_PI) * 
	    - mu(i) / (sig(i)*sig(i)) * exp(- 0.5 * mu(i)*mu(i)/(sig(i)*sig(i))) +
	    R::pnorm(mu(i) / sig(i), 0, 1, 1, 0) +
	    mu(i)/sig(i) * R::dnorm4(mu(i) / sig(i), 0, 1, 0);

	for(int j = 1; j <= (2 * l - 1); ++j) 
	{
	    dt_dmu += pow((-1.0), (j-1)) / j * (
		j * exp(
		    mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(-mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
		) + 
		- j * exp(
		    -mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
		) +
		- 1.0/sig(i) * exp(
		    mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::dnorm(-mu(i)/sig(i) - j*sig(i), 0, 1, 1)
		) +
		1.0/sig(i) * exp(
		    -mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::dnorm(mu(i)/sig(i) - j*sig(i), 0, 1, 1)
		)
	    );
	}
	
	for (int k = 0; k < mk; ++k) {
	    res(k) += dt_dmu * X(i, G(k));
	}
    }

    return res;
}


// [[Rcpp::export]]
arma::vec dt_ds(arma::mat X, arma::vec s, arma::vec mu, arma::vec sig, 
	arma::uvec G, int l) 
{
    // gradient wrt. s
    const int n = mu.n_rows;
    const int mk = G.n_rows;
    arma::vec res = arma::vec(mk, arma::fill::zeros);

    for (int i = 0; i < n; ++i) 
    {
	double dt_dsig = 0.0;

	dt_dsig = 1.0 / sqrt(2.0 * M_PI) * 
	    (1.0 + mu(i)*mu(i)/(sig(i)*sig(i))) *
	    exp(- 0.5 * mu(i)*mu(i) / (sig(i)*sig(i))) -
	    mu(i)*mu(i)/(sig(i)*sig(i)) * R::dnorm(mu(i) / sig(i), 0, 1, 0);

	for(int j = 1; j <= (2 * l - 1); ++j) 
	{
	    dt_dsig += pow((-1.0), (j-1)) / j * (
		j*j*sig(i) * exp(
		    mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(-mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
		) + 
		j*j*sig(i) * exp(
		    -mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::pnorm(mu(i)/sig(i) - j*sig(i), 0, 1, 1, 1)
		) +
		(mu(i)/(sig(i)*sig(i)) - j) * exp(
		    mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::dnorm(-mu(i)/sig(i) - j*sig(i), 0, 1, 1)
		) +
		(-mu(i)/(sig(i)*sig(i)) - j) * exp(
		    -mu(i)*j + 0.5*j*j*sig(i)*sig(i) + 
		    R::dnorm(mu(i)/sig(i) - j*sig(i), 0, 1, 1)
		)
	    );
	}
	
	for (int k = 0; k < mk; ++k) {
	    res(k) += dt_dsig / sig(i) * X(i, G(k)) * X(i, G(k)) * s(G(k));
	}
    }

    return res;
}



