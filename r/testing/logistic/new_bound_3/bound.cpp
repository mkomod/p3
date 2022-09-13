#include "RcppArmadillo.h"
#include <bitset>

// [[Rcpp::depends(RcppArmadillo)]]

double tll(arma::vec mu, arma::vec sig, int l);
arma::vec dt_dm(arma::mat X, arma::vec mu, arma::vec sig, arma::uvec G, int l);
arma::vec dt_ds(arma::mat X, arma::vec s, arma::vec mu, arma::vec sig,
	arma::uvec G, int l);


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

	    for (int j = 0; j < msize; ++j) 
	    {
		if ((i >> j) & 1) 
		{
		    mu += Xm.col(mid(j));
		    sig += Xs.col(mid(j));
		    prod_g *= g(mid(j));
		} 
		else 
		{
		    prod_g *= (1 - g(mid(j)));
		}
	    }
	    
	    tot += prod_g * tll(mu, sqrt(sig), l);
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

	    for (int j = 0; j < msize; ++j) 
	    {
		if ((i >> j) & 1) 
		{
		    mu += Xm.col(mid(j));
		    sig += Xs.col(mid(j));
		    prod_g *= g(mid(j));
		} 
		else 
		{
		    prod_g *= (1 - g(mid(j)));
		}
	    }
	    res += prod_g * dt_dm(X, mu, sqrt(sig), G, l);
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

	    for (int j = 0; j < msize; ++j) 
	    {
		if ((i >> j) & 1) 
		{
		    mu += Xm.col(mid(j));
		    sig += Xs.col(mid(j));
		    prod_g *= g(mid(j));
		} 
		else 
		{
		    prod_g *= (1 - g(mid(j)));
		}
	    }
	    res += prod_g * dt_ds(X, s, mu, sqrt(sig), G, l);
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
	
	for (arma::uword k : G) {
	    res(k) += dt_dmu * X(i, k);
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
	
	for (arma::uword k : G) {
	    res(k) += dt_dsig / sig(i) * X(i, k) * X(i, k) * s(k);
	}
    }

    return res;
}



