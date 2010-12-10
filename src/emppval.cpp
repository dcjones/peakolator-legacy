

#include "emppval.hpp"
#include "scanner.hpp"
#include "logger.h"

#include <gsl/gsl_statistics_double.h>


template <typename T> T sq( T x ) { return x*x; }



/*
 *                    The Generalized Extreme Value Distribution
 */

/* Distribution function */
template <typename T>
T gev_cdf( T q, const T& loc, const T& scale, const T& shape, bool lowertail = false )
{
    T p;
    q = (q-loc)/scale;
    if( shape == 0.0 ) p = exp(-exp(-q));
    else p = exp( -pow( max( 1.0 + shape * q, T(0.0) ), -1.0/shape ) );

    if( lowertail ) p = 1.0 - p;
    return p;
}



/* Log-liklihood */
template <typename T>
T gev_ll( const T& x, const T& mu, const T& sigma, const T& xi )
{
    T hx;
    if( sigma <= 0.0 ) return -HUGE_VAL;
    if( xi == 0.0 ) {
        hx = (mu - x) / sigma;
        return -log(sigma) + hx - exp(hx);
    }
    else {
        hx = 1.0 + xi * (x - mu) / sigma;
        if( hx <= 0  ) return -HUGE_VAL;
        return -log(sigma) - (1.0 + 1.0/xi) * log(hx) - pow( hx, -1.0/xi );
    }
}



/* Objective function used to maximize likelihood */
double gev_objf( unsigned int n, const double* params, double* grad, void* model )
{
    emppval* M = (emppval*)model;
    const double mu    = params[0];
    const double sigma = params[1];
    const double xi    = params[2];

    double x;
    double fx = 0.0;

    size_t i;
    for( i = 0; i < M->params->n_mc; i++ ) {
        x = M->qx_mc[i];
        fx += gev_ll( x, mu, sigma, xi );
    }

    log_printf( LOG_BLAB, "\tgev_objf: fx = %0.4e\n\tmu = %0.4e, sigma = %0.4e, xi = %0.4e\n",
                fx, mu, sigma, xi );

    return fx;
}





emppval::emppval( parameters* params  )
    : mu(NULL), sigma(NULL), xi(NULL)
{
    log_puts( LOG_MSG, "initializing p-value adjustment model ...\n" );
    log_indent();

    this->params  = params;
    this->n       = params->padj_n;
    this->spacing = params->padj_spacing;

    mu    = new mpfr_class[n];
    sigma = new mpfr_class[n];
    xi    = new mpfr_class[n];

    double mu_i, sigma_i, xi_i;

    qx_mc = new double[params->n_mc];

    /* dummy context */
    context ctx;
    scanner* M;

    int i;
    size_t j;
    for( i = 1; i < n; i++ ) {
        log_printf( LOG_MSG, "fitting length = %d ... ", i*spacing );

        /* generate examples */
        M = new scanner( params, &ctx );

        for( j = 0; j < params->n_mc; j++ ) {
            ctx.set_noise( params->dist, spacing*i );
            qx_mc[j] = mpfr_class(log( M->least_likely_interval( 0, i*spacing-1, 1.0 ).pval )).get_d();
        }
        delete M;
        M = NULL;

        /* fit */
        if( !fit_gev_to_emperical( &mu_i, &sigma_i, &xi_i ) ) {
            log_puts( LOG_ERROR, "GEVD MLE fitting failed! Please report/investigate.\n" );
            exit(1);
        }


        log_printf( LOG_BLAB,
                    "model fit for n = %ld.  loc = %0.8e, scale = %0.8e, shape = %0.8e\n",
                    i*spacing, mu_i, sigma_i, xi_i );


        mu[i]    = mu_i;
        sigma[i] = sigma_i;
        xi[i]    = xi_i;

        log_puts( LOG_MSG, "done.\n" );
    }


    /* XXX: for debugging */
#if 0
    /* report the last sample with fit values */
    FILE* evdlog = fopen( "evd.tab", "w" );
    for( j = 0; j < params->n_mc; j++ ) {
        fprintf( evdlog, "%.8e\n", qx_mc[j] );
    }
    fclose( evdlog );

    evdlog = fopen( "evd.params.tab", "w" );
    mpfr_fprintf( evdlog, "loc   <- %.8e\n", mu_i );
    mpfr_fprintf( evdlog, "scale <- %.8e\n", sigma_i );
    mpfr_fprintf( evdlog, "shape <- %.8e\n", xi_i );
    fclose( evdlog );
#endif
    /* XXX: end debugging output */



    delete[] qx_mc;

    log_unindent();
    log_puts( LOG_MSG, "done.\n" );
}



emppval::emppval( const emppval& padj )
{
    n       = padj.n;
    spacing = padj.spacing;

    mu    = new mpfr_class[n];
    sigma = new mpfr_class[n];
    xi    = new mpfr_class[n];

    int i;
    for( i = 0; i < n; i++ ) {
        mu[i]    = padj.mu[i];
        sigma[i] = padj.sigma[i];
        xi[i]    = padj.xi[i];
    }

    params = padj.params;
}




emppval::~emppval()
{
    delete[] mu;
    delete[] sigma;
    delete[] xi;
}


mpfr_class emppval::adjust( const mpfr_class& pval, pos len ) const
{
    /* find lerped parameters */
    mpfr_class mu_i, sigma_i, xi_i;
    if( len/spacing == 0 ) {
        mu_i    = mu[1];
        sigma_i = sigma[1];
        xi_i    = xi[1];
    }
    else if( len/spacing >= n-1 ) {
        mu_i    = mu[n-1];
        sigma_i = sigma[n-1];
        xi_i    = xi[n-1];
    }
    else {
        int i = len / spacing;
        double q = (double)len / (double)spacing;
        q = q - floor(q);

        mu_i    = mu[i]*(1-q)    + mu[i+1]*q;
        sigma_i = sigma[i]*(1-q) + sigma[i+1]*q;
        xi_i    = xi[i]*(1-q) + xi[i+1]*q;
    }

    /* adjust */
    mpfr_class padj = gev_cdf<mpfr_class>( log(pval), mu_i, sigma_i, xi_i );

    char *pval_str, *padj_str;
    mpfr_asprintf( &pval_str, "%.4Re", pval.get_mpfr_t() );
    mpfr_asprintf( &padj_str, "%.4Re", padj.get_mpfr_t() );

    log_printf( LOG_BLAB, "pvalue adjusted (len = %d): %s  -->  %s", len, pval_str, padj_str );

    free(pval_str);
    free(padj_str);

    /* Never adjust a p-value to make it more significant. This is mainly to
     * avoid strange tail behavior (due to imprecision of the fit) producing 0
     * p-values. */
    if( padj < pval ) padj = pval;

    return padj;
}


/* Maximimum likelihood fitting to the generalized extreme value distribution.
 * I make the simplifying assumption that xi != 0. So I fit with the constraint
 * xi > 0, then with xi < 0, and take the more likly parameters. */
bool emppval::fit_gev_to_emperical( double* mu, double* sigma, double* xi )
{
    nlopt_opt fsolve;

    double theta[3];
    double f_opt;

    double qx_mc_mean = gsl_stats_mean( qx_mc, 1, params->n_mc );
    double qx_mc_var  = gsl_stats_variance_m( qx_mc, 1, params->n_mc, qx_mc_mean );

    /* choose a reasonabl starting point */
    theta[1] = sqrt( 6.0 * qx_mc_var ) / M_PI;
    theta[0] = qx_mc_mean - 0.58 * theta[1];
    theta[2] = 0.0;


    const double TINY_VAL = 1e-20;
                                      /* mu       sigma      xi */
    const double lower_bounds[]  = { -HUGE_VAL, TINY_VAL,  -HUGE_VAL };

    const double step_size[]     = { 1.0, 1.0, 0.1 };


    /* solve 'low' case */

    fsolve = nlopt_create( NLOPT_LN_SBPLX, 3 );
    nlopt_set_max_objective( fsolve, gev_objf, (void*)this );

    /* bounds */
    nlopt_set_lower_bounds( fsolve, lower_bounds ); 

    /* step size */
    nlopt_set_initial_step( fsolve, step_size );

    /* stopping criteria */
    nlopt_set_ftol_rel( fsolve, 1e-8 );
    nlopt_set_maxeval(  fsolve, 1000 );

    /* optimize! */
    if( nlopt_optimize( fsolve, theta, &f_opt ) < 0 ) return false;

    nlopt_destroy(fsolve);

    *mu    = theta[0];
    *sigma = theta[1];
    *xi    = theta[2];

    return true;
}

