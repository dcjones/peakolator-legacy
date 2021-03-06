

#include "emppval.hpp"
#include "common.hpp"
#include "miscmath.hpp"
#include "scanner.hpp"
#include "logger.h"

#include <cmath>
#include <cstdio>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>

#include <nlopt.h>


template <typename T> T sq( T x ) { return x*x; }


/* context for maximum liklihood fitting of the gev model */
struct gev_fit
{
    gsl_matrix* xs;
    gsl_vector* ls;
};



double gev_objf( unsigned k, const double* theta, double* grad, void* extra )
{
    gev_fit* gf = (gev_fit*)extra;

    double l;
    double fx = 0.0;

    size_t i, j;
    for( i = 0; i < gf->xs->size1; i++ ) {
        l = gsl_vector_get( gf->ls, i );
        for( j = 0; j < gf->xs->size2; j++ ) {
            fx += ldgev( gsl_matrix_get( gf->xs, i, j ),
                         theta[0] + theta[1] * l,
                         theta[2],
                         theta[3] );
        }
    }

    log_printf( LOG_BLAB, "gev_objf: fx = %0.4e\n", fx );
    log_printf( LOG_BLAB, "\ttheta = (%0.2e, %0.2e, %0.2e, %0.2e)\n",
                theta[0], theta[1], theta[2], theta[3] );

    return fx;
}





emppval::emppval( parameters* params  )
{
    log_puts( LOG_MSG, "training emperical p-value model ...\n" );
    log_indent();

    gev_fit gf;

    /* the number of sizes being generated */
    const size_t n = 24;

    /* the number of repititions for each size */
    const size_t m = 50;

    /* lengths vector */
    gf.ls = gsl_vector_alloc( n );

    /* results for each monte carlo round */
    gf.xs = gsl_matrix_alloc( n, m );

    /* dummy context and scanner, for monte's sake */
    context ctx;
    scanner* M;

    log_puts( LOG_MSG, "sampling ... \n" );
    size_t i,j;
    pos l;
    double x;
    for( i = 0; i < n; i++ ) {

        /* generate noise of exponentialy increasing lengths */
        l = (pos)pow( 1.5, (double)(i+10) );
        gsl_vector_set( gf.ls, i, log((double)l) );

        M = new scanner( params, &ctx );
        for( j = 0; j < m; j++ ) {
            ctx.set_noise( params->dist, l );
            x = M->least_likely_interval( 0, l-1, 0.0 ).score;
            gsl_matrix_set( gf.xs, i, j, x );
        }
        delete M;
        M = NULL;
    }
    log_puts( LOG_MSG, "done.\n" );



    /* numerical optimization */

    log_puts( LOG_MSG, "optimizing ...\n" );
    nlopt_opt fsolve = nlopt_create( NLOPT_LN_SBPLX, 4 );

    //double theta[6]; [> ( mu0, mu1, sigma0, sigma1, xi0, xi1 ) <]
    double theta[4]; /* ( mu0, mu1, sigma0, xi0 ) */
    double f_opt;


    /* choose a reasonable starting point, using method of moments estimation */
    const double xs_mean = gsl_stats_mean( gf.xs->data, 1, n*m );
    const double xs_var  = gsl_stats_variance_m( gf.xs->data, 1, n*m, xs_mean );

    theta[2] = sqrt( 6.0 * xs_var ) / M_PI;
    theta[0] = xs_mean - 0.58 * theta[2];
    theta[1] = 0.0;
    theta[3] = 0.0;




    nlopt_set_max_objective( fsolve, gev_objf, (void*)&gf );
    
    /* bounds */
    const double TINY_VAL = 1e-20;
    const double lower_bounds[] = { -HUGE_VAL, -HUGE_VAL,
                                     TINY_VAL, -HUGE_VAL };
    const double upper_bounds[] = { HUGE_VAL, HUGE_VAL,
                                    HUGE_VAL, HUGE_VAL };

    nlopt_set_lower_bounds( fsolve, lower_bounds );
    nlopt_set_upper_bounds( fsolve, upper_bounds );


    /* step size */
    const double step_size[] = { 1.0, 1.0, 1.0, 0.1 };
    nlopt_set_initial_step( fsolve, step_size );

    /* stopping criteria */
    nlopt_set_ftol_rel( fsolve, 1e-12 );
    nlopt_set_maxeval( fsolve, 1000 );

    /* optimize! */
    int err = nlopt_optimize( fsolve, theta, &f_opt );

    if( err < 0 ) {
        failf( "Unable to fit emperical p-value model (err = %d)\n", err );
    }


    c_mu[0]    = theta[0];
    c_mu[1]    = theta[1];
    c_sigma[0] = theta[2];
    c_sigma[1] = 0.0;
    c_xi[0]    = theta[3];
    c_xi[1]    = 0.0;


    nlopt_destroy( fsolve );
    gsl_matrix_free( gf.xs );
    gsl_vector_free( gf.ls );

    log_puts( LOG_MSG, "done.\n" );
    log_unindent();
    log_puts( LOG_MSG, "done.\n" );
}



emppval::emppval( const emppval& padj )
{
    c_mu[0]    = padj.c_mu[0];    c_mu[1]    = padj.c_mu[1];
    c_sigma[0] = padj.c_sigma[0]; c_sigma[1] = padj.c_sigma[1];
    c_xi[0]    = padj.c_xi[0];    c_xi[1]    = padj.c_xi[1];

    params = padj.params;
}




emppval::~emppval()
{
}


double emppval::operator() ( double score, pos len ) const
{
    /* find lerped parameters */
    double mu, sigma, xi;

    mu    = c_mu[0]    + log((double)len) * c_mu[1];
    sigma = c_sigma[0] + log((double)len) * c_sigma[1];
    xi    = c_xi[0]    + log((double)len) * c_xi[1];

    if( sigma <= 0.0 ) sigma = 1e-20;


    /* adjust */
    double pval = pgev( score, mu, sigma, xi );

    log_printf( LOG_BLAB, "pvalue( %0.4e, %d ) = %0.4e\n", score, len, pval );

    /* Never adjust a p-value to make it more significant. This is mainly to
     * avoid strange tail behavior (due to imprecision of the fit) producing 0
     * p-values. */
    if( log(pval) < score ) pval = exp(score);

    return pval;
}


