
#include "dataset.hpp"
#include "context.hpp"
#include "logger.h"

#include <algorithm>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <nlopt.h>

#include <boost/math/distributions/negative_binomial.hpp>



using namespace std;

dataset::dataset(
        const char* ref_fn, const char* reads_fn,
        size_t bias_n, pos bias_L, pos bias_R )
{
    log_printf( LOG_MSG, "loading reads from %s ... ", reads_fn );

    if( ref_fn ) {
        bias = new sequencing_bias( ref_fn, reads_fn,
                                    bias_n, bias_L, bias_R );
    }
    else bias = NULL;

    reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        log_printf( LOG_ERROR, "Can't open bam file '%s'.", reads_fn );
        exit(1);
    }

    reads_index = bam_index_load( reads_fn );
    if( reads_index == NULL ) {
        log_printf( LOG_ERROR, "Can't open bam index '%s.bai'.", reads_fn );
        exit(1);
    }

    log_puts( LOG_MSG, "done.\n" );

    this->reads_fn = strdup(reads_fn);
}

dataset::dataset()
{
}

dataset* dataset::copy() const
{
    dataset* pd = new dataset();

    pd->bias = bias ? bias->copy() : NULL;

    pd->reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        log_printf( LOG_ERROR, "Can't open bam file '%s'.", reads_fn );
        exit(1);
    }

    pd->reads_index = bam_index_load( reads_fn );
    if( reads_index == NULL ) {
        log_printf( LOG_ERROR, "Can't open bam index '%s.bai'.", reads_fn );
        exit(1);
    }

    pd->reads_fn = strdup(reads_fn);

    return pd;
}

dataset::~dataset()
{
    bam_index_destroy(reads_index);
    samclose(reads_f);
    delete bias;
    free(reads_fn);
}


/* fit a negative binomial distribution to a number of training samples */
double nb_ll_f( unsigned int n, const double* rp, double* grad, void* params )
{
    gsl_matrix* CR = (gsl_matrix*)params;
    double r = rp[0];
    double p = rp[1];

    size_t i;
    double c_i;
    double r_i;
    double rr_i;

    double ll = 0.0;
    for( i = 0; i < CR->size1; i++ ) {

        c_i = gsl_matrix_get( CR, i, 0 );
        r_i = gsl_matrix_get( CR, i, 1 );
        rr_i = r*r_i;

        ll +=   rr_i * log(p)
           +    c_i  * log( 1.0 - p )
           +    gsl_sf_lngamma( rr_i + c_i )
           -    gsl_sf_lnfact( (unsigned int)c_i )
           -    gsl_sf_lngamma( rr_i );
    }

    return ll;
}



void dataset::fit_null_distr( interval_stack* is, double* r, double* p )
{
    log_puts( LOG_MSG, "fitting null distribution ... \n" );
    log_indent();


    /* STEP 1: Get counts and rates. */

    /* training data: two columns are:
     *      0: read counts, 1: rates */
    gsl_matrix* CR = gsl_matrix_alloc( is->size(), 2 );


    log_puts( LOG_MSG, "getting counts and rates ... " );
    context ctx;
    
    int u = 0;
    interval_stack::iterator i;
    for( i = is->begin(); i != is->end(); i++ ) {
        ctx.set( this, *i );

        gsl_matrix_set( CR, u, 0, max(1.0,(double)ctx.count( 0, ctx.length()-1 )) );
        gsl_matrix_set( CR, u, 1, (double)ctx.rate( 0, ctx.length()-1 ) );

        u++;
    }


    /* Get mean and variance of counts */
    gsl_vector_view C = gsl_matrix_column( CR, 0 );
    gsl_vector_view R = gsl_matrix_column( CR, 1 );

    double c_mean = gsl_stats_mean( C.vector.data, C.vector.stride, C.vector.size );
    double c_var  = gsl_stats_variance_m( C.vector.data, C.vector.stride,
                                          C.vector.size, c_mean );
    double r_mean = gsl_stats_mean( R.vector.data, R.vector.stride, R.vector.size );

    log_puts( LOG_MSG, "done.\n" );



    /* STEP 2: Initialize minimizer */

    const double lower[]     = { 1e-20, 1e-20 };
    const double upper[]     = { HUGE_VAL, 1.0 };
    const double step_size[] = { 1e-2, 1e-2 };


    nlopt_opt fmax = nlopt_create( NLOPT_LN_SBPLX, 2 );
    nlopt_set_max_objective( fmax, nb_ll_f, (void*)CR );
    nlopt_set_lower_bounds( fmax, lower );
    nlopt_set_upper_bounds( fmax, upper );
    nlopt_set_initial_step( fmax, step_size );
    nlopt_set_ftol_rel( fmax, 1e-6 );
    nlopt_set_maxeval( fmax, 5000 );

    double rp[2];

    /* initialize with essentially method of moments estimations */
    rp[0] = (c_mean*c_mean) / ((c_var - c_mean)*r_mean);
    rp[1] = rp[0] / (rp[0]*c_mean);


    log_puts( LOG_MSG, "fitting null distribution ... " );
    log_indent();

    double f_opt;
    nlopt_optimize( fmax, rp, &f_opt );

    log_unindent();
    log_puts( LOG_MSG, "done.\n" );


    *r = rp[0];
    *p = rp[1];

    gsl_matrix_free( CR );
    nlopt_destroy(fmax);

    log_unindent();
}

