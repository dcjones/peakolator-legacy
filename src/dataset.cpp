
#include "dataset.hpp"
#include "context.hpp"
#include "logger.h"
#include "miscmath.hpp"

#include <algorithm>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_math.h>

#include <nlopt.h>


using namespace std;

template <typename T> T sq( T x ) { return x*x; }

dataset::dataset( const char* reads_fn )
    : bias(NULL), bias_owner(false)
{
    log_printf( LOG_MSG, "loading reads from %s ...\n", reads_fn );
    log_indent();

    /* load SAM/BAM file */
    this->reads_fn = strdup(reads_fn);

    reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        failf( "Can't open bam file '%s'.", reads_fn );
    }

    reads_index = bam_index_load( reads_fn );
    if( reads_index == NULL ) {
        failf( "Can't open bam index '%s.bai'.", reads_fn );
    }

    bam_init_header_hash( reads_f->header );


    /* hash reads / make read counts */
    log_printf( LOG_MSG, "hashing ... \n" );
    log_indent();
    table_create( &T, n_targets() );
    bam1_t* read = bam_init1();
    size_t k = 0;

    while( samread( reads_f, read ) > 0 ) {
        k++;
        if( k % 1000000 == 0 ) {
            log_printf( LOG_MSG, "%zu reads\n", k );
        }
        table_inc( &T, read );
    }
    log_printf( LOG_MSG, "%zu reads\n", k );

    log_puts( LOG_MSG, "done.\n" );
    log_unindent();

    bam_destroy1( read );

    log_printf( LOG_MSG, "sorting ... " );
    read_counts_create( &counts, &T );
    log_printf( LOG_MSG, "done." );

    log_unindent();
    log_puts( LOG_MSG, "done.\n" );
}


void dataset::fit_sequence_bias( const char* ref_fn,
                                 size_t max_reads, pos L, pos R,
                                 double complexity_penalty,
                                 double offset_std )
{
    if( bias != NULL ) delete bias;
    bias = new sequencing_bias( ref_fn, &T, max_reads, L, R,
                                complexity_penalty, offset_std );
    bias_owner = true;
}


void dataset::load_sequence_bias( const char* ref_fn, const char* bias_fn )
{
    if( bias != NULL ) delete bias;
    bias = new sequencing_bias( ref_fn, bias_fn );
    bias_owner = true;
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
        failf( "Can't open bam file '%s'.", reads_fn );
    }

    pd->reads_index = bam_index_load( reads_fn );
    if( reads_index == NULL ) {
        failf( "Can't open bam index '%s.bai'.", reads_fn );
    }

    bam_init_header_hash( pd->reads_f->header );

    pd->reads_fn = strdup(reads_fn);

    table_copy( &pd->T, &T );
    read_counts_copy( &pd->counts, &counts );

    return pd;
}

dataset::~dataset()
{
    if( bias && bias_owner) delete bias;
    table_destroy( &T );
    read_counts_destroy( &counts );
    bam_index_destroy(reads_index);
    samclose(reads_f);
    free(reads_fn);
}


/* zero-inflated negative binomial log-likelihood function */
double zinb_ll_f( unsigned int n, const double* rp, double* grad, void* params )
{
    gsl_histogram* ksh = (gsl_histogram*)params;
    double r = rp[0];
    double p = rp[1];

    double ll = 0.0;

    size_t k;
    for( k = 1; k < ksh->n; k++ ) {
        if( ksh->bin[k] > 0 ) ll += ksh->bin[k] * lddnbinom( k, r, p );
    }

    log_printf( LOG_BLAB, "r = %0.4e, p = %0.4e, L = %0.8e\n",
                r, p, ll );

    return ll;
}


double nb_ll_f( unsigned int n, const double* rp, double* grad, void* params )
{
    gsl_histogram* ksh = (gsl_histogram*)params;
    double r = rp[0];
    double p = rp[1];

    double ll = 0.0;
    size_t k;
    for( k = 0; k < ksh->n; k++ ) {
        ll += ksh->bin[k] * ldnbinom( k, r, p );
    }

    log_printf( LOG_BLAB, "r = %0.4e, p = %0.4e, L = %0.4e\n",
                r, p, ll );

    return ll;
}


void dataset::fit_null_distr( interval_stack* is, double* r, double* p, double* a )
{
    log_puts( LOG_MSG, "fitting null distribution ... \n" );
    log_indent();

    const size_t max_k = 100000;
    uint64_t* ks = new uint64_t[max_k+1];

    memset( ks, 0, (max_k+1) * sizeof(uint64_t) );

    int32_t tid;
    interval_stack::iterator i;
    for( i = is->begin(); i != is->end(); i++ ) {
        tid = bam_get_tid( reads_f->header, i->seqname );
        read_count_occurances( &counts, tid, i->start, i->end, i->strand,
                               ks, max_k );
    }

    /* copy ks to a gsl histogram, for conveniance  */
    gsl_histogram* ksh = gsl_histogram_alloc( max_k+1 );
    gsl_histogram_set_ranges_uniform( ksh, 0.0, (double)(max_k+1) );



    /* for debugging purposes */
    size_t j;
    for( j = 0; j <= max_k; j++ ) {
        if( (j + 1) % 10 == 0 ) log_printf( LOG_BLAB, "\n" );
        log_printf( LOG_BLAB, "[%zu] %llu  ", j, ks[j] );
    }
    log_printf( LOG_BLAB, "\n" );



    for( j = 0; j <= max_k; j++ ) ksh->bin[j] = (double)ks[j];
    delete[] ks;


    const double TINY_VAL = 1e-10;

    const double lower[]     = { 1.0, TINY_VAL };
    const double upper[]     = { HUGE_VAL, 1.0 - TINY_VAL };
    const double step_size[] = { 0.5, 0.1 };
    
    nlopt_opt fmax = nlopt_create( NLOPT_LN_SBPLX, 2 );
    nlopt_set_max_objective( fmax, zinb_ll_f, (void*)ksh );
    nlopt_set_lower_bounds( fmax, lower );
    nlopt_set_upper_bounds( fmax, upper );
    nlopt_set_initial_step( fmax, step_size );
    nlopt_set_ftol_rel( fmax, 1e-12 );
    nlopt_set_maxeval( fmax, 1000 );

    double rp[2];
    /* initialize with the method of moments estimations */
    double ks_mean = gsl_histogram_mean( ksh );
    double ks_var  = sq( gsl_histogram_sigma( ksh ) );

    /* variance > mean, for a negative binomial distribution */
    ks_var = max( ks_var, ks_mean + 1e-6 );

    rp[0] = sq( ks_mean ) / (ks_var - ks_mean);
    rp[1] = rp[0] / (rp[0] + ks_mean);

    rp[0] = min( upper[0], max( lower[0], rp[0] ) );
    rp[1] = min( upper[1], max( lower[1], rp[1] ) );

    log_puts( LOG_MSG, "optimizing ... " );
    log_indent();

    double f_opt;
    nlopt_optimize( fmax, rp, &f_opt );


    *r = rp[0];
    *p = rp[1];

    *a = 0;
    for( j = 1; j <= max_k; j++ ) *a += ksh->bin[j];
    *a = (double)ksh->bin[0] / (*a + (double)ksh->bin[0]);

    log_unindent();
    log_printf( LOG_MSG, "done. (r = %0.4e, p = %0.4e, a = %0.4e\n", *r, *p, *a );

    gsl_histogram_free( ksh );
    log_unindent();
}







