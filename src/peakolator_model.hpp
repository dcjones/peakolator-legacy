/*
 *                                   _         _       _             
 *                  _ __   ___  __ _| | _____ | | __ _| |_ ___  _ __ 
 *                 | '_ \ / _ \/ _` | |/ / _ \| |/ _` | __/ _ \| '__|
 *                 | |_) |  __/ (_| |   < (_) | | (_| | || (_) | |   
 *                 | .__/ \___|\__,_|_|\_\___/|_|\__,_|\__\___/|_|   
 *                 |_|     
 *
 *                                      ~~~
 *
 *                          Finding enriched regions in 
 *                          high-throughput sequencing data.
 *
 *                                      ~~~
 *
 *                                  April, 2010
 *
 *                                 Daniel Jones
 *                           <dcjones@cs.washington.edu>
 */


#ifndef PEAKOLATOR_PEAKOLATOR_MODEL
#define PEAKOLATOR_PEAKOLATOR_MODEL

#include "common.hpp"
#include "peakolator_context.hpp"
#include "peakolator_interval.hpp"
#include "peakolator_dist.hpp"


#include <cstdlib>

#include "samtools/sam.h"



class peakolator_pval;


struct peakolator_parameters
{
    peakolator_parameters();
    peakolator_parameters( const peakolator_parameters& );
    ~peakolator_parameters();

    peakolator_parameters* copy() const;

    void rebuild_lookup( int m = 0, int n = 0 );
    void build_padj();

    /* adjusted p-value at which we consider a sub-interval significant */
    double alpha;

    /* negative binomial parameters */
    double r,p;

    /* Hard limits on the duration of segments. */
    /* IMPORTANT: for technical reasons, d_min > 2, must be true.
     * There is no upper limit for d_max, but d_min <= d_max must be true. */
    pos d_min; 
    pos d_max;

    /* number of samples to take to compute the emperical p-value */
    size_t n_mc;

    /* distribution lookup */
    peakolator_dist dist;

    /* pvalue model */
    pos padj_spacing;
    int padj_n;
    peakolator_pval* padj;

    private: void init_fsolve();
};



class peakolator_model
{
    public:
        peakolator_model( 
                peakolator_parameters* params,
                peakolator_context*    context );

        interval_stack* run();
        subinterval least_likely_interval   ( pos l, pos r, double alpha );

        ~peakolator_model();

        /* raw p-value for subinterval of duration d with read count x */
        mpfr_class QX( double r, rcount x );

        peakolator_context*    context;
        peakolator_parameters* params;

    private:
        void        trim_candidate( subinterval& I );

        peakolator_pval*       padj;
};




#endif

