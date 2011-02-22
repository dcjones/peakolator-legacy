
#include "scanner.hpp"
#include "emppval.hpp"
#include "common.hpp"
#include "logger.h"

#include <algorithm>
#include <set>

#include <cstdlib>
#include <cmath>
#include <ctime>

#include <gsl/gsl_math.h>
#include <nlopt.h>

#include "samtools/sam.h"


using namespace std;


/* construction */
scanner::scanner( 
                parameters* params,
                context*    ctx )
{
    this->ctx = ctx;
    this->params  = params;

    if( !this->params->dist.ready() ) {
        params->dist.build( params->r, params->p, params->a );
    }

}


scanner::~scanner()
{
}



double scanner::QX( rcount x, unsigned int z, unsigned int d )
{
    return params->dist.QX( x, z, d );
}




/* run the scanner, returning all significant sub-intervals */
interval_stack* scanner::run()
{
    log_printf( LOG_BLAB, "running scan of length %d\n", ctx->length() );

    /* in place of recursion, subintervals left to search */
    subinterval_stack unexplored;

    /* predictions */
    subinterval_stack predictions;

    /* subinterval being considered, end inclusive */
    subinterval S;

    /* least likely subinterval in S */
    subinterval S_min;

    unexplored.push_back( subinterval( 0, ctx->length()-1 ) );

    while( !unexplored.empty() ) {
        S = unexplored.back();
        unexplored.pop_back();


        S_min = least_likely_interval( S.start, S.end, log(params->alpha) );


        /* filter out the obviously non-significant intervals */
        if( S_min.start < 0 ||
            S_min.end < 0 ||
            S_min.score >= log(params->alpha) ) continue;

        if( params->padj ) {
            S_min.score = (*params->padj)( S_min.score, S.length() ); 
        }

        if( S_min.score < params->alpha ) {
        
            predictions.push_back( S_min );

            if( S_min.start - S.start >= params->d_min ) {
                unexplored.push_back( subinterval( S.start, S_min.start-1 ) );
            }

            if( S.end - S_min.end >= params->d_min ) {
                unexplored.push_back( subinterval( S_min.end+1, S.end ) );
            }
        }
    }

    log_printf( LOG_BLAB, "finished scan of length %d\n", ctx->length() );

    interval_stack* results = new interval_stack( predictions,
                                                  ctx->get_seqname(),
                                                  ctx->get_start(),
                                                  ctx->get_strand() );

    interval_stack::iterator i;
    for( i = results->begin(); i != results->end(); i++ ) {
        ctx->adjust_interval_by_coverage( *i );
    }

    return results;
}



void scanner::conditional_push_copy(
                 subinterval_bound_pqueue& q,
                 subinterval_bound& x,
                 double score_max )
{
    if( x.count > 0 &&
        x.min_length() <= params->d_max &&
        x.max_length() >= params->d_min )
    {
        x.score = QX( x.count,
                      ctx->min_zeros( x, params->d_min ),
                      ctx->min_duration( x, params->d_min ) );

        if( x.score < score_max ) q.push( new subinterval_bound(x) );
    }
}



subinterval scanner::least_likely_interval( pos i, pos j, double alpha )
{
    log_printf( LOG_MSG, "scanning %dnt...\n", j-i+1 );
    log_indent();
    clock_t t0, t1;
    t0 = clock();

    /* resort to brute force when there are this few options */
    /* NOTE: this must be at least 1 */
    const size_t epsilon = 50;

    /* candidate bounds */
    subinterval_bound_pqueue Q;

    /* temporary values for count and rate, resp. */
    double r;
    rcount c;
    unsigned int z;

    /* least likely interval discovered so far */
    /* (initialized to pval = alpha so we don't bother looking at anything higher) */
    subinterval S_min( -1, -1 );
    S_min.score = alpha;

    if( j-i+1 <= params->d_min ) {
        log_unindent();
        return S_min;
    }

    /* current and previous candidates being considered */
    subinterval S, S_prev;

    /* current bound being considered */
    subinterval_bound* B = new subinterval_bound();

    /* divisions */
    subinterval_bound C;
    pos mid1, mid2;
    double r1, r2;
    rcount c1, c2;
    unsigned int z1, z2;



    /* number of possible subintervals in B */
    size_t m;


    ctx->count_zeros( &c, &z, i, j );
    B->set( i, j, ctx->rate(i,j), c, z );
    B->score = QX( B->count,
                   ctx->min_zeros( *B, params->d_min ),
                   ctx->min_duration( *B, params->d_min ) );

    Q.push( B );



    while( !Q.empty() ) {
        B = Q.pop(); 


        /* because these are popped in order of increasing score lower bound,
         * we can halt as soon as we see something with a hopeless lower bound.  */
        if( B->score >= S_min.score ) break;

        m = B->subinterval_count( params->d_min, params->d_max );

        /* can we use B as a new bound? */
        if( B->max_length() <= params->d_max ) {
        
            S.score = QX( B->count, B->zeros, B->max_length() );

            if( S.score < S_min.score ) {
                S_min.start = B->I_min;
                S_min.end   = B->J_max;
                S_min.count = B->count;
                S_min.rate  = B->rate;
                S_min.score = S.score;
                log_printf( LOG_BLAB, "new min (A): (i,j,score) = (%d,%d,%0.4e)\n",
                            S_min.start, S_min.end, S_min.score );
            }
        }


        /* Case 0 (Brute Force): |subintervals of B| < epsilon */
        if( m <= epsilon ) {

            if( m == 0 ) {
                delete B;
                continue;
            }

            /* initial subinterval */
            S.start = max( B->I_min, B->J_min - (params->d_max-1) );
            S.end   = min( B->J_max, S.start + (params->d_max-1) );

            /* set proper rate/count for S */
            S.count = B->count;
            S.zeros = B->zeros;
            S.rate  = B->rate;

            if( S.start > B->I_min ) {
                ctx->count_zeros( &c, &z, B->I_min, S.start-1 );
                S.count -= c;
                S.zeros -= z;
                S.rate  -= ctx->rate ( B->I_min, S.start-1 );
            }

            if( S.end < B->J_max ) {
                ctx->count_zeros( &c, &z, S.end+1, B->J_max );
                S.count -= c;
                S.zeros -= z;
                S.rate  -= ctx->rate ( S.end+1, B->J_max );
            }


            /* walk the start forward */
            while( S.start <= B->I_max ) {
                /* efficiency trick: unless we are constrained by d_min and the
                 * interval end, don't bother checking subintervals that begin on a
                 * read count of 0. We can create a subinterval at least as
                 * significant by sliding left.
                 *
                 * NOTE: this only works when there is no prior distribution on
                 * prediction lengths.
                 * */
                if( S.start + params->d_min - 1 < B->J_max &&
                    ctx->count(S.start) == 0 ) {
                    S.zeros -= 1;
                    S.rate -= ctx->rate(S.start);
                    S.start++;
                    continue;
                }

                S_prev = S;
 
                /* walk the end backwards */
                while( S.end >= B->J_min && S.length() >= params->d_min ) {
                    ctx->count_zeros( &c, &z, S.end );
                    r = ctx->rate ( S.end );

                    /* efficiency trick, as above, same conditions apply */
                    if( S.length() > params->d_min && c == 0 ) {
                        S.end--;
                        S.rate  -= r;
                        S.zeros -= z;
                        continue;
                    }

                    S.score = QX( S.count, S.zeros, S.length() );

                    if( S.score < S_min.score ) {
                        S_min = S;

                        log_printf( LOG_BLAB, "new min (B): (i,j,score) = (%d,%d,%0.4e)\n",
                                    S.start, S.end, S_min.score );
                    }

                    S.end--;
                    S.count -= c;
                    S.zeros -= z;
                    S.rate  -= r;
                }

                /* reset and move S.start forward */
                S = S_prev;
                ctx->count_zeros( &c, &z, S.start );
                S.count -= c;
                S.zeros -= z;
                S.rate  -= ctx->rate (S.start);
                S.start++;

                if( S.end < B->J_max ) {
                    S.end++;
                    ctx->count_zeros( &c, &z, S.end );
                    S.count += c;
                    S.zeros += z;
                    S.rate  += ctx->rate (S.end);
                }
            }
        }



        /* Case 1 (Bisection): I = J */
        else if( B->equal_bounds() ) {
            mid1 = B->I_min + (B->I_max - B->I_min) / 2;

            ctx->count_zeros( &c1, &z1, B->I_min, mid1 );
            r1 = ctx->rate( B->I_min, mid1 );

            /* Left: ####|---- */
            C.set( B->I_min, mid1, r1, c1, z1 );
            conditional_push_copy( Q, C, S_min.score );


            /* Right: ----|#### */
            C.set( mid1+1, B->J_max, B->rate - r1, B->count - c1, B->zeros - z1 );
            conditional_push_copy( Q, C, S_min.score );


            /* Center: ####|#### */
            C.set( B->I_min, mid1, mid1+1, B->J_max, B->rate, B->count, B->zeros );
            conditional_push_copy( Q, C, S_min.score );
        }



        /* Case 2 (Double Bisection): I, J disjoint  */
        else if( B->disjoint_bounds() )
        {
            mid1 = B->I_min + (B->I_max - B->I_min) / 2;
            mid2 = B->J_min + (B->J_max - B->J_min) / 2;

            ctx->count_zeros( &c1, &z1, B->I_min, mid1 );
            r1 = ctx->rate ( B->I_min, mid1 );

            ctx->count_zeros( &c2, &z2, mid2, B->J_min );
            r2 = ctx->rate ( mid2, B->J_max );


            /* ####----|----####  */ 
            C.set( B->I_min, mid1, mid2, B->J_max, B->rate, B->count, B->zeros );
            conditional_push_copy( Q, C, S_min.score );


            /* ----####|####----  */ 
            C.set( mid1+1, B->I_max, B->J_min, mid2-1,
                   B->rate - r1 - r2,
                   B->count - c1 - c2,
                   B->zeros - z1 - z2 );
            conditional_push_copy( Q, C, S_min.score );


            /* ####----|####----  */ 
            if( c1 > 0 || mid2 - mid1 - 1 < params->d_min  ) {
                C.set( B->I_min, mid1, B->J_min, mid2-1,
                       B->rate - r2,
                       B->count - c2,
                       B->zeros - z2 );
                conditional_push_copy( Q, C, S_min.score );
            }

            /* ----####|----####  */ 
            if( c2 > 0 || mid2 - mid1 - 1 < params->d_min ) {
                C.set( mid1+1, B->I_max, mid2, B->J_max,
                       B->rate - r1,
                       B->count - c1,
                       B->zeros - z1 );
                conditional_push_copy( Q, C, S_min.score );
            }


        }


        /* Case 3: !! ERROR !! */
        else {
            fail( "Error: LLI bounds are niether equal nor disjoint. "
                  "Please investigate/report.\n" );
        }

        delete B;
    }

    log_unindent();
    t1 = clock();
    log_printf( LOG_BLAB, "finished in %0.2f seconds\n",
                          (double)(t1-t0)/(double)CLOCKS_PER_SEC );
    return S_min;
}




