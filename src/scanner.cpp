
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

    if( !this->params->dist.ready() ) params->dist.build( params->r, params->p );

}


scanner::~scanner()
{
}



mpfr_class scanner::QX( double r, rcount x )
{
    return params->dist.QX( r, x );
}




/* run the scanner, returning all significant sub-intervals */
interval_stack* scanner::run()
{
    log_printf( LOG_MSG, "running scan of length %d\n", ctx->length() );

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


        S_min = least_likely_interval( S.start, S.end, params->alpha );


        /* filter out the obviously non-significant intervals */
        if( S_min.pval >= params->alpha ) continue;

        if( params->padj ) {
            S_min.pval = params->padj->adjust( S_min.pval, S.length() ); 
        }

        if( S_min.pval < params->alpha ) {
        
            predictions.push_back( S_min );

            if( S_min.start - S.start >= params->d_min ) {
                unexplored.push_back( subinterval( S.start, S_min.start-1 ) );
            }

            if( S.end - S_min.end >= params->d_min ) {
                unexplored.push_back( subinterval( S_min.end+1, S.end ) );
            }
        }
    }

    log_printf( LOG_MSG, "finished scan of length %d\n", ctx->length() );

    return new interval_stack( predictions, ctx->seqname, ctx->start, ctx->strand );
}



void scanner::conditional_push_copy(
                 subinterval_bound_pqueue& q,
                 subinterval_bound& x,
                 const mpfr_class& p_max )
{
    if( x.count > 0 &&
        x.min_length() <= params->d_max &&
        x.max_length() >= params->d_min )
    {
        x.pval = QX(
                     ctx->min_rate( x, params->d_min ),
                     x.count );
        if( x.pval < p_max ) q.push( new subinterval_bound(x) );
    }
}



subinterval scanner::least_likely_interval( pos i, pos j, double alpha )
{
    log_printf( LOG_BLAB, "scanning %dnt...\n", j-i+1 );
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

    /* least likely interval discovered so far */
    /* (initialized to pval = alpha so we don't bother looking at anything higher) */
    subinterval S_min( -1, -1 );
    S_min.pval = alpha;

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


    /* number of possible subintervals in B */
    size_t m;

    B->set( i, j, ctx->rate(i,j), ctx->count(i,j) );
    B->pval = QX(  ctx->min_rate( *B, params->d_min ), B->count );

    Q.push( B );



    while( !Q.empty() ) {
        B = Q.pop(); 


        /* because these are popped in order of increasing p-value lower bound,
         * we can halt as soon as we see something with a hopeless lower bound.  */
        if( B->pval >= S_min.pval ) break;

        //char* tmp = mpfr_to_string( B->pval );
        //log_printf( "popped with p-val = %s", tmp );
        //free(tmp);

        m = B->subinterval_count( params->d_min, params->d_max );

        /* can we use B as a new bound? */
        if( B->max_length() <= params->d_max ) {
            S.pval = QX( B->rate, B->count );
            if( S.pval < S_min.pval ) {
                S_min.start = B->I_min;
                S_min.end   = B->J_max;
                S_min.count = B->count;
                S_min.rate  = B->rate;
                S_min.pval  = S.pval;
                char* tmp = mpfr_to_string( S_min.pval );
                log_printf( LOG_BLAB, "new min (A): (i,j,pval) = (%d,%d,%s)\n",
                            S_min.start, S_min.end, tmp );
                free(tmp);
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
            S.rate  = B->rate;

            if( S.start > B->I_min ) {
                S.count -= ctx->count( B->I_min, S.start-1 );
                S.rate  -= ctx->rate ( B->I_min, S.start-1 );
            }

            if( S.end < B->J_max ) {
                S.count -= ctx->count( S.end+1, B->J_max );
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
                    S.rate -= ctx->rate(S.start);
                    S.start++;
                    continue;
                }

                S_prev = S;
 
                /* walk the end backwards */
                while( S.end >= B->J_min && S.length() >= params->d_min ) {
                    c = ctx->count( S.end );
                    r = ctx->rate ( S.end );

                    /* efficiency trick, as above, same conditions apply */
                    if( S.length() > params->d_min && c == 0 ) {
                        S.end--;
                        S.rate -= r;
                        continue;
                    }

                    S.pval = QX( S.rate, S.count );

                    if( S.pval < S_min.pval ) {
                        S_min = S;

                        char* tmp = mpfr_to_string( S_min.pval );
                        log_printf( LOG_BLAB, "new min (B): (i,j,pval) = (%d,%d,%s)\n",
                                    S.start, S.end, tmp );
                        free(tmp);
                    }

                    S.end--;
                    S.count -= c;
                    S.rate  -= r;
                }

                /* reset and move S.start forward */
                S = S_prev;
                S.count -= ctx->count(S.start);
                S.rate  -= ctx->rate (S.start);
                S.start++;

                if( S.end < B->J_max ) {
                    S.end++;
                    S.count += ctx->count(S.end);
                    S.rate  += ctx->rate (S.end);
                }
            }
        }



        /* Case 1 (Bisection): I = J */
        else if( B->equal_bounds() ) {
            mid1 = B->I_min + (B->I_max - B->I_min) / 2;

            c1 = ctx->count( B->I_min, mid1 );
            r1 = ctx->rate ( B->I_min, mid1 );

            /* Left: ####|---- */
            C.set( B->I_min, mid1, r1, c1 );
            conditional_push_copy( Q, C, S_min.pval );


            /* Right: ----|#### */
            C.set( mid1+1, B->J_max, B->rate - r1, B->count - c1 );
            conditional_push_copy( Q, C, S_min.pval );


            /* Center: ####|#### */
            C.set( B->I_min, mid1, mid1+1, B->J_max, B->rate, B->count );
            conditional_push_copy( Q, C, S_min.pval );
        }



        /* Case 2 (Double Bisection): I, J disjoint  */
        else if( B->disjoint_bounds() )
        {
            mid1 = B->I_min + (B->I_max - B->I_min) / 2;
            mid2 = B->J_min + (B->J_max - B->J_min) / 2;

            c1 = ctx->count( B->I_min, mid1 );
            r1 = ctx->rate ( B->I_min, mid1 );

            c2 = ctx->count( mid2, B->J_max );
            r2 = ctx->rate ( mid2, B->J_max );


            /* ####----|----####  */ 
            C.set( B->I_min, mid1, mid2, B->J_max, B->rate, B->count );
            conditional_push_copy( Q, C, S_min.pval );


            /* ----####|####----  */ 
            C.set( mid1+1, B->I_max, B->J_min, mid2-1, B->rate - r1 - r2, B->count - c1 - c2 );
            conditional_push_copy( Q, C, S_min.pval );


            /* ####----|####----  */ 
            if( c1 > 0 || mid2 - mid1 - 1 < params->d_min  ) {
                C.set( B->I_min, mid1, B->J_min, mid2-1, B->rate - r2, B->count - c2 );
                conditional_push_copy( Q, C, S_min.pval );
            }

            /* ----####|----####  */ 
            if( c2 > 0 || mid2 - mid1 - 1 < params->d_min ) {
                C.set( mid1+1, B->I_max, mid2, B->J_max, B->rate - r1, B->count - c1 );
                conditional_push_copy( Q, C, S_min.pval );
            }


        }


        /* Case 3: !! ERROR !! */
        else {
            log_puts( LOG_ERROR, "Error: LLI bounds are niether equal nor disjoint. "
                      "Please investigate/report.\n" );
            exit(1);
        }

        delete B;
    }

    log_unindent();
    t1 = clock();
    log_printf( LOG_BLAB, "finished in %0.2f seconds\n",
                          (double)(t1-t0)/(double)CLOCKS_PER_SEC );
    return S_min;
}




