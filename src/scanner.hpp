
#ifndef PEAKOLATOR_PEAKOLATOR_SCANNER
#define PEAKOLATOR_PEAKOLATOR_SCANNER

#include "common.hpp"
#include "context.hpp"
#include "intervals.hpp"
#include "parameters.hpp"

#include <cstdlib>

#include "samtools/sam.h"



class emppval;


class scanner
{
    public:
        scanner( 
                parameters* params,
                context*    context );

        interval_stack* run();
        subinterval least_likely_interval( pos l, pos r, double alpha );

        ~scanner();


    private:
        /* raw p-value (score) for subinterval of duration d with read count x
         * and z positions with a count of zero. */
        double QX( rcount x, unsigned int z, unsigned int d );

        /* push a new bound onto the priority queue if it is valid
         * (i.e. duration within [d_min,d_max] and score lower bound does not rule it out.) */
        void conditional_push_copy(
                 subinterval_bound_pqueue& q,
                 subinterval_bound& x,
                 double score_max );

        context*    ctx;
        parameters* params;

        emppval* padj;
};




#endif

