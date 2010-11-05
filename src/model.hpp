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
#include "context.hpp"
#include "intervals.hpp"
#include "parameters.hpp"

#include <cstdlib>

#include "samtools/sam.h"



class emppval;


class model
{
    public:
        model( 
                parameters* params,
                context*    context );

        interval_stack* run();
        subinterval least_likely_interval   ( pos l, pos r, double alpha );

        ~model();

        /* raw p-value for subinterval of duration d with read count x */
        mpfr_class QX( double r, rcount x );

        context*    ctx;
        parameters* params;

    private:
        void     trim_candidate( subinterval& I );
        emppval* padj;
};




#endif

