
#ifndef PEAKOLATOR_PEAKOLATOR_DATASET
#define PEAKOLATOR_PEAKOLATOR_DATASET

#include "sequencing_bias.hpp"
#include "intervals.hpp"
#include "table.h"

#include "samtools/faidx.h"


class dataset
{
    public:
        dataset( const char* reads_fn, sequencing_bias* sb );
        dataset* copy() const;

        void fit_null_distr( interval_stack* is, double* r, double* p );
        void hash_reads( table* T, interval_stack* is );

        size_t n_targets() const { return reads_f->header->n_targets; }

        ~dataset();

        sequencing_bias* bias;

    private:
        dataset();

        /* reads file */
        samfile_t*   reads_f;
        bam_index_t* reads_index;
        char*        reads_fn;

        friend class context;
};



#endif
