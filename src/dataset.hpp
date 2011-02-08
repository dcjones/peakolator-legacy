
#ifndef PEAKOLATOR_PEAKOLATOR_DATASET
#define PEAKOLATOR_PEAKOLATOR_DATASET

#include "sequencing_bias.hpp"
#include "intervals.hpp"
#include "table.h"

#include "samtools/faidx.h"


class dataset
{
    public:
        dataset( const char* reads_fn );
        dataset* copy() const;

        void fit_sequence_bias( const char* ref_fn,
                                size_t max_reads, pos L, pos R,
                                double complexity_penalty = 1.0,
                                double offset_std = 10.0 );

        void load_sequence_bias( const char* ref_fn, const char* bias_fn );

        void fit_null_distr( interval_stack* is, double* r, double* p );

        size_t n_targets() const { return reads_f->header->n_targets; }

        ~dataset();

        sequencing_bias* bias;

    private:
        dataset();

        bool bias_owner;

        /* reads file */
        samfile_t*   reads_f;
        bam_index_t* reads_index;
        char*        reads_fn;

        /* indexed read counts */
        table T;
        read_counts counts;

        friend class context;
};



#endif
