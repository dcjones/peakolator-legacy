
#ifndef PEAKOLATOR_PEAKOLATOR_DATASET
#define PEAKOLATOR_PEAKOLATOR_DATASET

#include "bias_correction.hpp"
#include "peakolator_interval.hpp"

#include "samtools/faidx.h"


class peakolator_dataset
{
    public:
        peakolator_dataset( const char* fasta_fn, const char* bam_fn,
                            pos bias_L, pos bias_R, unsigned int bias_k );

        peakolator_dataset* copy() const;

        const sequencing_bias& get_bias() const;

        void fit_null_distr( interval_stack* is, double* r, double* p );

        ~peakolator_dataset();

    private:
        peakolator_dataset();
        sequencing_bias* bias;

        /* reads file */
        samfile_t*   reads_f;
        bam_index_t* reads_index;
        char*        reads_fn;

        friend class peakolator_context;
};



#endif
