
#ifndef PEAKOLATOR_BIAS_CORRECTION
#define PEAKOLATOR_BIAS_CORRECTION


#include "common.hpp"
#include "kmers.hpp"
#include "hash.h"

#include "samtools/sam.h"
#include "samtools/faidx.h"

class sequencing_bias
{
    public:
        sequencing_bias( const char* ref_fn,
                         const char* reads_fn,
                         size_t n, pos L, pos R,
                         double complexity_penalty = 1.0 );

        ~sequencing_bias();

        double* get_bias( const char* seqname, pos start, pos end, int strand );

        char* print_model_graph();

        sequencing_bias* copy() const;
        void clear();

    private:
        sequencing_bias();
        void build( const char* ref_fn,
                    const char* reads_fn,
                    size_t n, pos L, pos R,
                    double complexity_penality = 1.0 );


        void hash_reads( table* T, samfile_t* reads_fn,
                         size_t limit = 0 ) const;



        /* left and right sequence context */
        pos L, R;

        /* reference sequence */
        faidx_t* ref_f;
        char*    ref_fn;

        /* trained background (M0) and foreground (M1) models */
        motif* M0;
        motif* M1;

        static const double pseudocount;
};


#endif
