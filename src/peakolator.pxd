

from stdlib          cimport malloc, free
from sequencing_bias cimport *


cdef extern from "string.h":
    char* strdup( char* )


cdef extern from "samtools/bam.h":

    ctypedef struct bam_header_t:
        int           n_targets
        char**        target_name
        unsigned int* target_len
        pass

cdef extern from "samtools/sam.h":

    ctypedef struct samfile_t:
        bam_header_t* header
        pass

    samfile_t* samopen( char* fn, char* mode, void* aux )
    void samclose( samfile_t* fp )

cdef extern from "logger.h":
    void log_puts( int vl, char* msg )
    void log_printf( int vl, char* fmt, ... )
    void log_indent()
    void log_unindent()
    void log_verbosity( int )


cdef extern from "common.hpp":
    ctypedef long          pos
    ctypedef unsigned long rcount



cdef extern from "annotations.hpp":
    ctypedef struct rows:
        pass

    rows* parse_gtf( char* fn )
    rows* parse_bed( char* fn )

    void del_rows "delete" ( rows* )

    rows* get_constitutive_exons( rows*, bool )



cdef extern from "intervals.hpp":

    ctypedef struct c_interval "interval":
        void set( char* seqname, pos start, pos end, int strand )
        pos length()
        char* seqname
        pos start
        pos end
        int strand
        double score

    c_interval* new_interval "new interval" ()
    c_interval* copy_interval "new interval" ( c_interval )
    void del_interval "delete" ( c_interval* i )

    ctypedef struct c_interval_stack "interval_stack":
        int  empty()
        void pop_back()
        c_interval back()

    c_interval_stack* new_interval_stack "new interval_stack" ()
    void del_interval_stack "delete" ( c_interval_stack*  )

    void interval_stack_push ( c_interval_stack*, char* seqname, pos start, pos end, int strand )


cdef extern from "dataset.hpp":

    ctypedef struct c_dataset "dataset":
        c_dataset* copy()

        void fit_sequence_bias( char* ref_fn, size_t max_reads, pos L, pos R,
                                double complexity_penalty )
        void fit_null_distr( c_interval_stack* train, double* p, double* r )

        void hash_reads( table* T, c_interval_stack* I )
        size_t n_targets()

        c_sequencing_bias* bias



    c_dataset* new_dataset "new dataset" ( char* bam_fn )

    void del_dataset "delete" ( c_dataset* dataset )




cdef extern from "context.hpp":

    ctypedef struct c_context "context":
        void set( c_dataset* dataset, char* chrom, \
                  pos start, pos end, int strand )

        double rate()
        rcount count()

        pos length()

    c_context* new_context "new context" ()
    void del_context "delete" ( c_context* context )



cdef extern from "parameters.hpp":

    ctypedef struct c_parameters "parameters":
        c_parameters* copy()
        void rebuild_lookup( int m, int n )
        void build_padj()
        double alpha
        double r, p
        pos d_min
        pos d_max
        pass

    c_parameters* new_parameters \
                "new parameters" ()

    void del_parameters "delete" ( c_parameters* params )




cdef extern from "scanner.hpp":


    ctypedef struct c_scanner "scanner":
        c_interval_stack* run()


    c_scanner* new_scanner "new scanner" \
            ( c_parameters* params,
              c_context* context )

    void del_scanner "delete" ( c_scanner* )

