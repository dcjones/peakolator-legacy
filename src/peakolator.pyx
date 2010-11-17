
from sys import stderr
import re
import numpy as np


cdef extern from "stdlib.h":
    void free( void* )

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

    ctypedef struct mpfr_class:
        pass

    char* mpfr_to_string( mpfr_class )

cdef extern from "annotations.hpp":
    ctypedef struct rows:
        pass

    rows* parse_gtf( char* fn )
    rows* parse_bed( char* fn )

    void del_rows "delete" ( rows* )

    rows* get_constitutive_exons( rows*, bool )


cdef extern from "sequencing_bias.hpp":
    ctypedef struct c_sequencing_bias "sequencing_bias":
        double* get_bias( char* seqname, pos start, pos end, int strand )
        pass


cdef extern from "intervals.hpp":

    ctypedef struct c_interval "interval":
        void set( char* seqname, pos start, pos end, int strand )
        pos length()
        char* seqname
        pos start
        pos end
        int strand
        mpfr_class pval

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
        void fit_null_distr( c_interval_stack* train, double* p, double* r )
        c_sequencing_bias* bias



    c_dataset* new_dataset "new dataset" ( \
            char* fasta_fn, char* bam_fn, \
            pos bias_L, pos bias_R, unsigned int bias_k,
            bool count_dups, double q,
            char* training_seqname )

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


cdef extern from "model.hpp":

    ctypedef struct c_parameters "parameters":
        c_parameters* copy()
        void rebuild_lookup( int m, int n )
        void build_padj()
        double alpha
        double r, p
        pos d_min
        pos d_max
        size_t n_mc
        pos padj_spacing
        int padj_n
        pass

    c_parameters* new_parameters \
                "new parameters" ()

    void del_parameters "delete" ( c_parameters* params )



    ctypedef struct c_model "model":
        c_interval_stack* run()


    c_model* new_model "new model" \
            ( c_parameters* params,
              c_context* context )

    void del_model "delete" ( c_model* model )




## logger

LOG_ERROR=0
LOG_WARN =1
LOG_MSG  =2
LOG_BLAB =3


cdef class logger:
    def __cinit__( self ):
        pass

    def __dealloc__( self ):
        pass

    def indent( self ):
        log_indent()

    def unindent( self ):
        log_unindent()

    def write( self, msg, level = LOG_MSG ):
        log_puts( level, msg )

    def verbosity( self, level ):
        log_verbosity( level )



## annotations

cdef class annotations:
    cdef rows* rs

    def __cinit__( self ):
        self.rs = NULL


    def __dealloc__( self ):
        del_rows( self.rs )


    def clear( self ):
        del_rows( self.rs )
        self.rs = NULL

    def get_constitutive_exons( self, stranded = True ):
        cdef annotations b = annotations()
        b.rs = get_constitutive_exons( self.rs, stranded )
        return b

    def parse( self, fn ):
        if re.search( '\.bed$', fn ):   self.parse_gtf( fn )
        elif re.search( '\.gtf$', fn ): self.parse_bed( fn )

        self.parse_bed( fn )


    def parse_gtf( self, fn ):
        self.clear()
        self.rs = parse_gtf( fn )


    def parse_bed( self, fn ):
        self.clear()
        self.rs = parse_bed( fn )




## interval
cdef class interval:
    cdef c_interval* cthis

    def __cinit__( self, seqname = None, start = None, end = None, strand = None ):
        self.cthis = new_interval()
        self.set( seqname, start, end, strand )

    def __dealloc__( self ):
        del_interval( self.cthis )

    def set( self, seqname = None, start = None, end = None, strand = None ):
        if seqname != None:
            free(<void*>self.cthis.seqname)
            self.cthis.seqname = strdup(seqname)

        if start   != None: self.cthis.start   = start
        if end     != None: self.cthis.end     = end
        if strand  != None:
            assert strand in (-1,0,1)
            strand  = strand

    property seqname:
        def __get__( self ):
            if self.cthis.seqname != NULL:
                return self.cthis.seqname
            else:
                return None

        def __set__( self, seqname ): self.set( seqname = seqname )

    property start:
        def __get__( self ): return self.cthis.start
        def __set__( self, start ): self.set( start = start )

    property end:
        def __get__( self ): return self.cthis.end
        def __set__( self, end ): self.set( end = end )

    property strand:
        def __get__( self ): return self.cthis.strand
        def __set__( self, strand ): self.set( strand = strand )

    property pval:
        def __get__( self ):
            cdef char* s = mpfr_to_string( self.cthis.pval )
            py_s = str(s)
            free(<void*>s)

            return py_s



## dataset
cdef class dataset:
    cdef c_dataset* cthis

    def __cinit__( self, *args ):
        cdef char* fasta_fn_cstr    = NULL
        cdef char* training_seqname = NULL
        cdef double q = 0.1
        cdef bool count_dups = True

        if len(args) >= 5:
            (fasta_fn,bam_fn,bias_L,bias_R,bias_k) = args[:5]

            if fasta_fn is not None:
                fasta_fn_cstr = fasta_fn

            if len(args) > 5:
                count_dups = args[5]

            if len(args) > 6:
                q = args[6]

            if len(args) > 7:
                training_seqname = args[7]

            self.cthis = new_dataset(
                            fasta_fn_cstr, \
                            bam_fn, \
                            bias_L, bias_R, bias_k,
                            count_dups, q,
                            training_seqname )
        else:
            self.cthis = (<dataset>args[0]).cthis.copy()


    def get_bias( self, chrom, start, end, strand ):
        cdef double* c_ws
        c_ws = self.cthis.bias.get_bias( chrom, start, end, strand )

        cdef int i
        ws = np.empty( end - start + 1, dtype=np.double )
        for i in range(end-start+1):
            ws[i] = c_ws[i]

        return ws


    def __dealloc__( self ):
        del_dataset( self.cthis )

    def copy( self ):
        return dataset( self )

    def fit_null_distr( self, train ):

        # convert to an interval_stack
        cdef c_interval_stack* IS = new_interval_stack()

        for I in train:
            interval_stack_push( IS, I.seqname, I.start, I.end, I.strand )


        cdef double r, p
        self.cthis.fit_null_distr( IS, &r, &p )

        del_interval_stack( IS )

        return (float(r),float(p))





## context
cdef class context:
    cdef c_context* cthis

    def __cinit__( self ):
        self.cthis = new_context()

    def __dealloc__( self ):
        if self.cthis != NULL:
            del_context( self.cthis )

    def rate( self ):
        return self.cthis.rate()

    def count( self ):
        return self.cthis.count()

    def set( self, dataset d, chrom, start, end, strand ):
        self.cthis.set( d.cthis, chrom, start, end, strand )

    def length( self ):
        return self.cthis.length()


## parameters
cdef class parameters:
    cdef c_parameters* cthis

    def copy( self ):
        return parameters( self )

    def rebuild_lookup( self, m = 0, n = 0 ):
        self.cthis.rebuild_lookup( m, n )

    def build_padj( self ):
        self.cthis.build_padj()

    def __cinit__( self, parameters other=None ):
        if other:
            self.cthis = other.cthis.copy()
        else:
            self.cthis = new_parameters()

    def __dealloc__( self ):
        del_parameters( self.cthis )

    property alpha:
        def __get__(self):       return self.cthis.alpha
        def __set__(self,alpha): self.cthis.alpha = alpha

    property r:
        def __get__(self):   return self.cthis.r
        def __set__(self,r): self.cthis.r = r

    property p:
        def __get__(self):   return self.cthis.p
        def __set__(self,p): self.cthis.p = p

    property d_min:
        def __get__(self):      return self.cthis.d_min
        def __set__(self,d_min): self.cthis.d_min = d_min

    property d_max:
        def __get__(self):      return self.cthis.d_max
        def __set__(self,d_max): self.cthis.d_max = d_max

    property n_mc:
        def __get__(self):      return self.cthis.n_mc
        def __set__(self,n_mc): self.cthis.n_mc = n_mc

    property padj_spacing:
        def __get__(self):              return self.cthis.padj_spacing
        def __set__(self,padj_spacing): self.cthis.padj_spacing = padj_spacing

    property padj_n:
        def __get__(self):        return self.cthis.padj_n
        def __set__(self,padj_n): self.cthis.padj_n = padj_n




## model
cdef class model:
    cdef c_model* cthis

    def __cinit__( self, parameters prm, context ctx ):
        self.cthis = new_model( prm.cthis, ctx.cthis )

    def __dealloc__( self ):
        del_model( self.cthis )

    def run( self ):
        cdef c_interval_stack* IS = self.cthis.run()
        cdef interval i

        predictions = []

        # convert the output to python wrapped classes
        while not IS.empty():
            i = interval()
            del_interval( i.cthis )
            i.cthis = copy_interval( IS.back() )

            predictions.append( i )

            IS.pop_back()

        del_interval_stack( IS )
        return predictions


## extract the chromosome sizes from a BAM file header
def get_chrom_sizes_from_bam( bam_fn ):
    cdef samfile_t* fp
    fp = samopen( bam_fn, "rb", NULL )
    if fp == NULL:
        stderr.write( 'Can\'t open open BAM file "%s".', bam_fn )
        exit(1)

    chrom_sizes = {}
    for i in range(fp.header.n_targets):
        chrom_sizes[ fp.header.target_name[i] ] = fp.header.target_len[i]

    samclose(fp)

    return chrom_sizes


