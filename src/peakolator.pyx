
from sys import stderr
import re
import numpy as np

from sequencing_bias cimport *
include 'sequencing_bias.pxi'

from peakolator cimport *



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
            self.cthis.strand  = strand

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

    property score:
        def __get__( self ): return self.cthis.score
        def __set__( self, score ): self.cthis.score = score



## dataset
cdef class dataset:
    cdef c_dataset* cthis

    def __cinit__( self ):
        self.cthis = NULL

    def __dealloc__( self ):
        del_dataset( self.cthis )

    def init( self, reads_fn, sequencing_bias bias = None ):

        if self.cthis:
            del_dataset( self.cthis )

        self.cthis = new_dataset( reads_fn )

        if bias is not None:
            self.cthis.bias = bias.cthis


    def get_bias( self, chrom, start, end, strand ):
        if self.cthis.bias == NULL: return None

        cdef double* c_ws

        c_ws = self.cthis.bias.get_bias( chrom, start, end, strand )

        cdef int i
        ws = np.empty( end - start + 1, dtype=np.double )
        for i in range(end-start+1):
            ws[i] = c_ws[i]

        return ws


    def print_model_graph( self ):

        cdef char* graph_cstr

        graph_cstr = self.cthis.bias.print_model_graph()
        graph_str = str(graph_cstr)

        free(<void*>graph_cstr)

        return graph_str


    def copy( self ):
        cdef dataset ds = dataset()
        ds.cthis = self.cthis.copy()
        return ds


    def fit_sequence_bias( self, ref_fn, max_reads, L, R, complexity_penalty = 1.0 ):
        self.cthis.fit_sequence_bias( ref_fn, max_reads, L, R, complexity_penalty )


    def load_sequence_bias( self, ref_fn, bias_fn ):
        self.cthis.load_sequence_bias( ref_fn, bias_fn )


    def fit_null_distr( self, train ):

        # convert to an interval_stack
        cdef c_interval_stack* IS = new_interval_stack()

        for I in train:
            interval_stack_push( IS, I.seqname, I.start, I.end, I.strand )


        cdef double r, p, a
        self.cthis.fit_null_distr( IS, &r, &p, &a )

        del_interval_stack( IS )

        return (float(r),float(p),float(a))





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

    property a:
        def __get__(self):   return self.cthis.a
        def __set__(self,a): self.cthis.a = a

    property d_min:
        def __get__(self):      return self.cthis.d_min
        def __set__(self,d_min): self.cthis.d_min = d_min

    property d_max:
        def __get__(self):      return self.cthis.d_max
        def __set__(self,d_max): self.cthis.d_max = d_max



## scanner 
cdef class scanner:
    cdef c_scanner* cthis

    def __cinit__( self, parameters prm, context ctx ):
        self.cthis = new_scanner( prm.cthis, ctx.cthis )

    def __dealloc__( self ):
        del_scanner( self.cthis )

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



## a sequence motif representation using the bias correction machinery
cdef class motif:
    cdef c_sequencing_bias* cthis

    def __cinit__( self, ref_fn, positions, L = 5, R = 5 ):
        cdef table T

        tids = {}
        for (seqname, pos, strand) in positions:
            if seqname not in tids: tids[seqname] = len(tids)

        table_create( &T, len(tids) )
        table_inc_pos( &T, tids[seqname], pos, strand )

        cdef char** seq_names = <char**>malloc( len(tids) * sizeof(char*) )
        for (seqname,tid) in tids.iteritems():
            seq_names[tid] = strdup(seqname)
        T.seq_names = seq_names

        self.cthis = train_sequencing_bias2( ref_fn, &T, 0, L, R, 1.0, 10.0 )

        table_destroy( &T )

        for tid in xrange(len(tids)):
            free( <void*>seq_names[tid] )
        free(<void*>seq_names)


    def print_model_graph( self ):
        cdef char* c_dot_str
        c_dot_str = self.cthis.print_model_graph()
        dot_str = <str>c_dot_str
        stdlib.free( <void*>c_dot_str )

        return dot_str

    def __dealloc__( self ):
        del_sequencing_bias( self.cthis )



def lpdnbsum( x, r, p, d ):
    return c_lpdnbsum( x, r, p, d, False )

def lddnbsum( x, r, p, d ):
    return c_lddnbsum( x, r, p, d )



