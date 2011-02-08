
cimport stdlib

cdef class sequencing_bias:
    cdef c_sequencing_bias* cthis

    def __cinit__( self ):
        self.cthis = NULL


    def __dealloc__( self ):
        del_sequencing_bias( self.cthis )


    def load( self, ref_fn, model_fn ):
        if self.cthis:
            del_sequencing_bias( self.cthis )
        self.cthis = load_sequencing_bias( ref_fn, model_fn )


    def save( self, fn ):
        self.cthis.save_to_file( fn )


    def print_model_graph( self ):
        cdef char* c_dot_str
        c_dot_str = self.cthis.print_model_graph()
        dot_str = <str>c_dot_str
        stdlib.free( <void*>c_dot_str )

        return dot_str


    def train( self, ref_fn, reads_fn, n, L, R, whitelist = None,
               complexity_penalty = 1.0, offset_std = 10.0 ):

        cdef double c_complexity_penalty = complexity_penalty
        cdef double c_offset_std         = offset_std


        cdef c_interval_stack* IS = NULL
        cdef table T

        if whitelist is not None:
            IS = new_interval_stack()
            for I in whitelist:
                interval_stack_push( IS, I.seqname, I.start, I.end, I.strand )

            hash_reads( &T, reads_fn, IS )

            self.cthis = train_sequencing_bias2( ref_fn, &T, n, L, R,
                                                 c_complexity_penalty, c_offset_std )

            table_destroy( &T )
            del_interval_stack( IS )
        else:
            self.cthis = train_sequencing_bias( ref_fn, reads_fn, n, L, R,
                                                c_complexity_penalty, c_offset_std )





