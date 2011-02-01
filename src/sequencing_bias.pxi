
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


    def train( self, ref_fn, reads_fn, n, L, R, *etc ):

        cdef double c_complexity_penalty = 1.0
        if len(etc) > 0:
            c_complexity_penalty = etc[0]

        self.cthis = train_sequencing_bias( ref_fn, reads_fn, n, L, R,
                                            c_complexity_penalty )


