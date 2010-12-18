

cdef extern from "sequencing_bias.hpp":
    ctypedef struct c_sequencing_bias "sequencing_bias":
        double* get_bias( char* seqname, int start, int end, int strand )
        char* print_model_graph()
        void save_to_file( char* fn )
        pass

    c_sequencing_bias* load_sequencing_bias "new sequencing_bias" (
           char* ref_fn, char* reads_fn, char* model_fn )

    c_sequencing_bias* train_sequencing_bias "new sequencing_bias" (
           char* ref_fn, char* reads_fn, size_t n, int L, int R,
           int train_backwards, double complexity_penalty )

    void del_sequencing_bias "delete" ( c_sequencing_bias* )



