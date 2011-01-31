

cdef extern from "table.h":
    ctypedef struct hashed_value:
        pass

    ctypedef struct table:
        hashed_value** A
        size_t n
        size_t m
        size_t max_m
        size_t min_m
        char** seq_names

    void table_create( table* )
    void table_destroy( table* )

    void table_inc_pos( table*, int tid, int pos, unsigned int strand )


cdef extern from "sequencing_bias.hpp":
    ctypedef struct c_sequencing_bias "sequencing_bias":
        double* get_bias( char* seqname, int start, int end, int strand )
        char* print_model_graph()
        void save_to_file( char* fn )
        pass

    c_sequencing_bias* load_sequencing_bias "new sequencing_bias" (
           char* ref_fn, char* model_fn )

    c_sequencing_bias* train_sequencing_bias "new sequencing_bias" (
           char* ref_fn, char* reads_fn, size_t n, int L, int R,
           int train_backwards, double complexity_penalty )

    c_sequencing_bias* train_sequencing_bias2 "new sequencing_bias" (
           char* ref_fn, table* T, size_t n, int L, int R,
           int train_backwards, double complexity_penalty )

    void del_sequencing_bias "delete" ( c_sequencing_bias* )



