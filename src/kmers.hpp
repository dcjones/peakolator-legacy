
#ifndef PEAKOLATOR_KMERS
#define PEAKOLATOR_KMERS

#include "yaml-cpp/yaml.h"
#include <boost/cstdint.hpp>
#include <gsl/gsl_matrix.h>
#include <deque>
#include <set>

/* kmers are encoded in 16 bits, allowing for k <= 8 */
typedef boost::uint16_t kmer;


/*
 * A class that stores a matrix over kmers, of the form
 *
 *        AAA AAT AAC AAG ATA ...
 * pos 1  a11 a12 a13 ...
 * pos 2  a21 a22 a23 ...
 * pos 3  a31 a32 a33 ...
 *  ...   ...
 */
class kmer_matrix
{
    public:
        kmer_matrix( size_t n, size_t k );
        kmer_matrix( const kmer_matrix& );
        void operator=( const kmer_matrix& );

        void to_yaml( YAML::Emitter& out ) const;

        size_t n() const;
        size_t m() const;

        ~kmer_matrix();
        double get( size_t i, kmer K ) const;
        void   set( size_t i, kmer K, double x );
        void   inc( size_t i, kmer K, double x = 1.0 );
        void   setall( double x );
        void   setrow( size_t i, double x );



        /* normalize to turn each position into a proper distribution over kmers
         * */
        void dist_normalize();
        void dist_normalize_row( size_t i );
        void dist_marginalize( size_t i, size_t j );

        void dist_conditionalize( int effective_k = -1 );
        void dist_conditionalize_row( size_t i, int effective_k = -1 );
        void log_transform_row( size_t i, int effective_k = -1 );

        /* allows one row to be stored then reset, which is used when searching
         * for the optimal edge to add when training the model */
        void store_row( size_t i );
        void restore_stored_row();


    private:

        size_t k;
        gsl_matrix* A;

        gsl_vector* stored_row;
        size_t stored_row_index;
};




/*
 * Represent a sequence in 2bit encoding, extract kmers, etc. 
 */
class sequence
{
    public:
        sequence( const char* s, int meta = 0 );
        sequence( const sequence& );
        void operator=( const sequence& );
        ~sequence();

        kmer get( size_t i ) const;
        bool get( const bool* indexes, size_t maxn, kmer& K, size_t offset = 0 ) const;

        int meta;

    private:
        kmer* xs;
        size_t n;

        static const size_t kmer_max_k;
};





/* A 'bayesian' network representing sequence probability. */


class motif
{
    public:
        motif( size_t n, size_t k, int meta );
        motif( const motif& );
        ~motif();

        void to_yaml( YAML::Emitter& ) const;

        void add_edge( size_t i, size_t j, const std::deque<sequence*>* data );
        void remove_edge( size_t i, size_t j, const std::deque<sequence*>* data );
        void add_all_edges( const std::deque<sequence*>* data );

        double eval( const sequence&, size_t offset = 0 ) const;
        double eval_node( size_t i, const std::deque<sequence*>* data,
                          size_t offset = 0 ) const;

        size_t num_params() const;

        void store_row( size_t i );
        void restore_stored_row();
        char* print_model_graph( int offset = 0 );

        int meta; /* which subset of the training data to consider */

    private:

        size_t num_parents( size_t i ) const;
        bool has_edge( size_t i, size_t j );
        void set_edge( size_t i, size_t j, bool );


        void update_likelihood_column( gsl_matrix* L, size_t j,
                                       const std::deque<sequence*>* training_seqs );

        size_t n; /* number of positions */
        size_t k; /* maximum number of edges */
        kmer_matrix* P;

        bool* parents;

        static const double pseudocount;


        friend void train_motifs( motif& M0, motif& M1,
                                  const std::deque<sequence*>* training_seqs,
                                  size_t max_dep_dist, double complexity_penalty );

        friend void train_motifs_backwards( motif& M0, motif& M1,
                                            const std::deque<sequence*>* training_seqs,
                                            size_t max_dep_dist, double complexity_penalty );
};

void train_motifs( motif& M0, motif& M1,
                   const std::deque<sequence*>* training_seqs,
                   size_t max_dep_dist = 0, double complexity_penalty = 1.0 );

void train_motifs_backwards( motif& M0, motif& M1,
                             const std::deque<sequence*>* training_seqs,
                             size_t max_dep_dist = 0, double complexity_penalty = 1.0 );
                   



#endif


