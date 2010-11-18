
#ifndef PEAKOLATOR_KMERS
#define PEAKOLATOR_KMERS

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

        void dist_conditionalize();
        void dist_conditionalize_row( size_t i );

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
        sequence( const char* s );
        sequence( const sequence& );
        void operator=( const sequence& );
        ~sequence();

        kmer get( size_t i ) const;
        kmer get( const std::set<size_t>& ) const;


    private:
        kmer* xs;
        size_t n;

        static const size_t kmer_max_k;
};





/* A bayesian network representing sequence probability. */


class motif
{
    public:
        motif( size_t n, size_t k, const std::deque<sequence*>* data );
        ~motif();

        void add_edge( size_t i, size_t j );
        void rem_edge( size_t i, size_t j );

        double eval( const sequence& );

        size_t num_params() const;

        void store_row( size_t i );
        void restore_stored_row();

    private:
        size_t n;
        size_t k;
        kmer_matrix* P;

        const std::deque<sequence*>* data;
        std::set<size_t>* parents;

        static const double pseudocount;

        friend double motif_log_likelihood( motif& M0, motif& M1 );
        friend void train_motifs( motif& M0, motif& M1 );
};


void train_motifs( motif& M0, motif& M1 );
                   



#endif





