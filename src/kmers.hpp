
#ifndef PEAKOLATOR_KMERS
#define PEAKOLATOR_KMERS

#include <boost/cstdint.hpp>
#include <gsl/gsl_matrix.h>

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


        /* normalize to turn each position into a proper distribution over kmers
         * */
        void dist_normalize();
        void dist_marginalize( size_t i, size_t j );

        void dist_conditionalize();

    private:
        void dist_conditionalize_row( size_t i );

        size_t k;
        gsl_matrix* A;
};


#endif


