
#include "kmers.hpp"


kmer_matrix::kmer_matrix( size_t n, size_t k )
    : k(k)
{
    size_t four_to_k = 1<<(2*k);
    A = gsl_matrix_alloc( n, four_to_k );
}

kmer_matrix::kmer_matrix( const kmer_matrix& M )
{
    k = M.k;

    A = gsl_matrix_alloc( M.A->size1, M.A->size2 );
    gsl_matrix_memcpy( A, M.A );
}


void kmer_matrix::operator=( const kmer_matrix& M )
{
    k = M.k;

    if( A->size1 != M.A->size1 || A->size2 != M.A->size2 ) {
        gsl_matrix_free(A);
        A = gsl_matrix_alloc( M.A->size1,  M.A->size2 );
    }

    gsl_matrix_memcpy( A, M.A );
}


size_t kmer_matrix::n() const
{
    return A->size1;
}

size_t kmer_matrix::m() const
{
    return A->size2;
}


kmer_matrix::~kmer_matrix()
{
    gsl_matrix_free(A);
}


double kmer_matrix::get( size_t i, kmer K ) const
{
    return gsl_matrix_get( A, i, K );
}


void kmer_matrix::set( size_t i, kmer K, double x )
{
    gsl_matrix_set( A, i, K, x );
}


void kmer_matrix::inc( size_t i, kmer K, double x )
{
    gsl_matrix_set( A, i, K, gsl_matrix_get( A, i, K ) + x );
}


void kmer_matrix::setall( double x )
{
    gsl_matrix_set_all( A, x );
}


void kmer_matrix::dist_normalize()
{
    size_t i, K;
    double z;
    for( i = 0; i < A->size1; i++ ) {
        z = 0.0;
        for( K = 0; K < A->size2; K++ ) z += get( i, K );
        for( K = 0; K < A->size2; K++ ) set( i, K, get( i, K ) / z );
    }
}


void kmer_matrix::dist_conditionalize()
{
    size_t i;

    for( i = 0; i < A->size1; i++ ) {
        dist_conditionalize_row( i );
    }
}


void kmer_matrix::dist_conditionalize_row( size_t i )
{
    kmer L;
    kmer K;
    kmer nt;
    double z;

    for( L = 0; L < A->size2/4; L++ ) {
        K = L<<2;
        z = 0.0;
        for( nt = 0; nt < 4; nt++ ) {
            z += get( i, nt | K );
        }

        for( nt = 0; nt < 4; nt++ ) {
            set( i, nt | K, get( i, nt | K ) / z );
        }
    }
}


void kmer_matrix::dist_marginalize( size_t i, size_t j )
{
    /* TODO */
}




