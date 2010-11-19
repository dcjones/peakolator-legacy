
#include "kmers.hpp"
#include "common.hpp"
#include "logger.h"

#include <string.h>
#include <gsl/gsl_math.h>


kmer_matrix::kmer_matrix( size_t n, size_t k )
    : k(k)
{
    size_t four_to_k = 1<<(2*k);
    A = gsl_matrix_alloc( n, four_to_k );
    stored_row = gsl_vector_alloc( four_to_k );
}

kmer_matrix::kmer_matrix( const kmer_matrix& M )
{
    k = M.k;

    A = gsl_matrix_alloc( M.A->size1, M.A->size2 );
    stored_row = gsl_vector_alloc( M.stored_row->size );

    gsl_matrix_memcpy( A, M.A );
    gsl_vector_memcpy( stored_row, M.stored_row );
}


void kmer_matrix::operator=( const kmer_matrix& M )
{
    k = M.k;

    if( A->size1 != M.A->size1 || A->size2 != M.A->size2 ) {
        gsl_matrix_free(A);
        A = gsl_matrix_alloc( M.A->size1,  M.A->size2 );

        gsl_vector_free(stored_row);
        stored_row = gsl_vector_alloc( M.stored_row->size );
    }

    gsl_matrix_memcpy( A, M.A );
    gsl_vector_memcpy( stored_row, M.stored_row );
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
    gsl_vector_free(stored_row);
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

void kmer_matrix::setrow( size_t i, double x )
{
    size_t j;
    for( j = 0; j < A->size2; j++ ) {
        gsl_matrix_set( A, i, j, x );
    }
}


void kmer_matrix::dist_normalize()
{
    size_t i;
    for( i = 0; i < A->size1; i++ ) {
        dist_normalize_row( i );
    }
}

void kmer_matrix::dist_normalize_row( size_t i )
{
    double z = 0.0;
    kmer K;
    for( K = 0; K < A->size2; K++ ) z += get( i, K );
    for( K = 0; K < A->size2; K++ ) set( i, K, get( i, K ) / z );
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


void kmer_matrix::store_row( size_t i )
{
    gsl_matrix_get_row( stored_row, A, i );
    stored_row_index = i;
}


void kmer_matrix::restore_stored_row()
{
    gsl_matrix_set_row( A, stored_row_index, stored_row );
}


const size_t sequence::kmer_max_k = 4*sizeof(kmer);


sequence::sequence( const char* s )
    : xs(NULL), n(0)
{
    n = strlen(s);
    if( n > 0 ) {
        xs = new kmer[ n/kmer_max_k + 1 ];
        memset( xs, 0, (n/kmer_max_k + 1)*sizeof(kmer) );

        size_t i,j;
        for( i = 0; i < n; i++ ) {
            j = i / kmer_max_k;
            xs[j] = (xs[j] << 2) | nt2num( s[i] );
        }
    }
}

sequence::sequence( const sequence& s )
    : xs(NULL), n(0)
{
    n = s.n;
    if( n > 0 ) {
        xs = new kmer[ n/kmer_max_k + 1 ];
        memcpy( xs, s.xs, (n/kmer_max_k + 1)*sizeof(kmer) );
    }
}


void sequence::operator=( const sequence& s )
{
    delete[] xs;
    n = s.n;
    xs = new kmer[ n/kmer_max_k + 1 ];
    memcpy( xs, s.xs, (n/kmer_max_k + 1)*sizeof(kmer) );
}


sequence::~sequence()
{
    delete[] xs;
}

kmer sequence::get( size_t i ) const
{
    size_t j = i / (sizeof(kmer)/2);
    size_t l = i % (sizeof(kmer)/2);
    return (xs[j] >> l) & 0x3;
}


bool sequence::get( const bool* indexes, size_t maxn, kmer& K ) const
{
    bool nonempty = false;
    size_t j, l;
    K = 0;
    size_t i;
    for( i = 0; i < maxn; i++ ) {
        if( indexes[i] ) {
            nonempty = true;
            j = i / kmer_max_k;
            l = i % kmer_max_k;
            K = (K<<2) | ((xs[j] >> l) & 0x3);
        }
        i++;
    }

    return nonempty;
}


const double motif::pseudocount = 1;

motif::motif( size_t n, size_t k, const std::deque<sequence*>* data )
    : n(n), k(k), data(data)
{
    P = new kmer_matrix( n, k );
    P->setall( 1.0 );

    parents = new bool[n*n];
    memset( parents, 0, n*n*sizeof(bool) );
}

motif::~motif()
{
    delete[] parents;
    delete P;
}


double motif::eval( const sequence& seq )
{
    double l = 0.0; /* log likelihood */

    const size_t n = P->n();
    size_t i;
    kmer K;
    for( i = 0; i < n; i++ ) {
        if( !seq.get( parents + i*n, n, K ) ) continue;
        l += log( P->get( i, K ) );
    }

    return l;
}


size_t motif::num_params() const
{
    size_t N = 0;
    size_t i;
    for( i = 0; i < n; i++ ) {
        N += (1 << (2*num_parents(i))) - 1;
    }

    return N;
}

size_t motif::num_parents( size_t i ) const
{
    size_t j;
    size_t M = 0;
    for( j = 0; j <= i; j++ ) {
        if( parents[i*n+j] ) M++;
    }

    return M;
}

bool motif::has_edge( size_t i, size_t j )
{
    return parents[j*n+i];
}

void motif::set_edge( size_t i, size_t j, bool x )
{
    parents[j*n+i] = true;
}


void motif::store_row( size_t i )
{
    P->store_row( i );
}


void motif::restore_stored_row()
{
    P->restore_stored_row();
}



/* make an edge i --> j
 * That is, condition j on i. */
void motif::add_edge( size_t i, size_t j )
{
    if( i > j ) {
        log_printf( LOG_ERROR, "Invalid motif edge (%zu, %zu)\n", i, j );
        exit(1);
    }

    set_edge( i, j, true );

    P->setrow( j, 0.0 );
    size_t m = 1 << (2*num_parents(j));
    kmer K;
    for( K = 0; K < m; K++ ) {
        P->set( j, K, pseudocount );
    }


    std::deque<sequence*>::const_iterator seq;
    for( seq = data->begin(); seq != data->end(); seq++ ) {
        P->inc( j, (*seq)->get( parents[j*n] ) );
    }

    P->dist_normalize_row( j );
    P->dist_conditionalize_row( j );
}


double motif_log_likelihood( motif& M0, motif& M1 )
{
    double L, L0, L1;

    L = 0.0;
    
    std::deque<sequence*>::const_iterator i;
    for( i = M0.data->begin(); i != M0.data->end(); i++ ) {
        L0 = M0.eval( **i );
        L1 = M1.eval( **i );
        L += L0 - log( exp(L0) + exp(L1) );
    }

    for( i = M1.data->begin(); i != M1.data->end(); i++ ) {
        L0 = M0.eval( **i );
        L1 = M1.eval( **i );
        L += L1 - log( exp(L0) + exp(L1) );
    }

    return L;
}


void train_motifs( motif& M0, motif& M1 )
{
    log_puts( LOG_MSG, "training motifs ...\n" );
    log_indent();

    if( M0.n != M1.n ) {
        log_printf( LOG_ERROR, "Motif models of mismatching size. (%zu != %zu)\n", M0.n, M1.n );
        exit(1);
    }

    size_t i, j;

    double aic, aic_curr, aic_best;
    size_t j_best, i_best;

    //double n_obs    = M0.data->size() + M1.data->size();
    double n_params = M0.num_params() + M1.num_params();

    double L; /* log likelihood */

    L = motif_log_likelihood( M0, M1 );
    aic_curr = L - n_params;

    log_printf( LOG_MSG, "base aic = L - n = %0.4e - %f = %0.4e\n", L, n_params, aic_curr );

    while( true ) {
        aic_best = GSL_NEGINF;
        j_best = i_best = 0;

        log_puts( LOG_MSG, "trying edges ... \n" );
        log_indent();

        for( j = 0; j < M0.n; j++ ) {


            /* position j is not currently included in the model */
            if( !M0.has_edge( j, j ) ) {
                M0.store_row(j);
                M1.store_row(j);
                M0.add_edge( j, j );
                M1.add_edge( j, j );

                L        = motif_log_likelihood( M0, M1 );
                n_params = M0.num_params() + M1.num_params();
                aic      = L - n_params;
                log_printf( LOG_MSG, "edge (%zu, %zu): aic = L - n = %0.4e - %f = %0.4e\n",
                                      j, j, L, n_params, aic );


                if( aic > aic_best ) {
                    aic_best = aic;
                    i_best = j;
                    j_best = j;
                }
        
                M0.set_edge( j, j, false );
                M1.set_edge( j, j, false );
                M0.restore_stored_row();
                M1.restore_stored_row();

            }
            /* position j is included in the model */
            else {
                for( i = 0; i < j; i++ ) {
                    if( M0.has_edge( i, j ) ) {
                        log_printf( LOG_MSG, "edge (%zu, %zu): edge exists\n", i, j );
                        continue;
                    }

                    if( M0.num_parents(j) >= M0.k ) {
                        log_printf( LOG_MSG, "edge (%zu, %zu): %zu already has %zu parents\n",
                                    i, j, j, M0.num_parents(j) );
                        continue;
                    }

                    M0.store_row(j);
                    M1.store_row(j);
                    M0.add_edge( i, j );
                    M1.add_edge( i, j );

                    L        = motif_log_likelihood( M0, M1 );
                    n_params = M0.num_params() + M1.num_params();
                    aic      = L - n_params;
                    log_printf( LOG_MSG, "edge (%zu, %zu): aic = L - n = %0.4e - %f = %0.4e\n",
                                          i, j, L, n_params, aic );

                    if( aic > aic_best ) {
                        aic_best = aic;
                        i_best = i;
                        j_best = j;
                    }
            
                    M0.set_edge( i, j, false );
                    M1.set_edge( i, j, false );
                    M0.restore_stored_row();
                    M1.restore_stored_row();
                }
            }
        }

        log_printf( LOG_MSG, "best edge: (%zu,%zu) with aic = %0.4e\n",
                             i_best, j_best, aic_best );

        if( aic_best <= aic_curr || aic_best == GSL_NEGINF ) break;

        aic_curr = aic_best;
        M0.add_edge( i_best, j_best );
        M1.add_edge( i_best, j_best );

        log_unindent();
    }

    log_unindent();
    log_puts( LOG_MSG, "done.\n" );
}




