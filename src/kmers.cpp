
#include "kmers.hpp"
#include "common.hpp"
#include "logger.h"

#include <string>
#include <cstring>
#include <gsl/gsl_math.h>

kmer_matrix::kmer_matrix( const YAML::Node& node )
{
    size_t size1, size2;
    node["k"] >> k;

    node["size1"] >> size1;
    node["size2"] >> size2;

    A = gsl_matrix_alloc( size1, size2 );

    const YAML::Node& node_A = node["A"];
    size_t i;
    for( i = 0; i < size1*size2; i++ ) {
        node_A[i] >> A->data[i];
    }

    stored_row = gsl_vector_alloc( size2 );
}

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


void kmer_matrix::to_yaml( YAML::Emitter& out ) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "k";
    out << YAML::Value << k;
    out << YAML::Key   << "size1";
    out << YAML::Value << A->size1;
    out << YAML::Key   << "size2";
    out << YAML::Value << A->size2;
    out << YAML::Key   << "A";
    out << YAML::Flow;
    out << YAML::Value;
    out << YAML::BeginSeq;
    size_t i;
    for( i = 0; i < A->size1 * A->size2; i++ ) {
        out << A->data[i];
    }
    out << YAML::EndSeq;

    out << YAML::EndMap;
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


void kmer_matrix::dist_conditionalize( int effective_k )
{
    size_t i;

    for( i = 0; i < A->size1; i++ ) {
        dist_conditionalize_row( i, effective_k );
    }
}


void kmer_matrix::dist_conditionalize_row( size_t i, int effective_k  )
{
    if( effective_k <= 0 ) effective_k = A->size2;

    kmer L;
    kmer L_max = 1 << (2*(effective_k-1));
    kmer K;
    kmer nt;
    double z;

    for( L = 0; L < L_max; L++ ) {
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


void kmer_matrix::log_transform_row( size_t i, int effective_k )
{
    if( effective_k <= 0 ) effective_k = A->size2;
    size_t j_max = 1 << (2*effective_k);

    size_t j;
    for( j = 0; j < j_max; j++ ) {
        set( i, j, log( get( i, j ) ) );
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


sequence::sequence( const char* s, int meta )
    : meta(meta), xs(NULL), n(0)
{
    n = strlen(s);
    if( n > 0 ) {
        xs = new kmer[ n/kmer_max_k + 1 ];
        memset( xs, 0, (n/kmer_max_k + 1)*sizeof(kmer) );

        size_t i, block, offset;
        for( i = 0; i < n; i++ ) {
            block  = i / kmer_max_k;
            offset = i % kmer_max_k;
            xs[block] = xs[block] | (nt2num(s[i]) << (2*offset));
        }
    }
}

sequence::sequence( const sequence& s )
    : xs(NULL), n(0)
{
    meta = s.meta;
    n    = s.n;
    if( n > 0 ) {
        xs = new kmer[ n/kmer_max_k + 1 ];
        memcpy( xs, s.xs, (n/kmer_max_k + 1)*sizeof(kmer) );
    }
}


void sequence::operator=( const sequence& s )
{
    delete[] xs;
    meta = s.meta;
    n    = s.n;
    xs = new kmer[ n/kmer_max_k + 1 ];
    memcpy( xs, s.xs, (n/kmer_max_k + 1)*sizeof(kmer) );
}


sequence::~sequence()
{
    delete[] xs;
}

kmer sequence::get( size_t i ) const
{
    size_t block  = i / kmer_max_k;
    size_t offset = i % kmer_max_k;
    return (xs[block] >> (2*offset)) & 0x3;
}


bool sequence::get( const bool* indexes, size_t maxn, kmer& K, size_t seq_offset ) const
{
    bool nonempty = false;
    size_t block, offset;
    K = 0;
    size_t i;
    for( i = 0; i < maxn; i++ ) {
        if( indexes[i] ) {
            nonempty = true;
            block  = (i+seq_offset) / kmer_max_k;
            offset = (i+seq_offset) % kmer_max_k;
            K = (K<<2) | ((xs[block] >> (2*offset)) & 0x3);
        }
    }

    return nonempty;
}


const double motif::pseudocount = 1;

motif::motif( const YAML::Node& node )
{
    node["n"] >> n;
    node["k"] >> k;
    node["meta"] >> meta;

    parents = new bool[n*n];
    const YAML::Node& node_parents = node["parents"];
    size_t i;
    int parents_i;
    for( i = 0; i < n*n; i++ ) {
        node_parents[i] >> parents_i;
        parents[i] = (bool)parents_i;
    }

    P = new kmer_matrix( node["P"] );
}

motif::motif( size_t n, size_t k, int meta )
    : meta(meta), n(n), k(k)
{
    P = new kmer_matrix( n, k );
    P->setall( 0.0 );

    parents = new bool[n*n];
    memset( parents, 0, n*n*sizeof(bool) );
}


motif::motif( const motif& M )
{
    P = new kmer_matrix( *M.P );
    meta = M.meta;
    n    = M.n;
    k    = M.k;

    parents = new bool[n*n];
    memcpy( parents, M.parents, n*n*sizeof(bool) );
}



motif::~motif()
{
    delete[] parents;
    delete P;
}


void motif::to_yaml( YAML::Emitter& out ) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "n";
    out << YAML::Value << n;

    out << YAML::Key   << "k";
    out << YAML::Value << k;

    out << YAML::Key   << "meta";
    out << YAML::Value << meta;

    out << YAML::Key << "parents";
    out << YAML::Value;
    out << YAML::Flow << YAML::BeginSeq;
    size_t i;
    for( i = 0; i < n*n; i++ ) {
        out << (parents[i] ? 1 : 0);
    }
    out << YAML::EndSeq;

    out << YAML::Key   << "P";
    out << YAML::Value;
    P->to_yaml( out );

    out << YAML::EndMap;
}


double motif::eval( const sequence& seq, size_t offset ) const
{
    double p = 0.0;

    const size_t n = P->n();
    size_t i;
    kmer K;
    for( i = 0; i < n; i++ ) {
        if( !seq.get( parents + i*n, n, K, offset ) ) continue;
        p += P->get( i, K );
    }

    return p;
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
    parents[j*n+i] = x;
}


void motif::store_row( size_t i )
{
    P->store_row( i );
}


void motif::restore_stored_row()
{
    P->restore_stored_row();
}


char* motif::print_model_graph( int offset )
{
    std::string graph_str;
    char* tmp;
    int r;

    graph_str += "digraph {\n";
    graph_str += "splines=\"true\";\n";
    graph_str += "node [shape=\"box\"];\n";

    int i, j;
    for( j = 0; j < (int)n; j++ ) {
        r = asprintf( &tmp, "n%d [label=\"%d\",pos=\"%d,0\",style=\"%s\"];\n",
                            j, j - offset, j*100,
                            parents[j*n+j] ? "solid" : "dotted" 
                            );
        graph_str += tmp;
        free(tmp);
    }


    for( j = 0; j < (int)n; j++ ) {
        if( !parents[j*n+j] ) continue;

        for( i = 0; i < j; i++ ) {
            if( parents[j*n+i] ) {
                r = asprintf( &tmp, "n%d -> n%d;\n", i, j );
                graph_str += tmp;
                free(tmp);
            }

        }

    }

    graph_str += "}\n";
    return strdup(graph_str.c_str());
}


void motif::add_all_edges( const std::deque<sequence*>* data )
{
    size_t i, j;
    size_t n_parents;
    size_t m;
    kmer K;
    std::deque<sequence*>::const_iterator seq;

    for( j = 0; j < n; j++ ) {
        for( i = k-1 <= j ? j-(k-1) : 0; i <= j; i++ ) {
            set_edge( i, j, true );
        }
        P->setrow( j, 0.0 );
        n_parents = num_parents(j);
        m = 1 << (2*n_parents);

        for( K = 0; K < m; K++ ) {
            P->set( j, K, pseudocount );
        }

        for( seq = data->begin(); seq != data->end(); seq++ ) {
            if( (*seq)->meta == meta && (*seq)->get( parents + j*n, n, K ) ) {
                P->inc( j, K );
            }
        }

        P->dist_normalize_row( j );
        P->dist_conditionalize_row( j, n_parents );
        P->log_transform_row( j, n_parents );
    }
}


/* make an edge i --> j
 * That is, condition j on i. */
void motif::add_edge( size_t i, size_t j, const std::deque<sequence*>* data )
{
    if( i > j ) {
        log_printf( LOG_ERROR, "Invalid motif edge (%zu, %zu)\n", i, j );
        exit(1);
    }

    set_edge( i, j, true );

    P->setrow( j, 0.0 );
    size_t n_parents = num_parents(j);
    size_t m = 1 << (2*n_parents);
    kmer K;
    for( K = 0; K < m; K++ ) {
        P->set( j, K, pseudocount );
    }


    std::deque<sequence*>::const_iterator seq;
    for( seq = data->begin(); seq != data->end(); seq++ ) {
        if( (*seq)->meta == meta && (*seq)->get( parents + j*n, n, K ) ) {
            P->inc( j, K );
        }
    }

    P->dist_normalize_row( j );
    P->dist_conditionalize_row( j, n_parents );
    P->log_transform_row( j, n_parents );
}



void motif::remove_edge( size_t i, size_t j, const std::deque<sequence*>* data )
{
    if( i > j ) {
        log_printf( LOG_ERROR, "Invalid motif edge (%zu, %zu)\n", i, j );
        exit(1);
    }

    set_edge( i, j, false );

    P->setrow( j, 0.0 );
    size_t n_parents = num_parents(j);
    size_t m = 1 << (2*n_parents);
    kmer K;
    for( K = 0; K < m; K++ ) {
        P->set( j, K, pseudocount );
    }


    std::deque<sequence*>::const_iterator seq;
    for( seq = data->begin(); seq != data->end(); seq++ ) {
        if( (*seq)->meta == meta && (*seq)->get( parents + j*n, n, K ) ) {
            P->inc( j, K );
        }
    }

    P->dist_normalize_row( j );
    P->dist_conditionalize_row( j, n_parents );
    P->log_transform_row( j, n_parents );
}





double motif_log_likelihood( const motif& M0, const motif& M1,
                             const std::deque<sequence*>* training_seqs )
{
    double L, L0, L1;

    L = 0.0;
    
    std::deque<sequence*>::const_iterator i;
    for( i = training_seqs->begin(); i != training_seqs->end(); i++ ) {
        L0 = M0.eval( **i );
        L1 = M1.eval( **i );
        if( (*i)->meta == 0 ) {
            L += L0 - logaddexp( L0, L1 );
        }
        else if( (*i)->meta == 1 ) {
            L += L1 - logaddexp( L0, L1 );
        }
    }

    return L;
}



double eval_likelihood_vectors( const gsl_vector* l0, const gsl_vector* l1,
                                const int* meta )
{
    if( l0->size != l1->size ) {
        log_printf( LOG_ERROR, "Mismatching likelihood vectors.\n" );
        exit(1);
    }

    double l;
    const size_t n = l0->size;
    size_t i;

    l = 0.0;
    for( i = 0; i < n; i++ ) {
        l += (meta[i] == 0 ?  gsl_vector_get( l0, i ) : gsl_vector_get( l1, i ))
          -  logaddexp( gsl_vector_get( l0, i ), gsl_vector_get( l1, i ) );
    }

    return l;
}




void motif::update_likelihood_column( gsl_matrix* L, size_t j,
                                      const std::deque<sequence*>* training_seqs )
{
    size_t i;
    const size_t n = L->size1;
    const size_t m = L->size2;
    kmer K;

    gsl_vector_set_all( &gsl_matrix_column( L, j ).vector, 0.0 );

    for( i = 0; i < n; i++ ) {
        if( (*training_seqs)[i]->get( parents + j*m, m, K ) ) {
            gsl_matrix_set( L, i, j, P->get( j, K ) );
        }
    }
}





/* various information criterion to try */

/* Akaike Information Criterion */
double aic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return L - c*n_params*n_obs / (n_obs - n_params - 1.0);
}

/* Bayesian (Schwarz) Information Criterion */
double bic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return 2.0*L - c*n_params*log(n_obs);
}


/* Hannan-Quinn Information Criterion */
double hqic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return L - c*n_params*log(log(n_obs));
}



void train_motifs( motif& M0, motif& M1,
                   const std::deque<sequence*>* training_seqs,
                   size_t max_dep_dist, double complexity_penalty )
{

    log_puts( LOG_MSG, "training motifs (forwards) ...\n" );
    log_indent();

    if( M0.n != M1.n ) {
        log_printf( LOG_ERROR, "Motif models of mismatching size. (%zu != %zu)\n", M0.n, M1.n );
        exit(1);
    }

    double (*compute_ic)( double, double, double, double ) = aic;


    size_t i, j;
    size_t i_start;

    /* likelihood matrices: 
     * Lk_ij gives the likelihood of training example i on node j in model k.
     */
    gsl_matrix* L0 = gsl_matrix_calloc( training_seqs->size(), M0.n );
    gsl_matrix* L1 = gsl_matrix_calloc( training_seqs->size(), M0.n );

    /* likelihood vectors:
     * summed columns of L0, L1, maintained to minimize redundant computation.
     * */
    gsl_vector* l0 = gsl_vector_calloc( training_seqs->size() );
    gsl_vector* l1 = gsl_vector_calloc( training_seqs->size() );


    /* 0-1 vectors giving labeling each sequence as foreground or background */
    int* meta = new int[ training_seqs->size() ];

    for( i = 0; i < training_seqs->size(); i++ ) {
        meta[i] = (*training_seqs)[i]->meta;
    }


    /* backup likelihood matrix columns, to restore state after trying a new edge */
    gsl_vector* c0 = gsl_vector_alloc( training_seqs->size() );
    gsl_vector* c1 = gsl_vector_alloc( training_seqs->size() );
    


    /* keeping track of the optimal edge */
    double ic, ic_curr, ic_best;
    size_t j_best, i_best;

    /* parameters to compute information criterion */
    double n_obs    = training_seqs->size();
    double n_params = M0.num_params() + M1.num_params();


    /* log posterier probability */
    double l; 


    /* baseline ic */
    l = eval_likelihood_vectors( l0, l1, meta );
    ic_curr = compute_ic( l, n_obs, n_params, complexity_penalty );

    //log_printf( LOG_MSG, "l = %e, k = %0.0e, ic = %0.4e\n", l, n_params, ic_curr );


    size_t round = 0;

    while( true ) {
        round++;

        log_printf( LOG_MSG, "round %4d (ic = %0.4e) ", round, ic_curr );

        ic_best = GSL_NEGINF;
        j_best = i_best = 0;

        for( j = 0; j < M0.n; j++ ) {

            if( M0.has_edge( j, j ) ) {
                if( max_dep_dist == 0 || j <= max_dep_dist ) i_start = 0;
                else i_start = i - max_dep_dist;
            }
            else i_start = j;


            for( i = i_start; i <= j; i++ ) {
                log_puts( LOG_MSG, "." );

                if( M0.has_edge( i, j ) ) {
                    continue;
                }

                if( M0.num_parents(j) >= M0.k ) {
                    continue;
                }

                /* keep track of the old parameters to avoid retraining */
                M0.store_row(j);
                M1.store_row(j);

                /* keep track of the old likelihoods to avoid reevaluating */
                gsl_matrix_get_col( c0, L0, j );
                gsl_matrix_get_col( c1, L1, j );

                /* add edge */
                M0.add_edge( i, j, training_seqs );
                M1.add_edge( i, j, training_seqs );

                /* evaluate likelihoods for that column */
                M0.update_likelihood_column( L0, j, training_seqs );
                M1.update_likelihood_column( L1, j, training_seqs );

                /* update training example likelihoods */
                gsl_vector_sub( l0, c0 );
                gsl_vector_add( l0, &gsl_matrix_column( L0, j ).vector );

                gsl_vector_sub( l1, c1 );
                gsl_vector_add( l1, &gsl_matrix_column( L1, j ).vector );


                l        = eval_likelihood_vectors( l0, l1, meta );
                n_params = M0.num_params() + M1.num_params();
                ic       = compute_ic( l, n_obs, n_params, complexity_penalty );

                if( ic > ic_best ) {
                    ic_best = ic;
                    i_best = i;
                    j_best = j;
                }
        

                /* remove edge */
                M0.set_edge( i, j, false );
                M1.set_edge( i, j, false );

                /* restore previous parameters */
                M0.restore_stored_row();
                M1.restore_stored_row();

                /* restore previous likelihoods */
                gsl_vector_sub( l0, &gsl_matrix_column( L0, j ).vector );
                gsl_vector_add( l0, c0 );

                gsl_vector_sub( l1, &gsl_matrix_column( L1, j ).vector );
                gsl_vector_add( l1, c1 );

                gsl_matrix_set_col( L0, j, c0 );
                gsl_matrix_set_col( L1, j, c1 );
            }
        }

        log_puts( LOG_MSG, "\n" );

        if( ic_best <= ic_curr || ic_best == GSL_NEGINF ) break;

        ic_curr = ic_best;

        M0.add_edge( i_best, j_best, training_seqs );
        M1.add_edge( i_best, j_best, training_seqs );

        gsl_vector_sub( l0, &gsl_matrix_column( L0, j_best ).vector );
        gsl_vector_sub( l1, &gsl_matrix_column( L1, j_best ).vector );

        M0.update_likelihood_column( L0, j_best, training_seqs );
        M1.update_likelihood_column( L1, j_best, training_seqs );

        gsl_vector_add( l0, &gsl_matrix_column( L0, j_best ).vector );
        gsl_vector_add( l1, &gsl_matrix_column( L1, j_best ).vector );

    }



    delete[] meta;
    gsl_matrix_free( L0 );
    gsl_matrix_free( L1 );
    gsl_vector_free( l0 );
    gsl_vector_free( l1 );
    gsl_vector_free( c0 );
    gsl_vector_free( c1 );


    log_unindent();
}




void train_motifs_backwards( motif& M0, motif& M1,
                             const std::deque<sequence*>* training_seqs,
                             size_t max_dep_dist, double complexity_penalty )
{
    log_puts( LOG_MSG, "training motifs (backwards) ...\n" );
    log_indent();


    if( M0.n != M1.n ) {
        log_printf( LOG_ERROR, "Motif models of mismatching size. (%zu != %zu)\n", M0.n, M1.n );
        exit(1);
    }

    double (*compute_ic)( double, double, double, double ) = aic;

    size_t i, j;
    size_t i_last;

    /* likelihood matrices: 
     * Lk_ij gives the likelihood of training example i on node j in model k.
     */
    gsl_matrix* L0 = gsl_matrix_calloc( training_seqs->size(), M0.n );
    gsl_matrix* L1 = gsl_matrix_calloc( training_seqs->size(), M0.n );

    /* likelihood vectors:
     * summed columns of L0, L1, maintained to minimize redundant computation.
     * */
    gsl_vector* l0 = gsl_vector_calloc( training_seqs->size() );
    gsl_vector* l1 = gsl_vector_calloc( training_seqs->size() );


    /* 0-1 vectors giving labeling each sequence as foreground or background */
    int* meta = new int[ training_seqs->size() ];

    for( i = 0; i < training_seqs->size(); i++ ) {
        meta[i] = (*training_seqs)[i]->meta;
    }


    /* backup likelihood matrix columns, to restore state after trying a new edge */
    gsl_vector* c0 = gsl_vector_alloc( training_seqs->size() );
    gsl_vector* c1 = gsl_vector_alloc( training_seqs->size() );
    


    /* keeping track of the optimal edge */
    double ic, ic_curr, ic_best;
    size_t j_best, i_best;


    /* start with all edges */
    M0.add_all_edges( training_seqs );
    M1.add_all_edges( training_seqs );


    /* parameters to compute information criterion */
    double n_obs    = training_seqs->size();
    double n_params = M0.num_params() + M1.num_params();


    /* log posterier probability */
    double l; 

    /* baseline ic */
    l = eval_likelihood_vectors( l0, l1, meta );
    ic_curr = compute_ic( l, n_obs, n_params, complexity_penalty );

    size_t round = 0;

    while( true ) {
        round++;

        log_printf( LOG_MSG, "round %4d (ic = %0.4e) ", round, ic_curr );

        ic_best = GSL_NEGINF;
        j_best = i_best = 0;

        for( j = 0; j < M0.n; j++ ) {

            /* do not remove the position unless it has no edges */
            i_last = M0.num_parents(j) > 1 ? j-1 : j;

            for( i = 0; i <= i_last; i++ ) {
                log_puts( LOG_MSG, "." );

                if( !M0.has_edge( i, j ) ) continue;

                /* keep track of the old parameters to avoid retraining */
                M0.store_row(j);
                M1.store_row(j);

                /* keep track of the old likelihoods to avoid reevaluating */
                gsl_matrix_get_col( c0, L0, j );
                gsl_matrix_get_col( c1, L1, j );

                /* remove edge */
                M0.remove_edge( i, j, training_seqs );
                M1.remove_edge( i, j, training_seqs );

                /* evaluate likelihoods for that column */
                M0.update_likelihood_column( L0, j, training_seqs );
                M1.update_likelihood_column( L1, j, training_seqs );

                /* update training example likelihoods */
                gsl_vector_sub( l0, c0 );
                gsl_vector_add( l0, &gsl_matrix_column( L0, j ).vector );

                gsl_vector_sub( l1, c1 );
                gsl_vector_add( l1, &gsl_matrix_column( L1, j ).vector );

                /* evaluate likelihood / ic */
                l        = eval_likelihood_vectors( l0, l1, meta );
                n_params = M0.num_params() + M1.num_params();
                ic       = compute_ic( l, n_obs, n_params, complexity_penalty );

                if( ic > ic_best ) {
                    ic_best = ic;
                    i_best = i;
                    j_best = j;
                }

                /* replace edge */
                M0.set_edge( i, j, true );
                M1.set_edge( i, j, true );

                /* restore previous parameters */
                M0.restore_stored_row();
                M1.restore_stored_row();

                /* restore previous likelihoods */
                gsl_vector_sub( l0, &gsl_matrix_column( L0, j ).vector );
                gsl_vector_add( l0, c0 );

                gsl_vector_sub( l1, &gsl_matrix_column( L1, j ).vector );
                gsl_vector_add( l1, c1 );

                gsl_matrix_set_col( L0, j, c0 );
                gsl_matrix_set_col( L1, j, c1 );
            }
        }

        log_puts( LOG_MSG, "\n" );


        if( ic_best <= ic_curr || ic_best == GSL_NEGINF ) break;

        ic_curr = ic_best;

        M0.remove_edge( i_best, j_best, training_seqs );
        M1.remove_edge( i_best, j_best, training_seqs );

        gsl_vector_sub( l0, &gsl_matrix_column( L0, j_best ).vector );
        gsl_vector_sub( l1, &gsl_matrix_column( L1, j_best ).vector );

        M0.update_likelihood_column( L0, j_best, training_seqs );
        M1.update_likelihood_column( L1, j_best, training_seqs );

        gsl_vector_add( l0, &gsl_matrix_column( L0, j_best ).vector );
        gsl_vector_add( l1, &gsl_matrix_column( L1, j_best ).vector );
    }

    

    delete[] meta;
    gsl_matrix_free( L0 );
    gsl_matrix_free( L1 );
    gsl_vector_free( l0 );
    gsl_vector_free( l1 );
    gsl_vector_free( c0 );
    gsl_vector_free( c1 );


    log_unindent();
}


