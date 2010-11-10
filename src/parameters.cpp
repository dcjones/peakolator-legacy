
#include "parameters.hpp"
#include "emppval.hpp"


parameters::parameters()
    : alpha(0.01)
    , r(0.0)
    , p(1.0)
    , d_min(70)
    , d_max(100000)
    , n_mc(100)
    , padj_spacing(10000)
    , padj_n(20)
    , padj(NULL)
{
    mpfr_class::set_dprec(_mpfr_prec_);
}


parameters::parameters( const parameters& p )
    : alpha(p.alpha)
    , r(p.r)
    , p(p.p)
    , d_min(p.d_min)
    , d_max(p.d_max)
    , n_mc(p.n_mc)
    , padj_spacing(p.padj_spacing)
    , padj_n(p.padj_n)
{
    padj = p.padj ? new emppval(*p.padj) : NULL;
    mpfr_class::set_dprec(_mpfr_prec_);
}

parameters::~parameters()
{
    delete padj;
}

parameters* parameters::copy() const
{
    return new parameters( *this );
}


void parameters::rebuild_lookup( int m, int n )
{
    dist.build( r, p, m, n );
}

void parameters::build_padj()
{
    if( padj ) delete padj;
    padj = new emppval( this );
}

