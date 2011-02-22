
#include "parameters.hpp"
#include "emppval.hpp"


parameters::parameters()
    : alpha(0.01)
    , r(0.0)
    , p(1.0)
    , a(0.0)
    , d_min(70)
    , d_max(100000)
    , padj(NULL)
{
}


parameters::parameters( const parameters& p )
    : alpha(p.alpha)
    , r(p.r)
    , p(p.p)
    , a(p.a)
    , d_min(p.d_min)
    , d_max(p.d_max)
{
    padj = p.padj ? new emppval(*p.padj) : NULL;
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
    dist.build( r, p, a, m, n );
}

void parameters::build_padj()
{
    if( padj ) delete padj;
    padj = new emppval( this );
}


