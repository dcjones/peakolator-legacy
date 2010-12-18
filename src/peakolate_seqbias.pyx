
from sequencing_bias cimport *
import argparse

include 'sequencing_bias.pxi'

desc = \
'''
Train a statistical model for technical bias in sequencing data, given a
reference genome and reads in BAM format.
'''

def main():
    ap = argparse.ArgumentParser( description = desc, prog = 'peakolate-seqbias' )

    ap.add_argument( 'ref_fn',   metavar = 'reference.fa',
                     help = 'genome sequence against which reads were mapped' )

    ap.add_argument( 'reads_fn', metavar = 'reads.bam',
                     help = 'BAM file containing aligned reads'  )

    ap.add_argument( '-o', '--output', default='seqbias.yml', metavar='out.yml',
                     help = 'output the model to the given file '
                            '(default: "seqbias.yml")' )

    ap.add_argument( '-d', '--dot', default='seqbias.dot', metavar = 'out.dot',
                     help = 'output graph of the model in dot format '
                            '(default: "seqbias.dot")' )

    ap.add_argument( '-L', default = 20, type=int,
                     help = 'consider at most L nucleotides to the left of the '
                             'read start (default:20)' )

    ap.add_argument( '-R', default = 20, type=int,
                     help = 'consider at most R nucleotides to the right of the '
                             'read start (default: 20)' )

    ap.add_argument( '-n', default = 100000, type=int,
                     help = 'consider the first n reads when training '
                            '(default: 100000)' )

    ap.add_argument( '--complexity-penalty', metavar = 'C', type=float, default = 1.0,
                     help = 'reweight the information criterion to '
                            'increase/decrease the model complexity penalty '
                            '(default: 1.0, i.e. no weighting)')

    ap.add_argument( '--backward', default=False, action='store_true',
                     help = 'train using backward stepwise regression\n'
                            '(by default, forward regression is done)' )

    args = ap.parse_args()

    seqbias = sequencing_bias()
    seqbias.train( args.ref_fn, args.reads_fn, args.n, args.L, args.R,
                   args.backward, args.complexity_penalty )

    seqbias.save( args.output )

    dot_str = seqbias.print_model_graph()
    f = open( args.dot, 'w' )
    f.write( dot_str )
    f.close()


if __name__ == '__main__': main()


