#! /usr/bin/env python
import sys, os;
import itertools;
import numpy as np;
import jellyfish;

def allkmers(k):
    alph = ('A', 'C', 'G', 'T');
    return itertools.product(alph, repeat = k);

def kmercount(k, pos, chr = 21,
              fname_head = '../../data/GRCh38.p2.ch21/GRCh38.p2'):
    try:
        fname = '{head}.ch{chr}.{pos}.{k}.jf'.format(head = fname_head,
                                                     chr = chr,
                                                     pos = pos, 
                                                     k = k);
        qf = jellyfish.QueryMerFile(fname);
    except RuntimeError:
        raise;
    else:        
        # initialize with pseudo count
        # add 0.5 for smoothing
        # store data in doble quantity to use int vector
        c = np.ones((1 << (2 * k), 1), dtype = np.uint16);
        i = 0;
        for l in allkmers(k):
            c[i][0] += 2 * qf[jellyfish.MerDNA(''.join(l))];
            i += 1;
        # print c;
        print c.T
        print len(c);
        # print c.dot(c.T);
        return c;

def compute_q(k, pos_set, hic_file):
    with open(hic_file, 'r') as f:
        q = np.zeros(((1 << (2 * k)), (1 << (2 * k))), dtype = np.int64);
        for l in f:
            (si, sj, smij) = l[:-1].split();
            [i, j, mij] = [int(si), int(sj), float(smij)];
            if((i in pos_set) and (j in pos_set)):
                ci = kmercount(k, i);
                cj = kmercount(k, j);
                q + ci.dot(cj.T);
                #print '{0} {1} {2}'.format(i, j, mij);
        print q;
    return;

def main(argv):
    datadir='../../data/GM12878_conbined_1kb_intra_chr21_MAPQGE30';
    chr = 'chr21';
    res = '1kb';
    norm_str = 'KR';
    exp_str = 'KR';
    hic_file = ''.join([datadir, '/', chr, '_', res, '.', norm_str, '.', exp_str, '.dat']);
    pos_set = set(list(np.load('../pos.npz')['pos']));

    compute_q(5, pos_set, hic_file);

    kmercount(5, 5010000);
    return;

if __name__ == "__main__":
    main(sys.argv);
