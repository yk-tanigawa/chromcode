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
        fname = '{head}.ch{chr}.{pos}.jf'.format(head = fname_head,
                                                 chr = chr,
                                                 pos = pos);
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

def main(argv):
    kmercount(5, 5010000);
    return;

if __name__ == "__main__":
    main(sys.argv);
