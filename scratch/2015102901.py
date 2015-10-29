#! /usr/bin/env python
import sys, os;
import itertools;
sys.path.append('/bio/lib/python2.7/site-packages');
sys.path.append('/bio/lib/site-python');
import numpy as np;


def allkmers(k):
    alph = ('A', 'C', 'G', 'T');
    return itertools.product(alph, repeat = k);

def kmercount(k, fname):
    try:
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
        # print c.T
        # print len(c);
        return c;


def main(argv):
    try:
        jf = argv[1];
        k = argv[2];
    except IndexError:
        print 'usage: {prog} <Jellyfish DB file> <k>'.format(argv[0]);
        exit 1;
    else:
        c = kmercount(k, jf);
        print c;
        np.savetxt(jf + '.csv', c);
    return;

if __name__ == "__main__":
    main(sys.argv);
