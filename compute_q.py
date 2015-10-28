#! /usr/bin/env python
import sys, os;
import itertools;
from math import isnan;
import jellyfish;
sys.path.append('/bio/lib/python2.7/site-packages');
sys.path.append('/bio/lib/site-python');
import numpy as np;


def allkmers(k):
    alph = ('A', 'C', 'G', 'T');
    return itertools.product(alph, repeat = k);

def kmercount(k, pos, chr = 21,
              fname_head = '/data/yt/GRCh38.p2.ch21/GRCh38.p2'):
    try:
        fname = '{head}.ch{chr}.{pos}.fasta.{k}.jf'.format(head = fname_head,
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
        # print c.T
        # print len(c);
        return c;

def compute_q(k, binsize, pos_set, hic_file, outfile):
    with open(hic_file, 'r') as f:
        q = np.zeros(((1 << (2 * k)), (1 << (2 * k))), dtype = np.float64);
        scale = 0.5 / ((binsize - k + 1) * (binsize - k + 1));

        datapoints = 0;
        linenum = 0;

        for l in f:
            linenum += 1;
            (si, sj, smij) = l[:-1].split();
            [i, j, mij] = [int(si), int(sj), float(smij)];
            if((i in pos_set) and (j in pos_set) and (not isnan(mij))):
                datapoints += 1;

                ci = kmercount(k, i);
                cj = kmercount(k, j);
                #print 'min(ci.dot(cj.T)) = {0}'.format(ci.dot(cj.T));
                #print 'min(mij * scale * ci.dot(cj.T)) = {0}'.format(mij * scale * ci.dot(cj.T));
                #print 'i = {i}\tj = {j}\tmij = {mij}'.format(i = i, j = j, mij = mij);
                #print (mij * scale * ci.dot(cj.T)).min();

                q += (mij * scale * ci.dot(cj.T));

                print q;
                print '{space}# of data points: {i}, # of proceeded lines: {l}'.format(space = (' ' * 20),
                                                                                       i = datapoints,
                                                                                       l = linenum);
                #print 'min(q) = {0}'.format(q.min());
                #print '{0} {1} {2}'.format(i, j, mij);
        print q;
        np.savez(outfile, ccsum = q);
    return;

def main(argv):
    try:
        datadir = argv[1];
        chr = argv[2];
        res = argv[3];
        k = int(argv[4]);
        bin = int(argv[5]);
        norm_str = argv[6];
        exp_str = argv[7];
        posfile = argv[8];
        outfile = argv[9];
    except IndexError:
        print 'usage: {p} <data dir> <chr> <res> <k> <bin> <norm> <exp> <pos_file> <out_file>'.format(p = argv[0]);
        print '  data dir  : Hi-C data directory';
        print '  chr  : chromosome';
        print '  res  : resolution';
        print '  k    : length of k-mer';
        print '  bin  : size of bins';
        print '  norm : normalize method';
        print '  exp  : expected value method';
        print '  pos_file  : position list file';
        print '  out_file  : output file';
        print '-----'
        datadir='/data/yt/GM12878_combined/1kb_resolution_intrachromosomal/chr21/MAPQGE30'
        chr = 'chr21';
        res = '1kb';
        k = 3;
        bin = 1000;
        norm_str = 'KR';
        exp_str = 'KR';        
        posfile = './pos.npz'
        outfile = '/work2/yt/chromcode/ccsum.npz'
        print 'example: {prog} {d} {c} {r} {k} {b} {n} {e} {p} {o}'.format(prog = argv[0],
                                                                           d = datadir,
                                                                           c = chr,
                                                                           r = res,
                                                                           k = k,
                                                                           b = bin,
                                                                           n = norm_str,
                                                                           e = exp_str,
                                                                           p = posfile,
                                                                           o = outfile);
    else:
        hic_file = ''.join([datadir, '/', chr, '_', res, '.', norm_str, '.', exp_str, '.dat']);
        try:
            pos_set = set(list(np.load(posfile)['pos']));
        except KeyError:
            print 'wrong input of posfile: need to have a variable called pos';
            raise;
        else:
            compute_q(k, bin, pos_set, hic_file, outfile);
            # kmercount(5, 5010000);
            return;

if __name__ == "__main__":
    main(sys.argv);
