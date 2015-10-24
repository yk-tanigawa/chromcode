#! /usr/bin/env python

import os;
import sys;
import subprocess;
import numpy as np;
from collections import deque;
from datetime import datetime;
import itertools;

def test01():
    print subprocess.check_output('ls');

def get_ref(fp):
    # sequence header
    line = fp.readline();
    head = line[0 : len(line) - 1];

    # sequence body
    seq_buf = [];
    for line in fp:
        if(line.startswith('>')):
            return(head, seq);
        else:
            seq = line[0 : len(line) - 1];
            seq_buf.append(seq);
    seq = ''.join(seq_buf);
    return (head, seq);

def mkfasta_sub(head, seq, fname):
    fp = open(fname, 'w');

    fp.write(head);
    fp.write(seq);
    fp.write('\n');
    
    fp.close();
    return;

def mkfasta(head, seq, interval, fname_head, fnames):
    for i in xrange(len(seq) / interval + 1):
        subseq = seq[i * interval : (i + 1) * interval];    
        if('N' not in subseq and len(subseq) == interval):
            posstr = '{0}'.format(i * interval);
            subseq_head =  head + ' ' + posstr + '\n';
            subseq_fname = fname_head + '.' + posstr + '.fasta';
            #print subseq_head;    
            #print subseq;
            #print subseq_fname;
            mkfasta_sub(subseq_head, subseq, subseq_fname);
            fnames.append(subseq_fname);
    return;

def main(argv):
    try:
        # try to open fasta file
        fname = argv[1];
        fp = open(fname, 'r')
    except IndexError:
        print 'Usage: {prog_name} <FASTA FILE>'.format(prog_name = argv[0]);
    except IOError:
        print '{file} cannot be opened.'.format(file = fname);
    else:
        if(fname.endswith('.fasta')):
            fname_head = fname.split('/')[-1][:-6];
        elif(fname.endswith('.fa')):
            fname_head = fname.split('/')[-1][:-3];
        else:
            fname_head = fname.split('/')[-1];
        fpath = '/'.join(fname.split('/')[:-1]) + '/' + fname_head + '/';
        #print fpath;

        # read a fasta file
        (head, seq) = get_ref(fp);


        os.makedirs(fpath);
        fnames = [];
        mkfasta(head, seq, 1000, fpath + fname_head, fnames);

        fp.close;

        print '{0} FASTA files created.'.format(len(fnames));
        #for fn in fnames:
        #    print fn;

    return;

if __name__ == "__main__":
    main(sys.argv);
#    if(len(sys.argv) < 4):
#        sys.stderr.write("Some arguments are missing\n");
#        sys.exit(1);
#    else:
#        main(sys.argv);
