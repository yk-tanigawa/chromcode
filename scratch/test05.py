#! /usr/bin/env python

import os, sys, shlex, subprocess;
import numpy as np;
from collections import deque;
from datetime import datetime;
import itertools;

def test05(k):
    alph = ('A', 'C', 'G', 'T');
    for km in itertools.product(alph, repeat=k):
        kmer = ''.join(km);
        print kmer;

def test01():
    print subprocess.check_output('ls');

def test03():
    print 1;
    retcode = subprocess.check_call(["ls", "-l"], stdin = None);
    #retcode = subprocess.check_call(["sleep", "2"]);
    print retcode;
    return;

def test04():
    cmd1 = "ls -lt "
    cmd2 = "head -n 5 "
    p1 = subprocess.Popen(cmd1.strip().split(" "), stdout=subprocess.PIPE);
    p2 = subprocess.Popen(cmd2.strip().split(" "), stdin=p1.stdout, stdout=subprocess.PIPE);
    p1.stdout.close(); # Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate();
    print output[0];
    


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

def mkjfdb(fnames, fpath):
    for f in fnames:
        jfcmd = "jellyfish count -m 5 -s 100M -t 6 {input} -o {output}".format(input = f,
                                                                               output = fpath + get_fname_head(f) + '.jf');
        p = subprocess.Popen(jfcmd.strip().split(" "), stdout=subprocess.PIPE);
        output = p.communicate();
        if(output[1] != None):
            print output[1];
        #print jfcmd;


def get_fname_head(fname):
    if(fname.endswith('.fasta')):
        return fname.split('/')[-1][:-6];
    elif(fname.endswith('.fa')):
        return fname.split('/')[-1][:-3];
    else:
        return fname.split('/')[-1];
           

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
        fname_head = get_fname_head(fname);
        fpath = '/'.join(fname.split('/')[:-1]) + '/' + fname_head + '/';
        #print fpath;

        # read a fasta file
        (head, seq) = get_ref(fp);

        if(not os.path.exists(fpath)):
            os.makedirs(fpath);
        fnames = [];
        mkfasta(head, seq, 1000, fpath + fname_head, fnames);
        fp.close;

        mkjfdb(fnames, fpath);
        print '{0} FASTA files created.'.format(len(fnames));
        #for fn in fnames:
        #    print fn;

    return;

if __name__ == "__main__":
    test05(2);
    #main(sys.argv);

