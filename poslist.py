#! /usr/bin/env python
import sys, os;
sys.path.append('/bio/lib/python2.7/site-packages');
sys.path.append('/bio/lib/site-python');
import numpy as np;

def pos(genomefile, hicfile, outfile, output = False):        
    try:
        savedhic = np.load(hicfile);
        savedgen = np.load(genomefile);
        hic_set = set(list(savedhic['hic']));
        gen_set = set(list(savedgen['genome']));
    except IOError:
        print 'file not found';
        raise;
    except KeyError:
        print 'file exists, but parameter not found';
        raise;
    else: 
        pos = list(hic_set.intersection(gen_set));
        pos.sort();
        np.savez(outfile, pos = np.array(pos));
        if(output):
            print 'Hi-C data contains {0} positions'.format(len(hic_set));
            print 'Genome data contains {0} positions'.format(len(gen_set));
            print 'There are {0} common positions'.format(len(pos));
                 
def main(argv):
    try:
        genomefile = argv[1];
        hicfile = argv[2];
        outfile = argv[3];                          
    except IndexError:
        print 'usage: {p} <genome pos> <hic pos> <out pos>'.format(p = argv[0]);

        hicfile = 'hic.pos.npz';
        genomefile = 'genome.pos.npz';
        outfile = 'pos.npz';
        print 'example: {p} {g} {h} {o}'.format(p = argv[0],
                                                g = genomefile,
                                                h = hicfile,
                                                o = outfile);
        exit(1);
    else:
        pos(genomefile, hicfile, outfile, output = True);
        return;

if __name__ == "__main__":
    main(sys.argv);
