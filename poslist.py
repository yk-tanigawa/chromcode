#! /usr/bin/env python
import sys, os;
import numpy as np;

def pos(hicfile = 'hic.pos.npz', 
        genomefile = 'genome.pos.npz',
        outfile = 'pos.npz', 
        output = False):
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
        pos = np.array(list(hic_set.intersection(gen_set)));
        np.savez(outfile, pos = pos);
        if(output):
            print 'Hi-C data contains {0} positions'.format(len(hic_set));
            print 'Genome data contains {0} positions'.format(len(gen_set));
            print 'There are {0} common positions'.format(len(pos));
                 
def main(argv):
    pos(output = True);
    return;

if __name__ == "__main__":
    main(sys.argv);
